"""
This script generates a synthetic graph based on an input edge list and clustering information.
It uses the graph-tool library to create a stochastic block model (SBM) graph.

PARALLELIZED VERSION - Performance improvements:
- Parallel processing for key computation-intensive tasks
- Optimized data structures and algorithms
- Efficient memory management
- Performance monitoring and timing (with verbose flag)

Usage:
    python gen_SBM_parallel.py -f <edge_list_filepath> -c <clustering_filepath> -o <output_directory> [-j <num_jobs>] [-v]

Arguments:
    -f, --filepath: Path to the input edge list file.
    -c, --cluster_filepath: Path to the clustering file.
    -o, --output_directory: Directory to save the generated graph.
    -j, --jobs: Number of parallel jobs (default: number of CPU cores).
    -v, --verbose: Enable verbose output with performance statistics and progress information.

Example:
    python gen_SBM_parallel.py -f edge_list.txt -c clustering.txt -o output_directory -j 8 -v

Original source: https://github.com/illinois-or-research-analytics/lanne2_networks/blob/main/generate_synthetic_networks/gen_SBM.py
Original author: Lahari Anne (@lanne2)
"""

import pandas as pd
import numpy as np
import graph_tool.all as gt  # type: ignore[import]
import typer
import os
import networkit as nk
from scipy.sparse import lil_matrix, csr_matrix  # type: ignore[import]
import time
from typing import Dict, List, Set, Tuple
from pathlib import Path
import threading
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm  # For progress bars
import psutil  # type: ignore[import]
import queue


def create_logger(verbose: bool):
    """
    Create a logger function that only prints when verbose mode is enabled.
    
    Args:
        verbose: Whether to enable verbose output
        
    Returns:
        A logger function that conditionally prints
    """
    # Create a thread-safe print function
    print_lock = threading.Lock()
    
    def log(message: str, always: bool = False, return_verbose: bool = False):
        if return_verbose:
            return verbose
        if verbose or always:
            with print_lock:
                print(message)
    
    return log


def monitor_resources(logger):
    """
    Monitor and log current memory usage.
    
    Args:
        logger: The logger function to use
    """
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    logger(f"Memory usage: {memory_info.rss / (1024 * 1024):.2f} MB")


def read_graph(filepath: str, logger) -> Tuple[nk.Graph, Dict]:
    """
    Reads the graph from the given edge list file and returns a Networkit graph and a node mapping.

    Args:
        filepath: Path to the edge list file
        logger: Logger function for output
        
    Returns: 
        nk.Graph: Networkit graph object
        dict: Node mapping from string to integer
    """
    start_time = time.time()
    edgelist_reader = nk.graphio.EdgeListReader("\t", 0, directed=False, continuous=False)
    nk_graph = edgelist_reader.read(filepath)
    node_mapping = edgelist_reader.getNodeMap()
    logger(f"Graph reading completed in {time.time() - start_time:.2f} seconds")
    monitor_resources(logger)
    return nk_graph, node_mapping


def read_clustering(filepath: str, logger) -> pd.DataFrame:
    """
    Reads the clustering from the given file and returns a DataFrame with node IDs and cluster IDs.

    Args:
        filepath: Path to the clustering file
        logger: Logger function for output
        
    Returns:
        pd.DataFrame: DataFrame containing node IDs and cluster IDs
    """
    start_time = time.time()
    # Use low_memory=False to avoid unnecessary type checking
    cluster_df = pd.read_csv(filepath, sep="\t", header=None, names=["node_id", "cluster_name"], low_memory=False)
    
    # Use factorize which is more efficient than mapping unique values
    cluster_df['cluster_id'], _ = pd.factorize(cluster_df['cluster_name'])
    
    logger(f"Clustering reading completed in {time.time() - start_time:.2f} seconds")
    monitor_resources(logger)
    return cluster_df[['node_id', 'cluster_id']]


def process_edge_batch(edge_batch: List[Tuple], numerical_to_string_mapping: Dict, 
                       node_to_cluster_dict: Dict, num_clusters: int) -> csr_matrix:
    """
    Process a batch of edges to compute the inter-cluster probabilities.
    
    Args:
        edge_batch: List of edges to process
        numerical_to_string_mapping: Mapping from numerical node IDs to string IDs
        node_to_cluster_dict: Mapping from node IDs to cluster IDs
        num_clusters: Number of clusters in the graph
        
    Returns:
        A sparse matrix representing the probabilities for this batch
    """
    local_probs = lil_matrix((num_clusters, num_clusters), dtype=int)
    
    for edge in edge_batch:
        source = numerical_to_string_mapping.get(edge[0])
        target = numerical_to_string_mapping.get(edge[1])
        
        if source is not None and target is not None:
            source_cluster_idx = node_to_cluster_dict.get(source)
            target_cluster_idx = node_to_cluster_dict.get(target)
            
            if source_cluster_idx is not None and target_cluster_idx is not None:
                local_probs[source_cluster_idx, target_cluster_idx] += 1
                
                # Only add the reverse edge if it's different
                if source_cluster_idx != target_cluster_idx:
                    local_probs[target_cluster_idx, source_cluster_idx] += 1
    
    return local_probs.tocsr()


def get_probs_threaded(G_c: nk.Graph, node_mapping: Dict, cluster_df: pd.DataFrame, n_threads: int, logger) -> csr_matrix:
    """
    Compute the probability matrix using threads for the given graph and clustering.

    Args:
        G_c: Networkit graph object
        node_mapping: Node mapping from string to integer
        cluster_df: DataFrame containing node IDs and cluster IDs
        n_threads: Number of threads to use
        logger: Logger function for output
        
    Returns:
        csr_matrix: Sparse matrix representing the probabilities
    """
    start_time = time.time()
    
    # Pre-compute all mappings before the loop
    numerical_to_string_mapping = {v: int(k) for k, v in node_mapping.items()}
    
    # Create node_id to cluster_id mapping more efficiently
    node_to_cluster_dict = pd.Series(cluster_df['cluster_id'].values, 
                                     index=cluster_df['node_id']).to_dict()
    
    num_clusters = cluster_df['cluster_id'].nunique()
    
    # Get all edges from the graph
    all_edges = list(G_c.iterEdges())
    
    # Determine batch size based on number of edges and threads
    batch_size = max(1, len(all_edges) // (n_threads * 10))  # Create more batches than threads
    edge_batches = [all_edges[i:i + batch_size] for i in range(0, len(all_edges), batch_size)]
    
    logger(f"Processing {len(all_edges)} edges in {len(edge_batches)} batches using {n_threads} threads")
    
    # Process edge batches using thread pool
    results = []
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for batch in edge_batches:
            future = executor.submit(
                process_edge_batch, 
                batch, 
                numerical_to_string_mapping, 
                node_to_cluster_dict, 
                num_clusters
            )
            futures.append(future)
        
        # Show progress if verbose
        if logger(None, return_verbose=True):
            for future in tqdm(futures, desc="Computing probabilities"):
                results.append(future.result())
        else:
            for future in futures:
                results.append(future.result())
    
    # Combine results
    probs = results[0]
    for res in results[1:]:
        probs = probs + res
    
    logger(f"Probability matrix computation completed in {time.time() - start_time:.2f} seconds")
    monitor_resources(logger)
    return probs


def compute_node_degree(node_id: str, node_mapping: Dict, G_c: nk.Graph) -> int:
    """
    Compute the degree of a single node.
    
    Args:
        node_id: String identifier of the node
        node_mapping: Mapping from string node IDs to numerical IDs
        G_c: Networkit graph object
        
    Returns:
        Degree of the node
    """
    node_idx = node_mapping.get(node_id)
    if node_idx is not None:
        return G_c.degree(node_idx)
    return 0


def get_degree_sequence_threaded(cluster_df: pd.DataFrame, G_c: nk.Graph, 
                                node_mapping: Dict, n_threads: int, logger) -> List[int]:
    """
    Compute the degree sequence using threads.
    
    Args:
        cluster_df: DataFrame containing node IDs and cluster IDs
        G_c: Networkit graph object
        node_mapping: Mapping from string node IDs to numerical IDs
        n_threads: Number of threads to use
        logger: Logger function for output
        
    Returns:
        List of node degrees
    """
    start_time = time.time()
    
    # Pre-convert all node_ids to strings once
    node_ids = cluster_df['node_id'].astype(str).tolist()
    
    # Compute degrees using thread pool
    results = []
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {}
        for node_id in node_ids:
            future = executor.submit(compute_node_degree, node_id, node_mapping, G_c)
            futures[future] = node_id
        
        # Show progress if verbose
        if logger(None, return_verbose=True):
            for future in tqdm(futures, desc="Computing degrees"):
                results.append(future.result())
        else:
            for future in futures:
                results.append(future.result())
    
    logger(f"Degree sequence computation completed in {time.time() - start_time:.2f} seconds")
    monitor_resources(logger)
    return results


def process_edge_subset(edges: List[Tuple], node_idx_set: np.ndarray) -> Set[Tuple[int, int]]:
    """
    Process a subset of edges from the generated graph.
    
    Args:
        edges: List of edges to process
        node_idx_set: Array of node IDs
        
    Returns:
        Set of processed edges
    """
    edge_set = set()
    
    for edge in edges:
        source_idx, target_idx = edge[0], edge[1]
        # Ensure indices are valid
        if source_idx < len(node_idx_set) and target_idx < len(node_idx_set):
            source = node_idx_set[source_idx]
            target = node_idx_set[target_idx]
            # Store edges with source < target to avoid duplicates
            if source <= target:
                edge_set.add((source, target))
            else:
                edge_set.add((target, source))
    
    return edge_set


def extract_edges_threaded(N: gt.Graph, node_idx_set: np.ndarray, n_threads: int, logger) -> Set[Tuple[int, int]]:
    """
    Extract edges from the generated graph using threads.
    
    Args:
        N: Generated graph-tool graph
        node_idx_set: Array of node IDs
        n_threads: Number of threads to use
        logger: Logger function for output
        
    Returns:
        Set of edges in the graph
    """
    start_time = time.time()
    
    # Get all edges from the graph
    all_edges = list(N.get_edges())
    
    # Determine batch size based on number of edges and threads
    batch_size = max(1, len(all_edges) // (n_threads * 10))
    edge_batches = [all_edges[i:i + batch_size] for i in range(0, len(all_edges), batch_size)]
    
    logger(f"Processing {len(all_edges)} extracted edges in {len(edge_batches)} batches using {n_threads} threads")
    
    # Process edge batches using thread pool
    results = []
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for batch in edge_batches:
            future = executor.submit(process_edge_subset, batch, node_idx_set)
            futures.append(future)
        
        # Show progress if verbose
        if logger(None, return_verbose=True):
            for future in tqdm(futures, desc="Extracting edges"):
                results.append(future.result())
        else:
            for future in futures:
                results.append(future.result())
    
    # Combine results
    all_edges_set = set()
    for edge_set in results:
        all_edges_set.update(edge_set)
    
    logger(f"Edge extraction completed in {time.time() - start_time:.2f} seconds")
    logger(f"Extracted {len(all_edges_set)} unique edges")
    monitor_resources(logger)
    return all_edges_set


def save_generated_graph(edges_list: Set[Tuple[int, int]], out_edge_file: str, n_threads: int, logger) -> None:
    """
    Save the generated graph edges to a file.
    
    Args:
        edges_list: Set of edges in the graph
        out_edge_file: Output file path to save the edges
        n_threads: Number of threads to use
        logger: Logger function for output
    """
    start_time = time.time()
    
    # Convert to list for consistent ordering
    edges_list_sorted = sorted(edges_list)
    
    # Thread-safe queue for writing chunks
    write_queue = queue.Queue()
    
    # Function to format chunks of edges for writing
    def format_chunk(chunk_id, chunk):
        lines = [f"{source}\t{target}\n" for source, target in chunk]
        write_queue.put((chunk_id, ''.join(lines)))
    
    # Determine chunk size based on number of edges
    chunk_size = max(10000, len(edges_list_sorted) // n_threads)
    num_chunks = (len(edges_list_sorted) + chunk_size - 1) // chunk_size
    
    # Process chunks in parallel
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min(start_idx + chunk_size, len(edges_list_sorted))
            chunk = edges_list_sorted[start_idx:end_idx]
            executor.submit(format_chunk, i, chunk)
    
    # Write chunks to file in correct order
    with open(out_edge_file, 'w') as f:
        # Wait for all chunks to be processed
        chunks_written = 0
        processed_chunks = {}
        
        # Show progress if verbose
        progress = tqdm(total=num_chunks, desc="Writing file") if logger(None, return_verbose=True) else None
        
        while chunks_written < num_chunks:
            try:
                chunk_id, chunk_text = write_queue.get(timeout=0.1)
                processed_chunks[chunk_id] = chunk_text
                
                # Write chunks in order
                while chunks_written in processed_chunks:
                    f.write(processed_chunks[chunks_written])
                    del processed_chunks[chunks_written]
                    chunks_written += 1
                    if progress:
                        progress.update(1)
                
            except queue.Empty:
                # No chunks available yet, continue waiting
                continue
    
    if progress:
        progress.close()
    
    logger(f"Graph saving completed in {time.time() - start_time:.2f} seconds")
    monitor_resources(logger)


def main(
        edge_input: str = typer.Option(..., "--filepath", "-f"),
        cluster_input: str = typer.Option(..., "--cluster_filepath", "-c"),
        output_dir: str = typer.Option("", "--output_directory", "-o"),
        n_threads: int = typer.Option(-1, "--jobs", "-j"),
        verbose: bool = typer.Option(False, "--verbose", "-v")):
    """
    Main function to generate a synthetic graph based on the input edge list and clustering.

    Args:
        edge_input: Path to the edge list file.
        cluster_input: Path to the clustering file.
        output_dir: Output directory to save the generated graph.
        n_threads: Number of threads to use (-1 for all available cores).
        verbose: Whether to enable verbose output.
    """
    total_start_time = time.time()
    
    # Create a logger that respects the verbose flag
    logger = create_logger(verbose)
    
    # Set up thread count
    if n_threads <= 0:
        n_threads = os.cpu_count() or 4
        
    logger(f"Starting SBM graph generation with {n_threads} threads")
    monitor_resources(logger)
    
    # Use pathlib for more robust path handling
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Output filename
    out_edge_file = output_path / 'syn_sbm.tsv'

    logger(f"Processing input files: {edge_input} and {cluster_input}")
    
    # Read the graph and clustering
    G, node_mapping = read_graph(edge_input, logger)
    cluster_df = read_clustering(cluster_input, logger)

    # Get the inter-cluster probabilities and degree sequence using threads
    probs = get_probs_threaded(G, node_mapping, cluster_df, n_threads, logger)
    cluster_assignment = cluster_df['cluster_id'].to_numpy()
    out_deg_seq = get_degree_sequence_threaded(cluster_df, G, node_mapping, n_threads, logger)

    # Generate the synthetic graph
    sbm_start_time = time.time()
    logger("Generating SBM graph...")
    N = gt.generate_sbm(cluster_assignment, probs, out_degs=out_deg_seq, micro_ers=True, micro_degs=True)
    logger(f"SBM generation completed in {time.time() - sbm_start_time:.2f} seconds")
    monitor_resources(logger)

    # Get the edges of the generated graph using threads
    node_idx_set = cluster_df['node_id'].to_numpy()
    N_edge_list = extract_edges_threaded(N, node_idx_set, n_threads, logger)

    # Save the generated graph edges to a file
    save_generated_graph(N_edge_list, str(out_edge_file), n_threads, logger)
    
    # Always print the final completion message
    print(f"Generated graph saved to: {out_edge_file}")
    logger(f"Total execution time: {time.time() - total_start_time:.2f} seconds", always=True)
    monitor_resources(logger)


if __name__ == "__main__":
    typer.run(main)
    