"""
Hybrid SBM script that uses the original working strategy with new optimizations.
Based on the original plusEdges approach but with performance improvements.

Usage:
    python gen_SBM_degseq_hybrid.py -f <current_edges> -c <clustering> -ef <empirical_edges> -o <output_dir> [-j <threads>] [-v]
"""

import json
import os
import queue
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Set, Tuple

import graph_tool.all as gt
import numpy as np
import pandas as pd
import psutil
import typer
from scipy.sparse import coo_matrix
from tqdm import tqdm
import networkit as nk


def create_logger(verbose: bool):
    """Create a logger function that only prints when verbose mode is enabled."""
    print_lock = threading.Lock()
    
    def log(message: str, always: bool = False, return_verbose: bool = False):
        if return_verbose:
            return verbose
        if verbose or always:
            with print_lock:
                print(message)
    
    return log


def monitor_resources(logger):
    """Monitor and log current memory usage."""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    logger(f"Memory usage: {memory_info.rss / (1024 * 1024):.2f} MB")


def read_graph(input_network, logger):
    """Read the network using NetworkIt."""
    start_time = time.time()
    
    elr = nk.graphio.EdgeListReader('\t', 0, continuous=False, directed=False)
    graph = elr.read(input_network)
    graph.removeMultiEdges()
    graph.removeSelfLoops()
    node_mapping_dict = elr.getNodeMap()
    
    logger(f"Graph reading completed in {time.time() - start_time:.2f} seconds")
    logger(f"Loaded graph with {graph.numberOfNodes()} nodes and {graph.numberOfEdges()} edges")
    monitor_resources(logger)
    
    return graph, node_mapping_dict


def read_clustering(input_clustering, logger):
    """Read the clustering file."""
    start_time = time.time()
    
    cluster_df = pd.read_csv(input_clustering, sep="\t", header=None, names=[
                             "node_id", "cluster_name"], dtype=str)
    unique_values = cluster_df['cluster_name'].unique()
    value_map = {
        value: idx
        for idx, value in enumerate(unique_values)
    }
    cluster_df['cluster_id'] = cluster_df['cluster_name'].map(value_map)
    
    logger(f"Clustering reading completed in {time.time() - start_time:.2f} seconds")
    logger(f"Loaded {len(cluster_df)} nodes with {cluster_df['cluster_id'].nunique()} clusters")
    monitor_resources(logger)
    
    return cluster_df, value_map


def get_probs_optimized(G_c, node_mapping, cluster_df, logger, n_threads=None):
    """Optimized probability matrix computation."""
    start_time = time.time()
    
    if n_threads is None:
        n_threads = os.cpu_count() or 4
    
    numerical_to_string_mapping = {v: k for k, v in node_mapping.items()}
    cluster_ids, counts = np.unique(cluster_df['cluster_id'], return_counts=True)
    num_clusters = len(cluster_ids)
    
    node_to_cluster_dict = cluster_df.set_index('node_id')['cluster_id'].to_dict()

    # Collect all edges first
    edges_list = []
    for edge in G_c.iterEdges():
        source = numerical_to_string_mapping.get(edge[0])
        target = numerical_to_string_mapping.get(edge[1])
        source_cluster_idx = node_to_cluster_dict.get(source)
        target_cluster_idx = node_to_cluster_dict.get(target)
        if source_cluster_idx is not None and target_cluster_idx is not None:
            edges_list.append((source_cluster_idx, target_cluster_idx))
    
    # Count cluster pairs using Counter (optimized)
    pair_counts = Counter(edges_list)
    
    # Build COO matrix
    def chunk_dict(d, n_chunks):
        items = list(d.items())
        chunk_size = len(items) // n_chunks + 1
        return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]
    
    def process_chunk(chunk_pairs):
        local_rows, local_cols, local_data = [], [], []
        
        for (i, j), count in chunk_pairs:
            count = int(count)
            local_rows.extend([i, j])
            local_cols.extend([j, i])
            local_data.extend([count, count])
        
        return local_rows, local_cols, local_data
    
    chunks = chunk_dict(pair_counts, n_threads)
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        results = list(executor.map(process_chunk, chunks))
    
    all_rows, all_cols, all_data = [], [], []
    for rows, cols, data in results:
        all_rows.extend(rows)
        all_cols.extend(cols)
        all_data.extend(data)
    
    probs = coo_matrix(
        (all_data, (all_rows, all_cols)),
        shape=(num_clusters, num_clusters),
        dtype=np.int32
    ).tocsr()
    
    total_edges = probs.trace()//2 + (probs.sum() - probs.trace())//2
    logger(f"Probability matrix computation completed in {time.time() - start_time:.2f} seconds")
    logger(f"Number of edges as per probs matrix: {total_edges}")
    monitor_resources(logger)
    
    return probs


def save_new_edges_optimized(new_edges: Set[Tuple[str, str]], out_edge_file: str, n_threads: int, logger) -> None:
    """Optimized edge saving with threading."""
    start_time = time.time()
    
    edges_list_sorted = sorted(new_edges)
    
    write_queue = queue.Queue()
    
    def format_chunk(chunk_id, chunk):
        lines = [f"{source}\t{target}\n" for source, target in chunk]
        write_queue.put((chunk_id, ''.join(lines)))
    
    chunk_size = max(10000, len(edges_list_sorted) // n_threads)
    num_chunks = (len(edges_list_sorted) + chunk_size - 1) // chunk_size
    
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min(start_idx + chunk_size, len(edges_list_sorted))
            chunk = edges_list_sorted[start_idx:end_idx]
            executor.submit(format_chunk, i, chunk)
    
    with open(out_edge_file, 'w') as f:
        chunks_written = 0
        processed_chunks = {}
        
        progress = tqdm(total=num_chunks, desc="Writing file") if logger(None, return_verbose=True) else None
        
        while chunks_written < num_chunks:
            try:
                chunk_id, chunk_text = write_queue.get(timeout=0.1)
                processed_chunks[chunk_id] = chunk_text
                
                while chunks_written in processed_chunks:
                    f.write(processed_chunks[chunks_written])
                    del processed_chunks[chunks_written]
                    chunks_written += 1
                    if progress:
                        progress.update(1)
                
            except queue.Empty:
                continue
    
    if progress:
        progress.close()
    
    logger(f"New edges saving completed in {time.time() - start_time:.2f} seconds")
    monitor_resources(logger)


def main(
        edge_input: str = typer.Option(..., "--filepath", "-f"),
        cluster_input: str = typer.Option(..., "--cluster_filepath", "-c"), 
        empirical_edge_input: str = typer.Option(..., "--empirical_filepath", "-ef"),
        output_dir: str = typer.Option("", "--output_directory", "-o"),
        n_threads: int = typer.Option(-1, "--jobs", "-j"),
        verbose: bool = typer.Option(False, "--verbose", "-v")):
    """
    Generate edges using the original plusEdges strategy with optimizations.

    Args:
        edge_input: Path to current (synthetic) edge list file
        cluster_input: Path to clustering file
        empirical_edge_input: Path to empirical (target) edge list file
        output_dir: Output directory for new edges
        n_threads: Number of threads (-1 for all cores)
        verbose: Enable verbose output
    """
    total_start_time = time.time()
    
    logger = create_logger(verbose)
    
    if n_threads <= 0:
        n_threads = os.cpu_count() or 4
        
    logger(f"Starting hybrid SBM degree sequence gap filling with {n_threads} threads")
    monitor_resources(logger)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    out_edge_file = output_path / 'new_edges.tsv'

    try:
        # Step 1: Read empirical graph
        logger("Reading empirical graph...")
        emp_graph, emp_node_mapping = read_graph(empirical_edge_input, logger)
        
        # Step 2: Read clustering
        logger("Reading graph clustering...")
        cluster_df, cluster_mapping_dict = read_clustering(cluster_input, logger)
        
        # Create clustering dictionary
        clustering_dict = dict(zip(cluster_df['node_id'], cluster_df['cluster_id']))
        
        # Step 3: Get edge probs matrix for empirical network
        logger("Computing edge probability matrix for empirical network...")
        emp_clustered_nodes = [emp_node_mapping[v] for v in clustering_dict.keys() if v in emp_node_mapping]
        G_c = nk.graphtools.subgraphFromNodes(emp_graph, emp_clustered_nodes)
        emp_probs = get_probs_optimized(G_c, emp_node_mapping, cluster_df, logger, n_threads)
        
        # Step 4: Read current synthetic graph
        logger("Reading current synthetic graph...")
        graph, node_mapping = read_graph(edge_input, logger)
        
        clustered_nodes = [node_mapping[v] for v in clustering_dict.keys() if v in node_mapping]
        node_mapping_reversed = {u: str(v) for v, u in node_mapping.items()}
        
        # Step 5: Find nodes with available degrees (ORIGINAL LOGIC)
        logger("Finding nodes with available degrees...")
        start_time = time.time()
        
        emp_degrees = dict()
        sbm_degrees = dict()
        available_node_degrees = dict()
        
        for c_node in clustering_dict.keys():
            if c_node in emp_node_mapping and c_node in node_mapping:
                emp_node = emp_node_mapping[c_node]
                syn_node = node_mapping[c_node]
                emp_degrees[c_node] = emp_graph.degree(emp_node)
                sbm_degrees[c_node] = graph.degree(syn_node)
                deg_diff = emp_degrees[c_node] - sbm_degrees[c_node]
                if deg_diff > 0:
                    available_node_degrees[syn_node] = deg_diff
        
        logger(f"Available degree computation completed in {time.time() - start_time:.2f} seconds")
        
        if not available_node_degrees:
            logger("No nodes need additional edges!", always=True)
            with open(out_edge_file, 'w') as f:
                pass
            return
        
        # Step 6: Get edge probs matrix for synthetic network
        logger("Computing edge probability matrix for synthetic network...")
        G_c_wc = nk.graphtools.subgraphFromNodes(graph, clustered_nodes)
        sbm_wc_probs = get_probs_optimized(G_c_wc, node_mapping, cluster_df, logger, n_threads)
        
        # Step 7: Get probs matrix for remaining edges (ORIGINAL LOGIC)
        logger("Computing probability matrix for remaining edges...")
        probs = emp_probs - sbm_wc_probs
        # Set negative values to 0
        probs.data[probs.data < 0] = 0
        
        # Step 8: Prepare SBM generation (ORIGINAL LOGIC)
        logger("Preparing SBM generation...")
        total_available_nodes = len(available_node_degrees.keys())
        available_nodes_list = []
        available_degrees = []
        
        for node, degree in available_node_degrees.items():
            available_nodes_list.append(node)
            # CRITICAL: Apply original degree capping
            available_degrees.append(min(total_available_nodes - 1, degree))
        
        avail_nodes_block_assignment = []
        for syn_node in available_nodes_list:
            avail_nodes_block_assignment.append(clustering_dict[node_mapping_reversed[syn_node]])
        
        # Filter probability matrix (ORIGINAL LOGIC)
        clusters_to_be_considered = np.unique(avail_nodes_block_assignment)
        for i in range(probs.shape[0]):
            for ind in range(probs.indptr[i], probs.indptr[i + 1]):
                j = probs.indices[ind]
                if i not in clusters_to_be_considered:
                    probs.data[ind] = 0
                if j not in clusters_to_be_considered:
                    probs.data[ind] = 0
        
        # Validation (ORIGINAL LOGIC)
        total_avail_edges = probs.trace()//2 + (probs.sum() - probs.trace())//2
        total_degrees = sum(available_degrees)//2
        
        logger(f"SBM Validation:")
        logger(f"  - Total available nodes: {total_available_nodes}")
        logger(f"  - Total available edges: {total_avail_edges}")
        logger(f"  - Total degrees needed: {total_degrees}")
        
        # Step 9: Generate SBM (ORIGINAL LOGIC)
        if total_avail_edges > 0 and total_degrees > 0:
            logger("Generating SBM...")
            sbm_start_time = time.time()
            
            try:
                generated_graph = gt.generate_sbm(avail_nodes_block_assignment, probs, available_degrees, 
                                               directed=False)
                
                gt.remove_parallel_edges(generated_graph)
                gt.remove_self_loops(generated_graph)
                
                logger(f"SBM generation completed in {time.time() - sbm_start_time:.2f} seconds")
                logger(f"Total number of edges added: {generated_graph.num_edges()}")
                
                # Step 10: Extract and save new edges (OPTIMIZED)
                logger("Extracting new edges...")
                new_edges = set()
                for edge in generated_graph.iter_edges():
                    node1 = node_mapping_reversed[available_nodes_list[edge[0]]]
                    node2 = node_mapping_reversed[available_nodes_list[edge[1]]]
                    new_edges.add((node1, node2))
                
                # Save with optimized threading
                save_new_edges_optimized(new_edges, str(out_edge_file), n_threads, logger)
                
                print(f"Generated {len(new_edges)} new edges saved to: {out_edge_file}")
                
            except RuntimeError as e:
                logger(f"SBM generation failed: {e}", always=True)
                logger("Creating empty output file", always=True)
                with open(out_edge_file, 'w') as f:
                    pass
        else:
            logger("No valid edges can be generated", always=True)
            with open(out_edge_file, 'w') as f:
                pass
        
    except Exception as e:
        logger(f"Error occurred: {e}", always=True)
        import traceback
        traceback.print_exc()
        # Create empty output file on error
        with open(out_edge_file, 'w') as f:
            pass
    
    logger(f"Total execution time: {time.time() - total_start_time:.2f} seconds", always=True)
    monitor_resources(logger)


if __name__ == "__main__":
    typer.run(main)
