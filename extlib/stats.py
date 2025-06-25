import typer
import networkit as nk
import pandas as pd
import os
import time
import json
from typing import Dict

from hm01.graph import Graph, IntangibleSubgraph
from hm01.mincut import viecut

def from_existing_clustering(filepath) -> Dict[str, IntangibleSubgraph]:
    '''
    Load cluster assignments from a file and create IntangibleSubgraph objects.
    
    Args:
        filepath: Path to the clustering file (node_id cluster_id format)
        
    Returns:
        Dictionary mapping cluster IDs to IntangibleSubgraph objects containing nodes in that cluster.
        Only returns clusters with more than one node.
    '''
    # Initialize empty dictionary to store clusters
    clusters: Dict[str, IntangibleSubgraph] = {}
    
    # Parse the clustering file line by line
    with open(filepath) as f:
        for line in f:
            node_id, cluster_id = line.split()
            # Create a new cluster if it doesn't exist, or add node to existing cluster
            clusters.setdefault(
                cluster_id, IntangibleSubgraph([], cluster_id)
            ).subset.append(int(node_id))
    
    # Filter out singleton clusters (clusters with only one node)
    return {key: val for key, val in clusters.items() if val.n() > 1}

def main(
    input: str = typer.Option(..., "--input", "-i", help="Input graph file path (tab-separated edge list)"),
    existing_clustering: str = typer.Option(..., "--existing-clustering", "-e", help="Path to existing clustering file"),
    output: str = typer.Option("", "--output", "-o", help="Output file path for statistics (default: input_stats.csv)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose output")
): 
    '''
    Calculate graph statistics for clusters in a clustering.
    
    For each non-singleton cluster, calculates:
    - Number of nodes (n)
    - Number of edges (m)
    - Connectivity (minimum cut value)
    
    Results are saved as a CSV file.
    '''
    if verbose:
        print(f"Starting cluster analysis with inputs:")
        print(f"- Graph file: {input}")
        print(f"- Clustering file: {existing_clustering}")
        start_time = time.time()
    
    # Determine output file path if not specified
    if output == "":
        base, _ = os.path.splitext(existing_clustering)
        outfile = base + '_stats.csv'
    else:
        outfile = output
        
    if verbose:
        print(f"Output will be written to: {outfile}")
        print("Loading clustering file...")
        cluster_start_time = time.time()
    
    # Load the clusters from the clustering file
    cluster_dict = from_existing_clustering(existing_clustering)
    clusters = cluster_dict.values()
    
    if verbose:
        cluster_time = time.time() - cluster_start_time
        print(f"Loaded {len(clusters)} non-singleton clusters in {cluster_time:.2f} seconds")
        print("Extracting basic cluster properties...")
    
    # Extract cluster properties
    ids = [cluster.index for cluster in clusters]
    ns = [cluster.n() for cluster in clusters]
    
    if verbose:
        print(f"Loading graph from {input}...")
        graph_start_time = time.time()

    # Load full graph into NetworkIt Graph object
    edgelist_reader = nk.graphio.EdgeListReader("\t", 0)
    nk_graph = edgelist_reader.read(input)
    global_graph = Graph(nk_graph, "")
    
    if verbose:
        graph_time = time.time() - graph_start_time
        print(f"Loaded graph with {nk_graph.numberOfNodes()} nodes and {nk_graph.numberOfEdges()} edges in {graph_time:.2f} seconds")
        print("Calculating edge counts for each cluster...")
        edge_start_time = time.time()
    
    # Count edges in each cluster
    ms = [cluster.count_edges(global_graph) for cluster in clusters]
    
    if verbose:
        edge_time = time.time() - edge_start_time
        print(f"Edge counting completed in {edge_time:.2f} seconds")
        print("Realizing clusters as subgraphs...")
    
    # Convert IntangibleSubgraphs to actual graph objects
    clusters = [cluster.realize(global_graph) for cluster in clusters]
    
    if verbose:
        print("Calculating minimum cuts for each cluster (this may take some time)...")
        mincut_start_time = time.time()
    
    # Calculate minimum cut for each cluster
    mincut_results = [viecut(cluster) for cluster in clusters]
    mincuts = [result[-1] for result in mincut_results]
    
    if verbose:
        mincut_time = time.time() - mincut_start_time
        print(f"Minimum cut calculations completed in {mincut_time:.2f} seconds")
        print("Creating output dataframe...")
    
    # Create and save results dataframe
    df = pd.DataFrame(list(zip(ids, ns, ms, mincuts)),
        columns=['cluster', 'n', 'm', 'connectivity'])
    
    df.to_csv(outfile, index=False)
    
    if verbose:
        total_time = time.time() - start_time
        print(f"Results saved to {outfile}")
        print(f"Total processing time: {total_time:.2f} seconds")
        
        # Print summary statistics
        print("\nSummary Statistics:")
        print(f"Number of clusters analyzed: {len(clusters)}")
        print(f"Average nodes per cluster: {sum(ns)/len(ns):.2f}")
        print(f"Average edges per cluster: {sum(ms)/len(ms):.2f}")
        print(f"Average connectivity: {sum(mincuts)/len(mincuts):.2f}")
        print(f"Min connectivity: {min(mincuts)}")
        print(f"Max connectivity: {max(mincuts)}")


def entry_point():
    '''Entry point for command line interface'''
    typer.run(main)

if __name__ == "__main__":
    entry_point()
