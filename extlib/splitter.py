import pandas as pd
import typer
import time
from pathlib import Path

def main(
        edge_input: str = typer.Option(..., "--filepath", "-f"),
        cluster_input: str = typer.Option(..., "--cluster_filepath", "-c"),
        output_dir: str = typer.Option("", "--output_directory", "-o"),
        verbose: bool = typer.Option(False, "--verbose", "-v")):
    
    if verbose:
        print(f"Starting processing with inputs:\n- Edge file: {edge_input}\n- Cluster file: {cluster_input}\n- Output directory: {output_dir}")
        start_time = time.time()
    
    # Create output directory if it doesn't exist
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        if verbose:
            print(f"Output directory verified: {output_dir}")
    
    # Read the edge and cluster tsvs
    if verbose:
        print("Reading input files...")
        file_read_start = time.time()
    
    edge_df = pd.read_csv(edge_input, sep='\t', header=None)
    cluster_df = pd.read_csv(cluster_input, sep='\t', header=None)
    
    if verbose:
        file_read_time = time.time() - file_read_start
        print(f"Files read in {file_read_time:.2f} seconds")
        print(f"Edge file shape: {edge_df.shape}")
        print(f"Cluster file shape: {cluster_df.shape}")
    
    # Rename columns for clarity
    edge_df.columns = ['source', 'target']
    cluster_df.columns = ['node_id', 'cluster_id']
    
    if verbose:
        print("Identifying non-singleton clusters...")
        cluster_start = time.time()
    
    # Get the count of nodes in each cluster
    cluster_counts = cluster_df['cluster_id'].value_counts().reset_index()
    # Get all nodes that are in non-singleton clusters
    non_singleton_clusters = cluster_counts[cluster_counts['count'] > 1]['cluster_id'].tolist()
    
    if verbose:
        print(f"Found {len(non_singleton_clusters)} non-singleton clusters")
        nodes_start = time.time()
    
    non_singleton_nodes = cluster_df[cluster_df['cluster_id'].isin(non_singleton_clusters)]['node_id'].tolist()
    
    if verbose:
        nodes_time = time.time() - nodes_start
        print(f"Identified {len(non_singleton_nodes)} nodes in non-singleton clusters in {nodes_time:.2f} seconds")
        edges_start = time.time()
    
    # Get the edges that are between non-singleton clusters
    non_singleton_edges = edge_df[
        (edge_df['source'].isin(non_singleton_nodes)) & (edge_df['target'].isin(non_singleton_nodes))
    ]
    
    if verbose:
        print(f"Found {len(non_singleton_edges)} edges between non-singleton clusters")
    
    # Get the complement of the non-singleton edges
    singleton_edges = edge_df[
        ~((edge_df['source'].isin(non_singleton_nodes)) & (edge_df['target'].isin(non_singleton_nodes)))
    ]
    
    if verbose:
        edges_time = time.time() - edges_start
        print(f"Edge filtering completed in {edges_time:.2f} seconds")
        print(f"Found {len(singleton_edges)} singleton edges")
        clusters_start = time.time()
    
    # NEW: Split the clusterings based on which graph they belong to
    # Get all unique nodes in the non-singleton edges graph
    non_singleton_graph_nodes = set(pd.concat([non_singleton_edges['source'], non_singleton_edges['target']]))
    
    # Get all unique nodes in the singleton edges graph
    singleton_graph_nodes = set(pd.concat([singleton_edges['source'], singleton_edges['target']]))
    
    # Create the corresponding cluster files
    non_singleton_clusters_df = cluster_df[cluster_df['node_id'].isin(non_singleton_graph_nodes)]
    singleton_clusters_df = cluster_df[cluster_df['node_id'].isin(singleton_graph_nodes)]
    
    if verbose:
        clusters_time = time.time() - clusters_start
        print(f"Cluster splitting completed in {clusters_time:.2f} seconds")
        print(f"Non-singleton graph has {len(non_singleton_graph_nodes)} nodes and {len(non_singleton_clusters_df)} cluster assignments")
        print(f"Singleton graph has {len(singleton_graph_nodes)} nodes and {len(singleton_clusters_df)} cluster assignments")
        save_start = time.time()
    
    # Save edges as no-header tsvs
    non_singleton_edges_out = f"{output_dir}/non_singleton_edges.tsv"
    singleton_edges_out = f"{output_dir}/singleton_edges.tsv"
    non_singleton_edges.to_csv(non_singleton_edges_out, sep='\t', index=False, header=False)
    singleton_edges.to_csv(singleton_edges_out, sep='\t', index=False, header=False)
    
    # Save cluster assignments as no-header tsvs
    non_singleton_clusters_out = f"{output_dir}/non_singleton_clusters.tsv"
    singleton_clusters_out = f"{output_dir}/singleton_clusters.tsv"
    non_singleton_clusters_df.to_csv(non_singleton_clusters_out, sep='\t', index=False, header=False)
    singleton_clusters_df.to_csv(singleton_clusters_out, sep='\t', index=False, header=False)
    
    if verbose:
        save_time = time.time() - save_start
        total_time = time.time() - start_time
        print(f"Files saved in {save_time:.2f} seconds")
        print(f"Output files:")
        print(f"- {non_singleton_edges_out}")
        print(f"- {singleton_edges_out}")
        print(f"- {non_singleton_clusters_out}")
        print(f"- {singleton_clusters_out}")
        print(f"Total processing time: {total_time:.2f} seconds")

if __name__ == "__main__":
    typer.run(main)
    