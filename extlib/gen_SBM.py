"""
This script generates a synthetic graph based on an input edge list and clustering information.
It uses the graph-tool library to create a stochastic block model (SBM) graph.

The script reads the input edge list and clustering file, computes the probability matrix,
and generates a synthetic graph using the SBM model. The generated graph is saved to a specified output directory.

Usage:
    python gen_SBM.py -f <edge_list_filepath> -c <clustering_filepath> -o <output_directory>

Arguments:
    -f, --filepath: Path to the input edge list file.
    -c, --cluster_filepath: Path to the clustering file.
    -o, --output_directory: Directory to save the generated graph.

Example:
    python gen_SBM.py -f edge_list.txt -c clustering.txt -o output_directory

Adapted from: https://github.com/illinois-or-research-analytics/lanne2_networks/blob/main/generate_synthetic_networks/gen_SBM.py
Original author: Lahari Anne (@lanne2)
"""

import pandas as pd
import numpy as np
import graph_tool.all as gt # type: ignore
from graph_tool.all import * # type: ignore
import typer
import os
import networkit as nk
from scipy.sparse import dok_matrix # type: ignore


def read_graph(filepath):
    """
    Reads the graph from the given edge list file and returns a Networkit graph and a node mapping.

    Args:
        filepath (str): Path to the edge list file.

    Returns: 
        nk.Graph: Networkit graph object.
        dict: Node mapping from string to integer.
    """
    edgelist_reader = nk.graphio.EdgeListReader("\t", 0, directed=False, continuous=False)
    nk_graph = edgelist_reader.read(filepath)
    node_mapping = edgelist_reader.getNodeMap()
    return nk_graph, node_mapping

def read_clustering(filepath):
    """
    Reads the clustering from the given file and returns a DataFrame with node IDs and cluster IDs.

    Args:
        filepath (str): Path to the clustering file.

    Returns:
        pd.DataFrame: DataFrame containing node IDs and cluster IDs.
    """
    cluster_df = pd.read_csv(filepath, sep="\t", header=None, names=["node_id", "cluster_name"])
    unique_values = cluster_df["cluster_name"].unique()
    value_map = {value: idx for idx, value in enumerate(unique_values)}
    cluster_df['cluster_id'] = cluster_df['cluster_name'].map(value_map)
    return cluster_df[['node_id', 'cluster_id']]

def get_probs(G_c, node_mapping, cluster_df):
    """
    Computes the probability matrix for the given graph and clustering.

    Args:
        G_c (nk.Graph): Networkit graph object.
        node_mapping (dict): Node mapping from string to integer.
        cluster_df (pd.DataFrame): DataFrame containing node IDs and cluster IDs.

    Returns:
        dok_matrix: Sparse matrix representing the probabilities.
    """
    # Convert node mapping to a dictionary with integer keys and string values
    numerical_to_string_mapping = {v: int(k) for k, v in node_mapping.items()}
    cluster_ids, _ = np.unique(cluster_df['cluster_id'], return_counts=True)
    num_clusters = len(cluster_ids)
    
    # Create a mapping from node IDs to cluster IDs
    node_to_cluster_dict = cluster_df.set_index('node_id')['cluster_id'].to_dict()

    # Compute the probability matrix
    probs = dok_matrix((num_clusters, num_clusters), dtype=int)
    for edge in G_c.iterEdges():
        source = numerical_to_string_mapping.get(edge[0])
        target = numerical_to_string_mapping.get(edge[1])
        source_cluster_idx = node_to_cluster_dict.get(source)
        target_cluster_idx = node_to_cluster_dict.get(target)
        probs[source_cluster_idx, target_cluster_idx] += 1
        probs[target_cluster_idx, source_cluster_idx ] += 1

    probs = probs.tocsr(copy=False)
    return probs

def get_degree_sequence(cluster_df, G_c, node_mapping,step, non_singleton_components):
    """
    Computes the degree sequence for the given graph and clustering.

    Args:
        cluster_df (pd.DataFrame): DataFrame containing node IDs and cluster IDs.
        G_c (nk.Graph): Networkit graph object.
        node_mapping (dict): Node mapping from string to integer.
        step (int): Step number for the degree sequence calculation.
        non_singleton_components (list): List of non-singleton components.

    Returns:
        list: Degree sequence for the graph.
    """
    deg_seq = []
    if step ==3:
        for _, row in cluster_df.iterrows():
            deg_seq.append(G_c.degree(node_mapping.get(str(row['node_id']))))
    elif step ==4:
        for component in non_singleton_components:
            deg_seq_component = []
            for node in component:
                deg_seq_component.append(G_c.degree(node))
            deg_seq.append(deg_seq_component)
    return deg_seq

def save_generated_graph(edges_list, out_edge_file):
    """
    Saves the generated graph edges to a file.

    Args:
        edges_list (list): List of edges in the graph.
        out_edge_file (str): Output file path to save the edges.
    """
    edge_df = pd.DataFrame(edges_list, columns=['source', 'target'])
    edge_df.to_csv(out_edge_file, sep='\t', index=False, header=None)

def main(
        edge_input: str = typer.Option(..., "--filepath", "-f"),
        cluster_input: str = typer.Option(..., "--cluster_filepath", "-c"),
        output_dir: str = typer.Option("", "--output_directory", "-o")):
    """
    Main function to generate a synthetic graph based on the input edge list and clustering.

    Args:
        edge_input (str): Path to the edge list file.
        cluster_input (str): Path to the clustering file.
        output_dir (str): Output directory to save the generated graph.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Output filename
    out_edge_file = os.path.join(output_dir, f'syn_sbm.tsv')

    # Read the graph and clustering
    G, node_mapping = read_graph(edge_input)
    cluster_df = read_clustering(cluster_input)

    # Add in outliers - this is mainly for Stage 2 of RECCS
    new_cluster_id = np.max(cluster_df['cluster_id']) + 1
    node_mapping_reversed = {v: int(k) for k, v in node_mapping.items()}
    clustered_nodes_id_org = cluster_df['node_id'].to_numpy()
    nodes_set = set()
    for u in G.iterNodes():
        nodes_set.add(node_mapping_reversed.get(u))
    unclustered_nodes = nodes_set.difference(clustered_nodes_id_org)
    
    # Add in the outliers to the cluster_df
    unclustered_node_cluster_mapping = []
    for v in unclustered_nodes:
        row = {'node_id': v, 'cluster_id' : new_cluster_id}
        unclustered_node_cluster_mapping.append(row)
        new_cluster_id += 1
    cluster_df = cluster_df._append(unclustered_node_cluster_mapping, ignore_index=True)
    cluster_df = cluster_df.reset_index()

    # Get the inter-cluster probabilities and degree sequence
    probs = get_probs(G, node_mapping, cluster_df)
    cluster_assignment = cluster_df['cluster_id'].to_numpy()
    out_deg_seq = get_degree_sequence(cluster_df, G, node_mapping, 3, [])

    # Generate the synthetic graph
    N = gt.generate_sbm(cluster_assignment, probs, out_degs=out_deg_seq, micro_ers=True, micro_degs=True)

    # Get the edges of the generated graph
    N_edge_list = set()
    node_idx_set = cluster_df['node_id'].tolist()
    for edge in N.get_edges():
        source = node_idx_set[min(edge[0],edge[1])]
        target = node_idx_set[max(edge[1], edge[0])]
        N_edge_list.add((source, target))

    # Save the generated graph edges to a file
    save_generated_graph(N_edge_list, out_edge_file)

if __name__ == "__main__":
    typer.run(main)
