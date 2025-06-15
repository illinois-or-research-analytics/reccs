import json
import typer
import pandas as pd
import matplotlib.pyplot as plt # type: ignore
import numpy as np

def get_min_degree_mapping(degseq_json):
    """
    Extracts the minimum degree for each cluster from the degree sequence JSON.
    """
    min_degree_mapping = {}
    for cluster_id, degrees in degseq_json.items():
        if degrees:  # Ensure there are degrees to process
            min_degree_mapping[cluster_id] = min(degrees)
        else:
            min_degree_mapping[cluster_id] = None  # Handle empty degree lists
    return min_degree_mapping

def get_degseq_euclidean_distance(cluster_id, output_degrees, reference_degrees):
    """
    Calculates the Euclidean distance between the output and reference degree sequences for a given cluster.
    """
    if cluster_id not in output_degrees or cluster_id not in reference_degrees:
        return None  # Cluster ID not found in either sequence

    output_set = np.array(output_degrees[cluster_id])
    reference_set = np.array(reference_degrees[cluster_id])

    # Calculate the Euclidean distance
    distance = np.sqrt(np.sum((output_set - reference_set) ** 2))
    return distance

def run(
    stats_dir: str,
    degseq_dir: str,
    reference_stats_dir: str,
    reference_degseq_dir: str,
    sbm_stats_dir: str,
    sbm_degseq_dir: str
):
    """
    Checks validity of RECCS output stats
    """
    # Load the stats file csv
    stats_df = pd.read_csv(stats_dir)
    reference_stats_df = pd.read_csv(reference_stats_dir)
    sbm_stats_df = pd.read_csv(sbm_stats_dir)

    # Load the degree sequence json
    with open(degseq_dir, 'r') as f:
        degseq = json.load(f)

    with open(reference_degseq_dir, 'r') as f:
        reference_degseq = json.load(f)

    with open(sbm_degseq_dir, 'r') as f:
        sbm_degseq = json.load(f)

    # Create cluster_id : connectivity mapping from the stats dataframes
    stats_connectivity = dict(zip(stats_df['cluster'], stats_df['connectivity']))
    reference_connectivity = dict(zip(reference_stats_df['cluster'], reference_stats_df['connectivity']))

    # Get the minimum degree mapping for both current and reference degree sequences
    min_degree_mapping = get_min_degree_mapping(degseq)

    # Check if any disconnected clusters exist in the stats connectivity mapping
    disconnected_clusters = [cluster for cluster, connectivity in stats_connectivity.items() if connectivity == 0]
    if disconnected_clusters:
        return False

    # Check if the connectivity matches the minimum degree for each cluster

    for cluster_id, min_degree in min_degree_mapping.items():
        cluster_id = int(cluster_id)
        if cluster_id in stats_connectivity:
            connectivity = reference_connectivity[cluster_id]
            if connectivity > min_degree:
                return False
        else:
            False

    # Check if connectivity matches between RECCS and reference stats

    for cluster_id, connectivity in stats_connectivity.items():
        cluster_id = int(cluster_id)
        if cluster_id in reference_connectivity:
            reference_connectivity_value = reference_connectivity[cluster_id]
            if connectivity < reference_connectivity_value:
                return False
        else:
            return False
    
    return True