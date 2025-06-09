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

def main(
    stats_dir: str = typer.Option(..., "--stats_output", "-s", help="Output stats filename"),
    degseq_dir: str = typer.Option(..., "--degree_sequence", "-d", help="Output degree sequence filename"),
    reference_stats_dir: str = typer.Option(..., "--reference_stats", "-rs", help="Reference stats filename"),
    reference_degseq_dir: str = typer.Option(..., "--reference_degree_sequence", "-rd", help="Reference degree sequence filename"),
    sbm_stats_dir: str = typer.Option(..., "--sbm_stats", "-ss", help="SBM stats filename"),
    sbm_degseq_dir: str = typer.Option(..., "--sbm_degree_sequence", "-sd", help="SBM degree sequence filename"),
    plot_output_path: str = typer.Option("degree_sequence_distances_boxplot.png", "--plot_output", "-p", help="Path to save the plot output")
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
    reference_min_degree_mapping = get_min_degree_mapping(reference_degseq)

    # Check if any disconnected clusters exist in the stats connectivity mapping
    disconnected_clusters = [cluster for cluster, connectivity in stats_connectivity.items() if connectivity == 0]
    if disconnected_clusters:
        print(f"FAIL: Disconnected clusters found in stats connectivity mapping: {disconnected_clusters}")
    else:
        print("PASS: No disconnected clusters found in stats connectivity mapping!")

    # Check if the connectivity matches the minimum degree for each cluster
    connectivity_mismatch = False
    for cluster_id, min_degree in min_degree_mapping.items():
        cluster_id = int(cluster_id)
        if cluster_id in stats_connectivity:
            connectivity = reference_connectivity[cluster_id]
            if connectivity > min_degree:
                print(f"FAIL: Mismatch for cluster {cluster_id}: "
                      f"Connectivity {connectivity} does not match minimum degree {min_degree}.")
                connectivity_mismatch = True
        else:
            print(f"FAIL: Cluster {cluster_id} not found in stats connectivity mapping.")

    if not connectivity_mismatch:
        print("PASS: Minimum degree matches connectivity for all clusters!")

    # Get a histogram of degseq distances between SBM and reference, and RECCS and reference
    sbm_distances = []
    reccs_distances = []

    for cluster_id in degseq.keys():
        sbm_distance = get_degseq_euclidean_distance(cluster_id, sbm_degseq, reference_degseq)
        reccs_distance = get_degseq_euclidean_distance(cluster_id, degseq, reference_degseq)

        if sbm_distance is not None:
            sbm_distances.append(sbm_distance)
        if reccs_distance is not None:
            reccs_distances.append(reccs_distance)
    
    # Plot the histograms on the same figure
    plt.figure(figsize=(10, 6))
    
    # Create box plot comparing SBM and RECCS distances
    data_to_plot = [sbm_distances, reccs_distances]
    labels = ['SBM Distances', 'RECCS Distances']
    
    plt.boxplot(data_to_plot, labels=labels, showfliers=False)
    plt.title('Comparison of Degree Sequence Distances')
    plt.ylabel('Euclidean Distance')
    plt.grid(True, alpha=0.3)
    plt.savefig(plot_output_path)

if __name__ == "__main__":
    typer.run(main)
    