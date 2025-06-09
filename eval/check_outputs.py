import json
import typer
import pandas as pd


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

def verify_no_degree_deficit(output_degseq, reference_degseq):
    """
    Verifies that the output degree sequence does not have a deficit compared to the reference degree sequence.
    """
    output_degseq = list(sorted(output_degseq))
    reference_degseq = list(sorted(reference_degseq))

    for o, r in zip(output_degseq, reference_degseq):
        if o < r:
            print(f"FAIL: Degree deficit found: output {o} is less than reference {r}.")
            return False
    
    return True

def main(
    stats_dir: str = typer.Option(..., "--stats_output", "-s", help="Output stats filename"),
    degseq_dir: str = typer.Option(..., "--degree_sequence", "-d", help="Output degree sequence filename"),
    reference_stats_dir: str = typer.Option(..., "--reference_stats", "-rs", help="Reference stats filename"),
    reference_degseq_dir: str = typer.Option(..., "--reference_degree_sequence", "-rd", help="Reference degree sequence filename")
):
    """
    Checks validity of RECCS output stats
    """
    # Load the stats file csv
    stats_df = pd.read_csv(stats_dir)
    reference_stats_df = pd.read_csv(reference_stats_dir)

    # Load the degree sequence json
    with open(degseq_dir, 'r') as f:
        degseq = json.load(f)

    with open(reference_degseq_dir, 'r') as f:
        reference_degseq = json.load(f)

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

    # Verify that the output degree sequence does not have a deficit compared to the reference degree sequence
    for cluster_id, output_degrees in degseq.items():
        if cluster_id in reference_degseq:
            reference_degrees = reference_degseq[cluster_id]
            if not verify_no_degree_deficit(output_degrees, reference_degrees):
                print(f"FAIL: Degree deficit found for cluster {cluster_id}.")
                return
        else:
            print(f"FAIL: Cluster {cluster_id} not found in reference degree sequence.")
            return
    print("PASS: No degree deficit found in the output degree sequence compared to the reference.")

if __name__ == "__main__":
    typer.run(main)