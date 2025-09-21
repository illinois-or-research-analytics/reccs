import json
import typer
import pandas as pd
import matplotlib.pyplot as plt # type: ignore

def main(
    edgelist_output: str = typer.Option(..., "--edgelist", "-e", help="Output edgelist filename"),
    sbm_edgelist_output: str = typer.Option(..., "--sbm_edgelist", "-se", help="SBM edgelist filename"),
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
    # Load the edgelist output
    try:
        edgelist_df = pd.read_csv(edgelist_output, header=None, sep='\t')
        print(edgelist_df.head())
        if edgelist_df.shape[1] != 2:
            raise ValueError("Edgelist output must have exactly two columns for source and target nodes.")
    except Exception as e:
        print(f"FAIL: Could not load edgelist output file '{edgelist_output}'. Error: {e}")
        return
    
    print(f"PASS: Edgelist output file '{edgelist_output}' loaded successfully.")

    # Check if the edgelist contains any self-loops
    if (edgelist_df[0] == edgelist_df[1]).any():
        print("FAIL: Edgelist contains self-loops.")
    else:
        print("PASS: No self-loops found in the edgelist.")

    # Check if the edgelist contains any duplicate edges
    if edgelist_df.duplicated().any():
        print("FAIL: Edgelist contains duplicate edges.")
    else:
        print("PASS: No duplicate edges found in the edgelist.")

    # Do the same checks for the SBM edgelist output
    try:
        sbm_edgelist_df = pd.read_csv(sbm_edgelist_output, header=None, sep='\t')
        print(sbm_edgelist_df.head())
        if sbm_edgelist_df.shape[1] != 2:
            raise ValueError("SBM edgelist output must have exactly two columns for source and target nodes.")
    except Exception as e:
        print(f"FAIL: Could not load SBM edgelist output file '{sbm_edgelist_output}'. Error: {e}")
        return
    print(f"PASS: SBM edgelist output file '{sbm_edgelist_output}' loaded successfully.")
    # Check if the SBM edgelist contains any self-loops
    if (sbm_edgelist_df[0] == sbm_edgelist_df[1]).any():
        print("FAIL: SBM edgelist contains self-loops.")
    else:
        print("PASS: No self-loops found in the SBM edgelist.")

    # Check if the SBM edgelist contains any duplicate edges
    if sbm_edgelist_df.duplicated().any():
        print("FAIL: SBM edgelist contains duplicate edges.")
    else:
        print("PASS: No duplicate edges found in the SBM edgelist.")

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

    # Check if any disconnected clusters exist in the stats connectivity mapping
    disconnected_clusters = [cluster for cluster, connectivity in stats_connectivity.items() if connectivity == 0]
    if disconnected_clusters:
        print(f"FAIL: Disconnected clusters found in stats connectivity mapping: {disconnected_clusters}")
    else:
        print("PASS: No disconnected clusters found in stats connectivity mapping!")

    # Check if connectivity matches between RECCS and reference stats
    connectivity_mismatch = False
    for cluster_id, connectivity in stats_connectivity.items():
        cluster_id = int(cluster_id)
        if cluster_id in reference_connectivity:
            reference_connectivity_value = reference_connectivity[cluster_id]
            if connectivity < reference_connectivity_value:
                print(f"FAIL: RECCS connectivity {connectivity} does not match reference connectivity {reference_connectivity_value}.")
                connectivity_mismatch = True
        else:
            print(f"FAIL: Cluster {cluster_id} not found in reference stats connectivity mapping.")
    if not connectivity_mismatch:
        print("PASS: Connectivity matches between RECCS and reference stats for all clusters!")


    # Plot the histograms on the same figure
    plt.figure(figsize=(10, 6))

    ref_degrees = reference_degseq
    reccs_degrees = degseq
    
    # Create box plot comparing SBM and RECCS distances
    data_to_plot = [ref_degrees, reccs_degrees]
    labels = ['Reference Degrees', 'RECCS Degrees']

    plt.boxplot(data_to_plot, labels=labels)
    plt.title('Comparison of Degree Sequences')
    plt.ylabel('Euclidean Distance')
    plt.grid(True, alpha=0.3)
    plt.savefig(plot_output_path)

if __name__ == "__main__":
    typer.run(main)
    