import pandas as pd
import typer
from tqdm import tqdm
import multiprocessing as mp
import os
from functools import partial


def process_cluster(edgelist, nodes, cluster_id, output_dir):
    """Process a single cluster and save its edges to a file."""
    # Filter edges where both source and target are in the current cluster
    cluster_df = edgelist[
        (edgelist["source"].isin(nodes)) & (edgelist["target"].isin(nodes))
    ]
    
    # Save to file
    output_file = f"{output_dir}/cluster_{cluster_id}.tsv"
    cluster_df.to_csv(output_file, sep="\t", index=False, header=False)
    
    return cluster_id


def main(
    input_file: str = typer.Option(..., "--edgelist", "-e", help="Path to the input TSV edgelist"),
    cluster_file: str = typer.Option(..., "--clustering", "-c", help="Column name for the cluster TSV"),
    output_dir: str = typer.Option(..., "--output-directory", "-o", help="Directory to save the output TSV files"),
    n_processes: int = typer.Option(None, "--processes", "-p", help="Number of processes to use (defaults to CPU count)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose output"),
):
    """
    Splits a TSV edgelist into separate cluster subgraphs based on a specified cluster column.
    Uses multiprocessing to speed up the filtering and saving operations.
    """
    # Function to print only if verbose is enabled
    def vprint(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
            
    # Create tqdm iterator that only shows progress if verbose is enabled
    def vtqdm(*args, **kwargs):
        if verbose:
            return tqdm(*args, **kwargs)
        else:
            return args[0]  # Return the iterable without progress tracking
    
    vprint("Reading input files...")
    edgelist = pd.read_csv(input_file, sep="\t", header=None)
    edgelist.columns = ["source", "target"]

    clustering = pd.read_csv(cluster_file, sep="\t", header=None)
    clustering.columns = ["node", "cluster"]

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Get a mapping of cluster IDs to arrays of nodes
    vprint("Creating cluster mappings...")
    cluster_mapping = clustering.groupby("cluster")["node"].apply(list).to_dict()
    
    # Determine number of processes
    if n_processes is None:
        n_processes = mp.cpu_count()
    vprint(f"Using {n_processes} processes for parallel execution")
    
    # Process clusters in parallel
    vprint("Processing and saving clusters in parallel...")
    process_func = partial(process_cluster, edgelist, output_dir=output_dir)
    
    with mp.Pool(processes=n_processes) as pool:
        # Map the process_cluster function to each (cluster_id, nodes) pair
        items = [(nodes, cluster_id) for cluster_id, nodes in cluster_mapping.items()]
        
        # Process in parallel with progress bar (only if verbose)
        results = list(vtqdm(
            pool.starmap(process_func, items),
            total=len(items),
            desc="Processing clusters"
        ))
    
    vprint(f"Successfully processed {len(results)} clusters")


if __name__ == "__main__":
    typer.run(main)
    