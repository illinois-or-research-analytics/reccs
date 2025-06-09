import json
import typer
import pandas as pd


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

    print(stats_df.shape)
    print(reference_stats_df.shape)
    print(len(degseq))
    print(len(reference_degseq))
    

if __name__ == "__main__":
    typer.run(main)