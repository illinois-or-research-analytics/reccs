"""
EC-SBM graph generator for a single reference network.

Wraps the EC-SBM pipeline (Tsitsulin et al.) which generates a synthetic
network that matches the empirical community structure, outlier structure,
and degree sequence. The pipeline runs four stages:

  1. Preprocessing: remove outlier nodes (clean_outlier.py)
  2. Stage 1:       generate synthetic clustered subnetwork (generate_clustered.py)
  3. Stage 2:       generate + combine outlier subnetwork (generate_outliers.py,
                    combine_networks.py)
  4. Stage 3:       degree correction + combine (correct_degree.py,
                    combine_networks.py)

Final output edge list: <tmpdir>/ecsbm+o+e/edge.tsv  → copied to -o <output.tsv>

The EC-SBM source lives in the ec-sbm/ subdirectory alongside this script.
Activate the ec-sbm conda environment before running:

    conda activate ec-sbm
    python3 gen_ec_sbm.py -f <edge_list.tsv> -c <clustering.tsv> -o <output.tsv>

Dependencies: the ec-sbm/ subpackage (numpy, scipy, networkx — see ec-sbm/env.yml)
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

EC_SBM_DIR = Path(__file__).parent / "ec-sbm" / "ec-sbm"


def run(cmd: list[str], cwd: Path) -> None:
    print("  $", " ".join(str(c) for c in cmd))
    subprocess.check_call([sys.executable] + [str(c) for c in cmd], cwd=str(cwd))


def parse_args():
    parser = argparse.ArgumentParser(description="EC-SBM graph generator")
    parser.add_argument('-f', '--edgelist',   required=True,
                        help='Input edge list TSV (source\\ttarget)')
    parser.add_argument('-c', '--clustering', required=True,
                        help='Input clustering TSV (node\\tcommunity)')
    parser.add_argument('-o', '--output',     required=True,
                        help='Output edge list TSV path')
    return parser.parse_args()


def main():
    args = parse_args()

    edgelist   = Path(args.edgelist).resolve()
    clustering = Path(args.clustering).resolve()
    output     = Path(args.output).resolve()

    if not edgelist.is_file():
        raise FileNotFoundError(f"Edge list not found: {edgelist}")
    if not clustering.is_file():
        raise FileNotFoundError(f"Clustering not found: {clustering}")
    if not EC_SBM_DIR.is_dir():
        raise FileNotFoundError(
            f"EC-SBM source not found at {EC_SBM_DIR}.\n"
            "Make sure the ec-sbm/ directory is present alongside this script."
        )

    with tempfile.TemporaryDirectory(prefix="ecsbm_") as tmpdir:
        tmp = Path(tmpdir)

        # ------------------------------------------------------------------
        # Preprocessing: remove outliers
        # ------------------------------------------------------------------
        print("Preprocessing: removing outliers ...")
        run([
            EC_SBM_DIR / "clean_outlier.py",
            "--input-edgelist",  edgelist,
            "--input-clustering", clustering,
            "--output-folder",   tmp / "emp_wo_o",
        ], cwd=EC_SBM_DIR)

        # ------------------------------------------------------------------
        # Stage 1: generate synthetic clustered subnetwork
        # ------------------------------------------------------------------
        print("Stage 1: generating synthetic clustered subnetwork ...")
        run([
            EC_SBM_DIR / "generate_clustered.py",
            "--input-edgelist",   tmp / "emp_wo_o" / "edge.tsv",
            "--input-clustering", tmp / "emp_wo_o" / "com.tsv",
            "--output-folder",    tmp / "ecsbm",
        ], cwd=EC_SBM_DIR)

        # ------------------------------------------------------------------
        # Stage 2: generate outlier subnetwork and combine
        # ------------------------------------------------------------------
        print("Stage 2: generating outlier subnetwork ...")
        run([
            EC_SBM_DIR / "generate_outliers.py",
            "--input-edgelist",   edgelist,
            "--input-clustering", clustering,
            "--output-folder",    tmp / "ecsbm+o",
        ], cwd=EC_SBM_DIR)

        print("Stage 2: combining clustered + outlier subnetworks ...")
        run([
            EC_SBM_DIR / "combine_networks.py",
            "--input-edgelist-1", tmp / "ecsbm"   / "edge.tsv",
            "--input-edgelist-2", tmp / "ecsbm+o" / "outlier_edge.tsv",
            "--input-clustering", tmp / "ecsbm"   / "com.tsv",
            "--output-folder",    tmp / "ecsbm+o",
        ], cwd=EC_SBM_DIR)

        # ------------------------------------------------------------------
        # Stage 3: degree correction and combine
        # ------------------------------------------------------------------
        print("Stage 3: degree correction ...")
        run([
            EC_SBM_DIR / "correct_degree.py",
            "--input-edgelist",  tmp / "ecsbm+o"   / "edge.tsv",
            "--ref-edgelist",    edgelist,
            "--ref-clustering",  clustering,
            "--output-folder",   tmp / "ecsbm+o+e",
        ], cwd=EC_SBM_DIR)

        print("Stage 3: combining with degree-corrected edges ...")
        run([
            EC_SBM_DIR / "combine_networks.py",
            "--input-edgelist-1", tmp / "ecsbm+o"   / "edge.tsv",
            "--input-edgelist-2", tmp / "ecsbm+o+e" / "degcorr_edge.tsv",
            "--input-clustering", tmp / "ecsbm+o"   / "com.tsv",
            "--output-folder",    tmp / "ecsbm+o+e",
        ], cwd=EC_SBM_DIR)

        # ------------------------------------------------------------------
        # Copy final output
        # ------------------------------------------------------------------
        final = tmp / "ecsbm+o+e" / "edge.tsv"
        if not final.is_file():
            raise RuntimeError(f"EC-SBM pipeline finished but output not found: {final}")

        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(final, output)

    print(f"Done. Output written to {output}")


if __name__ == '__main__':
    main()
