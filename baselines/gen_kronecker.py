"""
Kronecker graph generator for a single reference network.

Uses NetworKit's RmatGenerator with the SNAP KronFit binary for fitting.
RmatGenerator implements the stochastic Kronecker / R-MAT model
(Chakrabarti et al., SDM 2004; Leskovec et al., JMLR 2010).

NetworKit's fit() method shells out to the SNAP kronfit binary. The SNAP
source is included as a git submodule at baselines/snap_kronfit/. Build it
once before running this script:

    cd baselines/snap_kronfit/examples/kronfit && make -j4

Usage:
    python3 gen_kronecker.py -f <edge_list.tsv> -c <clustering.tsv> -o <output.tsv>

The clustering file is accepted for interface consistency with reccs/reccspp
but is not used by the Kronecker generator, which is community-agnostic.

Dependencies: pip install networkit
"""

import argparse
import csv
import math
from pathlib import Path
from typing import Optional

import networkit as nk


# ---------------------------------------------------------------------------
# Locate SNAP kronfit binary
# ---------------------------------------------------------------------------

KRONFIT_BUILD_DIR = Path(__file__).parent / "snap_kronfit"


def ensure_kronfit(kronfit_path: Optional[str]) -> str:
    """Return path to a working kronfit binary."""
    if kronfit_path and Path(kronfit_path).is_file():
        return kronfit_path

    binary = KRONFIT_BUILD_DIR / "examples" / "kronfit" / "kronfit"
    if binary.is_file():
        return str(binary)

    raise RuntimeError(
        f"kronfit binary not found at {binary}.\n"
        "Build it with:\n"
        "  cd baselines/snap_kronfit/examples/kronfit && make -j4"
    )


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_graph(edge_path: str) -> tuple[nk.Graph, list]:
    raw_nodes: set = set()
    raw_edges: set[tuple] = set()
    with open(edge_path) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2:
                continue
            u, v = row[0].strip(), row[1].strip()
            if u == v:
                continue
            raw_nodes.add(u)
            raw_nodes.add(v)
            raw_edges.add((min(u, v), max(u, v)))

    node_list = sorted(raw_nodes)
    node2id = {n: i for i, n in enumerate(node_list)}
    N = len(node_list)

    G = nk.Graph(N, directed=False)
    for u, v in raw_edges:
        G.addEdge(node2id[u], node2id[v])

    return G, node_list


def save_graph(G: nk.Graph, node_list: list, n_orig: int,
               output_path: str) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for u, v in G.iterEdges():
            # Map synthetic integer IDs back to original labels
            u_label = node_list[u % n_orig]
            v_label = node_list[v % n_orig]
            if u_label != v_label:
                writer.writerow([u_label, v_label])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Kronecker (R-MAT) graph generator via NetworKit")
    parser.add_argument('-f', '--edgelist', required=True,
                        help='Input edge list TSV (source\\ttarget)')
    parser.add_argument('-c', '--clustering', required=True,
                        help='Clustering TSV (accepted for interface consistency, unused)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output edge list TSV path')
    parser.add_argument('--scale', type=int, default=1,
                        help='Node count scale factor passed to RmatGenerator.fit (default: 1)')
    parser.add_argument('--kronfit-iters', type=int, default=50,
                        help='KronFit iterations (default: 50)')
    parser.add_argument('--kronfit-path', type=str, default=None,
                        help='Path to kronfit binary (built automatically if not provided)')
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading graph from {args.edgelist} ...")
    G, node_list = load_graph(args.edgelist)
    N = G.numberOfNodes()
    M = G.numberOfEdges()
    print(f"  {N} nodes, {M} edges")

    if N == 0 or M == 0:
        print("Empty graph — writing empty output.")
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text("")
        return

    kronfit_bin = ensure_kronfit(args.kronfit_path)
    nk.generators.RmatGenerator.setPaths(kronfit_bin)

    print(f"Fitting R-MAT/Kronecker model (kronfit=True, iters={args.kronfit_iters}) ...")
    gen = nk.generators.RmatGenerator.fit(
        G,
        scale=args.scale,
        kronfit=True,
        iterations=args.kronfit_iters,
    )

    print("Generating synthetic graph ...")
    syn_G = gen.generate()
    print(f"  Generated {syn_G.numberOfNodes()} nodes, {syn_G.numberOfEdges()} edges")

    print(f"Saving to {args.output} ...")
    save_graph(syn_G, node_list, N, args.output)
    print("Done.")


if __name__ == '__main__':
    main()
