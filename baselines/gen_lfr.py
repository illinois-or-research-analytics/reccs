"""
LFR benchmark graph generator for a single reference network.

Uses NetworKit's LFRGenerator fitted to the input graph. The empirical
clustering is passed directly via setPartition, so the community structure
matches the ground-truth clustering supplied — consistent with how
reccs/reccspp take the clustering as input.

Degree sequence and mixing parameter mu are derived from the input graph.

Usage:
    python3 gen_lfr.py -f <edge_list.tsv> -c <clustering.tsv> -o <output.tsv>

Dependency: pip install networkit
"""

import argparse
import csv
from pathlib import Path

import networkit as nk


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

    return G, node_list, node2id


def load_partition(clustering_path: str, node2id: dict, N: int) -> nk.Partition:
    partition = nk.Partition(N)
    partition.allToSingletons()

    cluster_label2id: dict[str, int] = {}
    next_cluster_id = 0

    with open(clustering_path) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2:
                continue
            node_str, cluster_str = row[0].strip(), row[1].strip()
            if node_str not in node2id:
                continue
            node_int = node2id[node_str]
            if cluster_str not in cluster_label2id:
                cluster_label2id[cluster_str] = next_cluster_id
                next_cluster_id += 1
            partition.moveToSubset(cluster_label2id[cluster_str], node_int)

    return partition


def save_graph(G: nk.Graph, node_list: list, n_orig: int,
               output_path: str) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for u, v in G.iterEdges():
            u_label = node_list[u % n_orig]
            v_label = node_list[v % n_orig]
            if u_label != v_label:
                writer.writerow([u_label, v_label])


# ---------------------------------------------------------------------------
# Compute mu from empirical graph + partition
# ---------------------------------------------------------------------------

def compute_per_node_mu(G: nk.Graph, partition: nk.Partition) -> list[float]:
    """
    Per-node mu clamped so that internal degree <= community_size - 1.
    For node u with degree d in community of size s:
      max internal degree = min(d, s-1)
      min mu = 1 - (s-1)/d  (when all possible intra edges are used)
    We use the empirical mu per node but floor it at that minimum.
    """
    # Community sizes
    community_sizes: dict[int, int] = {}
    for u in G.iterNodes():
        c = partition.subsetOf(u)
        community_sizes[c] = community_sizes.get(c, 0) + 1

    mu_per_node = []
    for u in G.iterNodes():
        deg = G.degree(u)
        if deg == 0:
            mu_per_node.append(0.5)
            continue
        inter = sum(1 for v in G.iterNeighbors(u)
                    if partition.subsetOf(v) != partition.subsetOf(u))
        empirical_mu = inter / deg
        # Minimum feasible mu: internal degree can be at most community_size - 1
        c_size = community_sizes.get(partition.subsetOf(u), 1)
        max_internal = min(deg, c_size - 1)
        min_mu = 1.0 - max_internal / deg if deg > 0 else 0.0
        mu_per_node.append(max(empirical_mu, min_mu))

    return mu_per_node


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="LFR graph generator via NetworKit")
    parser.add_argument('-f', '--edgelist', required=True,
                        help='Input edge list TSV (source\\ttarget)')
    parser.add_argument('-c', '--clustering', required=True,
                        help='Clustering TSV (node\\tcommunity)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output edge list TSV path')
    parser.add_argument('--scale', type=int, default=1,
                        help='Node count scale factor (default: 1)')
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading graph from {args.edgelist} ...")
    G, node_list, node2id = load_graph(args.edgelist)
    N = G.numberOfNodes()
    M = G.numberOfEdges()
    print(f"  {N} nodes, {M} edges")

    if N == 0 or M == 0:
        print("Empty graph — writing empty output.")
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text("")
        return

    print(f"Loading clustering from {args.clustering} ...")
    partition = load_partition(args.clustering, node2id, N)
    n_communities = partition.numberOfSubsets()
    print(f"  {n_communities} communities")

    mu_per_node = compute_per_node_mu(G, partition)
    avg_mu = sum(mu_per_node) / len(mu_per_node)
    print(f"  Empirical mu (per-node avg): {avg_mu:.4f}")

    print("Fitting LFR generator to input graph ...")
    gen = nk.generators.LFRGenerator.fit(G, scale=args.scale)

    # Override community structure with ground-truth clustering and per-node mu
    gen.setPartition(partition)
    gen.setMu(mu_per_node)

    print("Generating synthetic graph ...")
    gen.run()
    syn_G = gen.getGraph()
    print(f"  Generated {syn_G.numberOfNodes()} nodes, {syn_G.numberOfEdges()} edges")

    print(f"Saving to {args.output} ...")
    save_graph(syn_G, node_list, N, args.output)
    print("Done.")


if __name__ == '__main__':
    main()
