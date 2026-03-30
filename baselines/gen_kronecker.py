"""
Kronecker graph generator for a single reference network.

KronFit: fits a 2x2 initiator matrix to the input graph via gradient ascent
on the Kronecker log-likelihood (Leskovec et al., JMLR 2010).
KronGen: samples a stochastic Kronecker graph from the fitted matrix.

Usage:
    python3 gen_kronecker.py -f <edge_list.tsv> -c <clustering.tsv> -o <output.tsv>

The clustering file is accepted for interface consistency with reccs/reccspp
but is not used by the Kronecker generator, which is community-agnostic.
"""

import argparse
import csv
import math
import random
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_graph(edge_path: str) -> tuple[set[tuple], set]:
    edges: set[tuple] = set()
    nodes: set = set()
    with open(edge_path) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2:
                continue
            u, v = row[0].strip(), row[1].strip()
            if u == v:
                continue
            edges.add((min(u, v), max(u, v)))
            nodes.add(u)
            nodes.add(v)
    return edges, nodes


def save_graph(edges: list[tuple], output_path: str) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for u, v in edges:
            writer.writerow([u, v])


# ---------------------------------------------------------------------------
# KronFit: gradient ascent on Kronecker log-likelihood
# ---------------------------------------------------------------------------

def compute_log_likelihood(theta: np.ndarray, adj: np.ndarray, k: int) -> float:
    """
    Approximate Kronecker log-likelihood.
    P_ij = product of theta[b_i(t), b_j(t)] over t=0..k-1,
    where b_i(t) is the t-th bit of node i's binary representation.
    Uses the full n x n probability matrix (feasible for n <= ~2048).
    """
    n = adj.shape[0]
    # Build probability matrix via iterated Kronecker product
    prob = theta.copy()
    for _ in range(k - 1):
        prob = np.kron(prob, theta)
    prob = np.clip(prob, 1e-10, 1 - 1e-10)
    # Log-likelihood: sum over all pairs
    ll = np.sum(adj * np.log(prob) + (1 - adj) * np.log(1 - prob))
    return float(ll)


def kronfit(edges: set[tuple], nodes: set, k: int,
            n_iters: int = 100, lr: float = 0.01,
            rng: random.Random = None) -> np.ndarray:
    """
    Fit a 2x2 Kronecker initiator matrix to the graph via gradient ascent.

    Uses a sampled gradient estimator: for each observed/non-observed edge,
    the gradient of log P_ij w.r.t. theta[a,b] is:
        grad = (A_ij - P_ij) / (theta[a,b] * P_ij / theta[a,b]) -- simplified
    We use the closed-form per-entry gradient of the log-likelihood.
    """
    if rng is None:
        rng = random.Random(42)

    # Map nodes to integer IDs 0..N-1
    node_list = sorted(nodes)
    node2id = {n: i for i, n in enumerate(node_list)}
    N = len(node_list)

    # Build adjacency matrix
    adj = np.zeros((N, N), dtype=np.float32)
    for u, v in edges:
        i, j = node2id[u], node2id[v]
        adj[i, j] = 1.0
        adj[j, i] = 1.0

    # Initialise theta with a community-like seed
    theta = np.array([[0.9, 0.5], [0.5, 0.2]], dtype=np.float64)

    # Gradient ascent with momentum
    velocity = np.zeros_like(theta)
    momentum = 0.9

    for iteration in range(n_iters):
        # Build full Kronecker probability matrix
        prob = theta.copy()
        for _ in range(k - 1):
            prob = np.kron(prob, theta)
        prob_full = np.clip(prob, 1e-10, 1 - 1e-10)

        # Truncate or pad to match N x N
        Nk = prob_full.shape[0]
        if Nk >= N:
            prob_used = prob_full[:N, :N]
            adj_used = adj
        else:
            prob_used = prob_full
            adj_used = adj[:Nk, :Nk]

        # Compute gradient of log-likelihood w.r.t. each theta[a,b]
        # d LL / d theta[a,b] = sum_{ij} (A_ij/P_ij - (1-A_ij)/(1-P_ij)) * dP_ij/d_theta[a,b]
        # dP_ij / d theta[a,b] = P_ij * (count of (a,b) at level t for pair ij) / theta[a,b]
        # We approximate this cheaply: for each (i,j), decompose into k bit-pairs and accumulate.
        grad = np.zeros((2, 2), dtype=np.float64)
        n_used = prob_used.shape[0]

        residual = (adj_used - prob_used) / (prob_used * (1 - prob_used))  # (n_used, n_used)

        for t in range(k):
            stride = 2 ** t
            for a in range(2):
                for b in range(2):
                    # Mask of (i,j) pairs where bit t of i is a and bit t of j is b
                    # bit t of integer x = (x // stride) % 2
                    i_idx = np.array([i for i in range(n_used)
                                      if (i // stride) % 2 == a])
                    j_idx = np.array([j for j in range(n_used)
                                      if (j // stride) % 2 == b])
                    if len(i_idx) == 0 or len(j_idx) == 0:
                        continue
                    # Contribution: P_ij / theta[a,b] * residual[i,j]
                    sub_prob = prob_used[np.ix_(i_idx, j_idx)]
                    sub_res = residual[np.ix_(i_idx, j_idx)]
                    if theta[a, b] > 1e-12:
                        grad[a, b] += np.sum(sub_prob / theta[a, b] * sub_res)

        # Update with momentum, project to (0,1)
        velocity = momentum * velocity + lr * grad
        theta = np.clip(theta + velocity, 0.01, 0.99)

        if (iteration + 1) % 10 == 0:
            ll = compute_log_likelihood(theta, adj_used, k)
            print(f"  iter {iteration+1:3d}/{n_iters}  ll={ll:.2f}  theta={theta.flatten().round(4)}")

    return theta


# ---------------------------------------------------------------------------
# KronGen: sample a stochastic Kronecker graph
# ---------------------------------------------------------------------------

def krongen(theta: np.ndarray, k: int, node_list: list,
            directed: bool = False) -> list[tuple]:
    """
    Sample a stochastic Kronecker graph with 2^k nodes using the fast
    edge-by-edge sampler (Leskovec & Faloutsos, KDD 2010).

    Returns edges as (node_label, node_label) pairs using node_list labels,
    truncated / mapped back to the original node IDs.
    """
    n = 2 ** k
    n_orig = len(node_list)

    # Expected number of edges ≈ (sum of theta)^k
    expected_edges = int(round(np.sum(theta) ** k))

    edges = []
    seen = set()

    for _ in range(expected_edges * 2):  # oversample then deduplicate
        u, v = 0, 0
        for _ in range(k):
            # Flatten theta to choose one of 4 cells proportionally
            flat = theta.flatten()
            flat = flat / flat.sum()
            cell = np.random.choice(4, p=flat)
            row, col = divmod(cell, 2)
            u = u * 2 + row
            v = v * 2 + col

        if u == v:
            continue
        if not directed:
            u, v = min(u, v), max(u, v)
        if (u, v) in seen:
            continue
        seen.add((u, v))

        # Map synthetic integer IDs back to original node labels (wrap if n > n_orig)
        u_label = node_list[u % n_orig]
        v_label = node_list[v % n_orig]
        if u_label != v_label:
            edges.append((u_label, v_label))

        if len(edges) >= expected_edges:
            break

    return edges


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Kronecker graph generator")
    parser.add_argument('-f', '--edgelist', required=True,
                        help='Input edge list TSV (source\\ttarget)')
    parser.add_argument('-c', '--clustering', required=True,
                        help='Clustering TSV (accepted for interface consistency, unused)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output edge list TSV path')
    parser.add_argument('--iters', type=int, default=100,
                        help='KronFit gradient ascent iterations (default: 100)')
    parser.add_argument('--lr', type=float, default=0.005,
                        help='KronFit learning rate (default: 0.005)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    return parser.parse_args()


def main():
    args = parse_args()
    np.random.seed(args.seed)
    random.seed(args.seed)

    print(f"Loading graph from {args.edgelist} ...")
    edges, nodes = load_graph(args.edgelist)
    N = len(nodes)
    M = len(edges)
    print(f"  {N} nodes, {M} edges")

    if N == 0:
        print("Empty graph — writing empty output.")
        save_graph([], args.output)
        return

    # Choose k such that 2^k >= N
    k = max(1, math.ceil(math.log2(N)))
    print(f"  Using k={k} (2^k={2**k} nodes in Kronecker graph)")

    print(f"Running KronFit ({args.iters} iterations) ...")
    node_list = sorted(nodes)
    theta = kronfit(edges, nodes, k,
                    n_iters=args.iters, lr=args.lr,
                    rng=random.Random(args.seed))
    print(f"Fitted theta:\n{theta}")

    print("Sampling Kronecker graph ...")
    syn_edges = krongen(theta, k, node_list)
    print(f"  Generated {len(syn_edges)} edges")

    print(f"Saving to {args.output} ...")
    save_graph(syn_edges, args.output)
    print("Done.")


if __name__ == '__main__':
    main()
