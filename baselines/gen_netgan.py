"""
NetGAN graph generator for a single reference network.

Implements the NetGAN architecture from:
    Bojchevski et al., "NetGAN: Generating Graphs via Random Walks", ICML 2018.

The generator learns to produce random walks that are indistinguishable from
real random walks on the input graph, then reconstructs a graph by scoring
edges with walk co-occurrence frequency.

Usage:
    python3 gen_netgan.py -f <edge_list.tsv> -c <clustering.tsv> -o <output.tsv>

The clustering file is accepted for interface consistency with reccs/reccspp
but is not used by NetGAN, which is community-agnostic.

Dependencies (installed automatically if missing):
    torch
"""

# Dependency: pip install torch
# e.g.: conda activate gt2 && pip install torch

import argparse
import csv
import random
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_graph(edge_path: str) -> tuple[list[tuple[int, int]], dict, list]:
    raw_edges: set[tuple] = set()
    raw_nodes: set = set()
    with open(edge_path) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2:
                continue
            u, v = row[0].strip(), row[1].strip()
            if u == v:
                continue
            raw_edges.add((min(u, v), max(u, v)))
            raw_nodes.add(u)
            raw_nodes.add(v)
    node_list = sorted(raw_nodes)
    node2id = {n: i for i, n in enumerate(node_list)}
    edges = [(node2id[u], node2id[v]) for u, v in raw_edges]
    return edges, node2id, node_list


def save_graph(edges: list[tuple], node_list: list, output_path: str) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for u, v in edges:
            writer.writerow([node_list[u], node_list[v]])


# ---------------------------------------------------------------------------
# Random walk sampling
# ---------------------------------------------------------------------------

def build_adjacency(edges: list[tuple[int, int]], N: int) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(N)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def sample_random_walks(adj: list[list[int]], N: int,
                        rw_len: int, n_walks: int,
                        rng: random.Random) -> np.ndarray:
    """Sample random walks of length rw_len from the graph."""
    walks = []
    nodes_with_edges = [i for i in range(N) if len(adj[i]) > 0]
    if not nodes_with_edges:
        return np.zeros((n_walks, rw_len), dtype=np.int64)
    for _ in range(n_walks):
        start = rng.choice(nodes_with_edges)
        walk = [start]
        cur = start
        for _ in range(rw_len - 1):
            nbrs = adj[cur]
            if nbrs:
                cur = rng.choice(nbrs)
            else:
                cur = rng.choice(nodes_with_edges)
            walk.append(cur)
        walks.append(walk)
    return np.array(walks, dtype=np.int64)


# ---------------------------------------------------------------------------
# NetGAN model: LSTM-based generator and discriminator (WGAN-GP)
# ---------------------------------------------------------------------------

class Generator(nn.Module):
    """LSTM generator that produces sequences of node IDs as random walks."""

    def __init__(self, n_nodes: int, rw_len: int,
                 z_dim: int = 16, hidden_dim: int = 64):
        super().__init__()
        self.n_nodes = n_nodes
        self.rw_len = rw_len
        self.z_dim = z_dim
        self.hidden_dim = hidden_dim

        self.z_to_hidden = nn.Linear(z_dim, hidden_dim)
        self.lstm = nn.LSTM(n_nodes, hidden_dim, batch_first=True)
        self.output = nn.Linear(hidden_dim, n_nodes)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """
        Args:
            z: (batch, z_dim) noise vector
        Returns:
            logits: (batch, rw_len, n_nodes) soft walk
        """
        batch = z.size(0)
        h0 = self.z_to_hidden(z).unsqueeze(0)          # (1, batch, hidden)
        c0 = torch.zeros_like(h0)

        # Start token: uniform over all nodes
        inp = torch.full((batch, 1, self.n_nodes),
                         1.0 / self.n_nodes, device=z.device)

        logits_list = []
        h, c = h0, c0
        for _ in range(self.rw_len):
            out, (h, c) = self.lstm(inp, (h, c))        # out: (batch,1,hidden)
            logit = self.output(out)                     # (batch,1,n_nodes)
            logits_list.append(logit)
            inp = torch.softmax(logit, dim=-1)
        return torch.cat(logits_list, dim=1)             # (batch, rw_len, n_nodes)

    def sample_walks(self, z: torch.Tensor) -> torch.Tensor:
        """Sample discrete walks via argmax. Returns (batch, rw_len) int tensor."""
        with torch.no_grad():
            logits = self.forward(z)
        return logits.argmax(dim=-1)


class Discriminator(nn.Module):
    """LSTM discriminator that scores a sequence of one-hot node visits."""

    def __init__(self, n_nodes: int, rw_len: int, hidden_dim: int = 64):
        super().__init__()
        self.lstm = nn.LSTM(n_nodes, hidden_dim, batch_first=True)
        self.score = nn.Linear(hidden_dim, 1)

    def forward(self, walks_onehot: torch.Tensor) -> torch.Tensor:
        """
        Args:
            walks_onehot: (batch, rw_len, n_nodes) float tensor
        Returns:
            scores: (batch,)
        """
        _, (h, _) = self.lstm(walks_onehot)
        return self.score(h.squeeze(0)).squeeze(-1)


# ---------------------------------------------------------------------------
# WGAN-GP training
# ---------------------------------------------------------------------------

def gradient_penalty(disc: Discriminator,
                     real: torch.Tensor,
                     fake: torch.Tensor,
                     device: torch.device) -> torch.Tensor:
    batch = real.size(0)
    alpha = torch.rand(batch, 1, 1, device=device)
    interp = (alpha * real + (1 - alpha) * fake).requires_grad_(True)
    d_interp = disc(interp)
    grads = torch.autograd.grad(
        outputs=d_interp, inputs=interp,
        grad_outputs=torch.ones_like(d_interp),
        create_graph=True, retain_graph=True)[0]
    gp = ((grads.norm(2, dim=(1, 2)) - 1) ** 2).mean()
    return gp


def train_netgan(adj: list[list[int]], N: int,
                 rw_len: int = 16,
                 z_dim: int = 16,
                 hidden_dim: int = 64,
                 batch_size: int = 128,
                 n_iter: int = 5000,
                 n_critic: int = 5,
                 lr: float = 1e-4,
                 lambda_gp: float = 10.0,
                 device: torch.device = torch.device('cpu'),
                 rng: random.Random = None) -> tuple[Generator, int]:

    if rng is None:
        rng = random.Random(42)

    gen = Generator(N, rw_len, z_dim, hidden_dim).to(device)
    disc = Discriminator(N, rw_len, hidden_dim).to(device)

    opt_g = torch.optim.Adam(gen.parameters(), lr=lr, betas=(0.5, 0.9))
    opt_d = torch.optim.Adam(disc.parameters(), lr=lr, betas=(0.5, 0.9))

    def real_walks_batch(bs: int) -> torch.Tensor:
        walks = sample_random_walks(adj, N, rw_len, bs, rng)
        onehot = np.zeros((bs, rw_len, N), dtype=np.float32)
        for b in range(bs):
            for t in range(rw_len):
                onehot[b, t, walks[b, t]] = 1.0
        return torch.tensor(onehot, device=device)

    for iteration in range(1, n_iter + 1):
        # --- Critic steps ---
        for _ in range(n_critic):
            z = torch.randn(batch_size, z_dim, device=device)
            real = real_walks_batch(batch_size)
            with torch.no_grad():
                fake_logits = gen(z)
                fake = torch.softmax(fake_logits, dim=-1)

            d_real = disc(real).mean()
            d_fake = disc(fake).mean()
            gp = gradient_penalty(disc, real, fake.detach(), device)
            loss_d = d_fake - d_real + lambda_gp * gp

            opt_d.zero_grad()
            loss_d.backward()
            opt_d.step()

        # --- Generator step ---
        z = torch.randn(batch_size, z_dim, device=device)
        fake_logits = gen(z)
        fake = torch.softmax(fake_logits, dim=-1)
        loss_g = -disc(fake).mean()

        opt_g.zero_grad()
        loss_g.backward()
        opt_g.step()

        if iteration % 500 == 0:
            print(f"  iter {iteration:5d}/{n_iter}  "
                  f"D={loss_d.item():.3f}  G={loss_g.item():.3f}")

    return gen, N


# ---------------------------------------------------------------------------
# Graph reconstruction from sampled walks
# ---------------------------------------------------------------------------

def reconstruct_graph(gen: Generator,
                      n_nodes: int,
                      n_samples: int,
                      z_dim: int,
                      target_edges: int,
                      device: torch.device) -> list[tuple[int, int]]:
    """
    Score each potential edge by co-occurrence frequency in generated walks,
    then threshold to match the target edge count.
    """
    scores = np.zeros((n_nodes, n_nodes), dtype=np.float32)
    batch = 256

    n_batches = max(1, n_samples // batch)
    for _ in range(n_batches):
        z = torch.randn(batch, z_dim, device=device)
        walks = gen.sample_walks(z).cpu().numpy()  # (batch, rw_len)
        for walk in walks:
            for t in range(len(walk) - 1):
                u, v = int(walk[t]), int(walk[t + 1])
                if u != v:
                    scores[u, v] += 1.0
                    scores[v, u] += 1.0

    # Symmetrise and take upper triangle
    scores = (scores + scores.T) / 2.0
    np.fill_diagonal(scores, 0.0)

    # Pick top-k edges
    rows, cols = np.triu_indices(n_nodes, k=1)
    edge_scores = scores[rows, cols]
    if target_edges >= len(edge_scores):
        mask = edge_scores > 0
    else:
        threshold = np.partition(edge_scores, -target_edges)[-target_edges]
        mask = edge_scores >= threshold

    edges = [(int(rows[i]), int(cols[i]))
             for i in np.where(mask)[0]
             if rows[i] != cols[i]]
    return edges


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="NetGAN graph generator")
    parser.add_argument('-f', '--edgelist', required=True,
                        help='Input edge list TSV (source\\ttarget)')
    parser.add_argument('-c', '--clustering', required=True,
                        help='Clustering TSV (accepted for interface consistency, unused)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output edge list TSV path')
    parser.add_argument('--rw-len', type=int, default=16,
                        help='Random walk length (default: 16)')
    parser.add_argument('--z-dim', type=int, default=16,
                        help='Noise vector dimension (default: 16)')
    parser.add_argument('--hidden-dim', type=int, default=64,
                        help='LSTM hidden size (default: 64)')
    parser.add_argument('--batch-size', type=int, default=128,
                        help='Training batch size (default: 128)')
    parser.add_argument('--n-iter', type=int, default=5000,
                        help='Training iterations (default: 5000)')
    parser.add_argument('--n-critic', type=int, default=5,
                        help='Discriminator steps per generator step (default: 5)')
    parser.add_argument('--lr', type=float, default=1e-4,
                        help='Learning rate (default: 1e-4)')
    parser.add_argument('--n-samples', type=int, default=5000,
                        help='Walks to sample for graph reconstruction (default: 5000)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    return parser.parse_args()


def main():
    args = parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    rng = random.Random(args.seed)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    print(f"Loading graph from {args.edgelist} ...")
    edges, node2id, node_list = load_graph(args.edgelist)
    N = len(node_list)
    M = len(edges)
    print(f"  {N} nodes, {M} edges")

    if N == 0 or M == 0:
        print("Empty graph — writing empty output.")
        save_graph([], node_list, args.output)
        return

    adj = build_adjacency(edges, N)

    print(f"Training NetGAN ({args.n_iter} iterations, rw_len={args.rw_len}) ...")
    gen, _ = train_netgan(
        adj=adj, N=N,
        rw_len=args.rw_len,
        z_dim=args.z_dim,
        hidden_dim=args.hidden_dim,
        batch_size=args.batch_size,
        n_iter=args.n_iter,
        n_critic=args.n_critic,
        lr=args.lr,
        device=device,
        rng=rng,
    )

    print(f"Reconstructing graph from {args.n_samples} sampled walks ...")
    syn_edges = reconstruct_graph(
        gen=gen, n_nodes=N,
        n_samples=args.n_samples,
        z_dim=args.z_dim,
        target_edges=M,
        device=device,
    )
    print(f"  Generated {len(syn_edges)} edges")

    print(f"Saving to {args.output} ...")
    save_graph(syn_edges, node_list, args.output)
    print("Done.")


if __name__ == '__main__':
    main()
