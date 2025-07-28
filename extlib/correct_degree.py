"""
Original author: The-Anh Vu Le
"""

import argparse
import logging
import heapq
import time
import csv
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set, Tuple

import pandas as pd
import numpy as np


class DegreeCorrector:
    """Handles degree sequence correction for networks."""
    
    def __init__(self):
        # Node ID mappings
        self.node_id_to_int: Dict[str, int] = {}
        self.node_int_to_id: Dict[int, str] = {}
        
        # Cluster ID mappings  
        self.cluster_id_to_int: Dict[str, int] = {}
        self.cluster_int_to_id: Dict[int, str] = {}
        
        # Node-cluster assignments
        self.node_to_cluster: Dict[int, int] = {}
        self.cluster_to_nodes: Dict[int, Set[int]] = defaultdict(set)
        
        # Network structure
        self.original_neighbors: Dict[int, Set[int]] = defaultdict(set)
        self.existing_neighbors: Dict[int, Set[int]] = defaultdict(set)
        
        # Degree information
        self.target_degrees: np.ndarray = None
        self.outliers: Set[int] = set()
        
    def add_node(self, node_id: str) -> int:
        """Add a node and return its integer ID."""
        if node_id not in self.node_id_to_int:
            node_int = len(self.node_id_to_int)
            self.node_id_to_int[node_id] = node_int
            self.node_int_to_id[node_int] = node_id
            return node_int
        return self.node_id_to_int[node_id]
    
    def add_cluster(self, cluster_id: str) -> int:
        """Add a cluster and return its integer ID."""
        if cluster_id not in self.cluster_id_to_int:
            cluster_int = len(self.cluster_id_to_int)
            self.cluster_id_to_int[cluster_id] = cluster_int
            self.cluster_int_to_id[cluster_int] = cluster_id
            return cluster_int
        return self.cluster_id_to_int[cluster_id]
    
    def load_clustering(self, clustering_file: Path) -> None:
        """Load node clustering information."""
        logging.info("Loading clustering information...")
        start = time.perf_counter()
        
        with open(clustering_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for node_id, cluster_id in reader:
                node_int = self.add_node(node_id)
                cluster_int = self.add_cluster(cluster_id)
                
                self.node_to_cluster[node_int] = cluster_int
                self.cluster_to_nodes[cluster_int].add(node_int)
        
        elapsed = time.perf_counter() - start
        logging.info(f"Loaded clustering: {elapsed:.2f}s")
    
    def load_original_network(self, edgelist_file: Path) -> None:
        """Load original network and compute target degrees."""
        logging.info("Loading original network...")
        start = time.perf_counter()
        
        with open(edgelist_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for src_id, tgt_id in reader:
                # Handle nodes not in clustering (outliers)
                src_int = self.add_node(src_id)
                tgt_int = self.add_node(tgt_id)
                
                if src_int not in self.node_to_cluster:
                    self.outliers.add(src_int)
                if tgt_int not in self.node_to_cluster:
                    self.outliers.add(tgt_int)
                
                self.original_neighbors[src_int].add(tgt_int)
                self.original_neighbors[tgt_int].add(src_int)
        
        # Create singleton clusters for outliers
        for outlier_int in self.outliers:
            cluster_int = len(self.cluster_id_to_int)
            cluster_id = str(cluster_int)  # Use integer as string ID
            self.add_cluster(cluster_id)
            
            self.node_to_cluster[outlier_int] = cluster_int
            self.cluster_to_nodes[cluster_int].add(outlier_int)
        
        # Compute target degrees
        num_nodes = len(self.node_int_to_id)
        self.target_degrees = np.zeros(num_nodes, dtype=int)
        for node_int, neighbors in self.original_neighbors.items():
            self.target_degrees[node_int] = len(neighbors)
        
        elapsed = time.perf_counter() - start
        logging.info(f"Loaded original network: {elapsed:.2f}s")
    
    def load_existing_network(self, edgelist_file: Path) -> None:
        """Load existing network and update remaining degrees."""
        logging.info("Loading existing network...")
        start = time.perf_counter()
        
        with open(edgelist_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for src_id, tgt_id in reader:
                src_int = self.node_id_to_int[src_id]
                tgt_int = self.node_id_to_int[tgt_id]
                
                # Skip if edge already exists
                if tgt_int in self.existing_neighbors[src_int]:
                    continue
                
                self.existing_neighbors[src_int].add(tgt_int)
                self.existing_neighbors[tgt_int].add(src_int)
                
                # Reduce target degrees
                self.target_degrees[src_int] = max(0, self.target_degrees[src_int] - 1)
                self.target_degrees[tgt_int] = max(0, self.target_degrees[tgt_int] - 1)
        
        elapsed = time.perf_counter() - start
        logging.info(f"Loaded existing network: {elapsed:.2f}s")
    
    def correct_degrees(self) -> Set[Tuple[int, int]]:
        """Perform degree correction and return new edges."""
        logging.info("Starting degree correction...")
        start = time.perf_counter()
        
        # Build available nodes with remaining degree
        available_nodes = {}
        for node_int in range(len(self.node_int_to_id)):
            remaining_degree = self.target_degrees[node_int]
            if remaining_degree > 0:
                available_nodes[node_int] = remaining_degree
        
        # Create max-heap for processing nodes by degree
        max_heap = [(-degree, node) for node, degree in available_nodes.items()]
        heapq.heapify(max_heap)
        
        new_edges = set()
        nodes_processed = 0
        edges_added = 0
        
        logging.info(f'Processing {len(max_heap)} nodes with remaining degree')
        
        while max_heap and available_nodes:
            # Get node with highest remaining degree
            _, current_node = heapq.heappop(max_heap)
            
            if current_node not in available_nodes:
                continue
            
            # Find available non-neighbors
            existing_neighbors = self.existing_neighbors.get(current_node, set())
            forbidden = existing_neighbors | {current_node}
            
            available_candidates = [
                node for node in available_nodes.keys() 
                if node not in forbidden
            ]
            
            # Add edges up to remaining degree or available candidates
            remaining_degree = available_nodes[current_node]
            max_edges = min(remaining_degree, len(available_candidates))
            
            for i in range(max_edges):
                if not available_candidates:
                    break
                    
                target_node = available_candidates.pop()
                
                # Add edge
                new_edges.add((current_node, target_node))
                
                # Update neighbors
                self.existing_neighbors[current_node].add(target_node)
                self.existing_neighbors[target_node].add(current_node)
                
                # Update available degrees
                available_nodes[target_node] -= 1
                if available_nodes[target_node] == 0:
                    del available_nodes[target_node]
                
                edges_added += 1
            
            # Remove processed node
            del available_nodes[current_node]
            nodes_processed += 1
            
            if nodes_processed % 1000 == 0:
                logging.info(f'Processed {nodes_processed} nodes, added {edges_added} edges')
        
        elapsed = time.perf_counter() - start
        logging.info(f'Degree correction completed: {elapsed:.2f}s')
        logging.info(f'Processed {nodes_processed} nodes, added {edges_added} edges')
        
        return new_edges
    
    def save_edges(self, edges: Set[Tuple[int, int]], output_file: Path) -> None:
        """Save edges to file."""
        logging.info("Saving corrected edges...")
        start = time.perf_counter()
        
        edge_data = [
            (self.node_int_to_id[src], self.node_int_to_id[tgt])
            for src, tgt in edges
        ]
        
        df = pd.DataFrame(edge_data, columns=['src_id', 'tgt_id'])
        df.to_csv(output_file, sep='\t', index=False, header=False)
        
        elapsed = time.perf_counter() - start
        logging.info(f"Saved {len(edges)} edges: {elapsed:.2f}s")


def setup_logging(output_dir: Path) -> None:
    """Setup logging configuration."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / 'degcorr_run.log'
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler()
        ]
    )


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Degree sequence correction for networks')
    parser.add_argument('--input-edgelist', type=str, required=True,
                       help='Path to existing network edgelist')
    parser.add_argument('--ref-edgelist', type=str, required=True,
                       help='Path to reference network edgelist')
    parser.add_argument('--ref-clustering', type=str, required=True,
                       help='Path to reference clustering file')
    parser.add_argument('--output-folder', type=str, required=True,
                       help='Output directory for results')
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()
    
    # Convert to Path objects
    existing_edgelist = Path(args.input_edgelist)
    reference_edgelist = Path(args.ref_edgelist)
    reference_clustering = Path(args.ref_clustering)
    output_dir = Path(args.output_folder)
    
    # Setup logging
    setup_logging(output_dir)
    
    # Log configuration
    logging.info('Degree Sequence Correction')
    logging.info(f'Reference network: {reference_edgelist}')
    logging.info(f'Reference clustering: {reference_clustering}')
    logging.info(f'Existing network: {existing_edgelist}')
    logging.info(f'Output folder: {output_dir}')
    
    # Initialize corrector
    corrector = DegreeCorrector()
    
    # Load data
    corrector.load_clustering(reference_clustering)
    corrector.load_original_network(reference_edgelist)
    corrector.load_existing_network(existing_edgelist)
    
    # Perform correction
    new_edges = corrector.correct_degrees()
    
    # Save results
    output_file = output_dir / 'degcorr_edge.tsv'
    corrector.save_edges(new_edges, output_file)
    
    logging.info('Degree correction completed successfully')


if __name__ == '__main__':
    main()
