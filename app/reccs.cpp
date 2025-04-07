#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "io/graph_io.h"
#include "utils/subgraph_extractor.h"
#include "algorithm/sbm.h"


int main(int argc, char* argv[]) {
    // Read command line arguments
    bool verbose = false;
    std::string edgelist_file;
    std::string clustering_file;
    for (int i = 1; i < argc; ++i) {
        // Read verbose flag
        if (std::string(argv[i]) == "-v" || std::string(argv[i]) == "--verbose") {
            verbose = true;
        }

        // Read clustering file
        if (std::string(argv[i]) == "-c" || std::string(argv[i]) == "--clustering") {
            if (i + 1 < argc) {
                clustering_file = argv[++i];
            } else {
                std::cerr << "Error: No clustering file provided." << std::endl;
                return 1;
            }
        }

        // Read edgelist file
        if (std::string(argv[i]) == "-e" || std::string(argv[i]) == "--edgelist") {
            if (i + 1 < argc) {
                edgelist_file = argv[++i];
            } else {
                std::cerr << "Error: No edgelist file provided." << std::endl;
                return 1;
            }
        }
    }

    // Read the edgelist and clustering files
    if (verbose) {
        std::cout << "Reading graph and clustering..." << std::endl;
    }

    // Use graph_io to read the graph and clustering
    graph_io io;
    
    // Read the edgelist
    Graph<int> graph = io.load_graph_from_tsv(edgelist_file);
    if (!graph.node_count()) {
        std::cerr << "Error opening edgelist file: " << edgelist_file << std::endl;
        return 1;
    }
    
    // Read clustering
    auto clustering = io.load_clustering_from_tsv(clustering_file);
    if (!clustering.size()) {
        std::cerr << "Error opening clustering file: " << clustering_file << std::endl;
        return 1;
    }

    // Print statistics
    if (verbose) {
        std::cout << "Number of nodes: " << graph.node_count() << std::endl;
        std::cout << "Number of edges: " << graph.edge_count() << std::endl;
        std::cout << "Number of clusters: " << clustering.get_num_clusters() << std::endl;

        std::cout << std::endl;
        std::cout << "Removing singletons..." << std::endl;
    }

    // Fetch singletons
    std::vector<int> singletons = clustering.get_singletons();

    // Extract subgraph
    Graph<int> clustered_subgraph = SubgraphExtractor::extract_clustered_subgraph(graph, clustering);

    // Print clustered subgraph statistics
    if (verbose) {
        std::cout << "Number of nodes in subgraph: " << clustered_subgraph.node_count() << std::endl;
        std::cout << "Number of edges in subgraph: " << clustered_subgraph.edge_count() << std::endl;
        std::cout << "Number of singletons: " << singletons.size() << std::endl;
        std::cout << "Number of clusters remaining: " << clustering.get_num_clusters() << std::endl;
    }

    // Start computing the Stochastic Block Model
    StochasticBlockModel sbm(clustered_subgraph, clustering);
    if (verbose) {
        std::cout << std::endl;
        std::cout << "Computing Stochastic Block Model..." << std::endl;
    }

    // Generate graph and print the graph statistics
    std::vector<int> block_assignments = clustering.get_block_assignments();
    Graph<int> sbm_graph = sbm.generate_graph(clustered_subgraph.node_count(), block_assignments);
    sbm_graph.remove_floating_nodes();
    sbm_graph.remove_self_loops();
    sbm_graph.remove_parallel_edges();
    if (verbose) {
        std::cout << "Graph generated from SBM." << std::endl;
        std::cout << "Number of nodes in SBM graph: " << sbm_graph.node_count() << std::endl;
        std::cout << "Number of edges in SBM graph: " << sbm_graph.edge_count() << std::endl;
    }

    // Create a clustering for the SBM graph
    Clustering sbm_clustering;
    for (auto node_it = sbm_graph.begin(); node_it != sbm_graph.end(); ++node_it) {
        int node_id = *node_it;
        int cluster_id = block_assignments[node_id];
        sbm_clustering.add_node_to_cluster(node_id, cluster_id);
    }
    if (verbose) {
        std::cout << "Clustering created for SBM graph." << std::endl;
    }

    return 0;
}
