#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "io/graph_io.h"
#include "utils/subgraph_extractor.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <edgelist_tsv> <clustering_tsv>" << std::endl;
        return 1;
    }

    std::string edgelist_file = argv[1];
    std::string clustering_file = argv[2];

    std::cout << "Reading graph and clustering..." << std::endl;

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
    std::cout << "Number of nodes: " << graph.node_count() << std::endl;
    std::cout << "Number of edges: " << graph.edge_count() << std::endl;
    std::cout << "Number of clusters: " << clustering.get_num_clusters() << std::endl;

    std::cout << std::endl;
    std::cout << "Removing singletons..." << std::endl;

    // Fetch singletons
    std::vector<int> singletons = clustering.get_singletons();

    // Extract subgraph
    Graph<int> clustered_subgraph = SubgraphExtractor::extract_clustered_subgraph(graph, clustering);
    std::cout << "Number of nodes in subgraph: " << clustered_subgraph.node_count() << std::endl;
    std::cout << "Number of edges in subgraph: " << clustered_subgraph.edge_count() << std::endl;
    std::cout << "Number of singletons: " << singletons.size() << std::endl;
    std::cout << "Number of clusters remaining: " << clustering.get_num_clusters() << std::endl;

    return 0;
}
