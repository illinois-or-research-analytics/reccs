#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <igraph/igraph.h>

#include "io/graph_io.h"
#include "utils/subgraph_extractor.h"
#include "algorithm/sbm.h"


int main(int argc, char* argv[]) {
    // Initialize the igraph library
    igraph_i_set_attribute_table(&igraph_cattribute_table);

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

    // Read graph
    auto graph = io.load_graph_from_tsv(edgelist_file);
    
    // Read clustering
    auto clustering = io.load_clustering_from_tsv(clustering_file);
    if (!clustering.size()) {
        std::cerr << "Error opening clustering file: " << clustering_file << std::endl;
        return 1;
    }

    // Print graph and clustering information
    if (verbose) {
        std::cout << "Graph loaded with " << igraph_vcount(&graph) << " vertices and "
                  << igraph_ecount(&graph) << " edges." << std::endl;
        std::cout << "Clustering loaded with " << clustering.size() << " clusters." << std::endl;
    }

    return 0;
}
