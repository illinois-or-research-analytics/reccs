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
#include "algorithm/sbm_gt.h"

#include <networkit/graph/Graph.hpp>
#include <networkit/io/EdgeListReader.hpp>


#include <iostream>
#include <igraph.h>
#include <string>

int main(int argc, char* argv[]) {
    // Initialize the igraph library
    igraph_set_attribute_table(&igraph_cattribute_table);

    // Read command line arguments
    bool verbose = false;
    std::string edgelist_file;
    std::string clustering_file;
    std::string method = "gt"; // Default method
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

        // Optional argument for sbm generation method
        if (std::string(argv[i]) == "-sbm" || std::string(argv[i]) == "--sbm") {
            if (i + 1 < argc) {
                // Read the method argument
                method = argv[++i];
                if (method != "gt" && method != "igraph") {
                    std::cerr << "Error: Invalid method. Use 'gt' or 'igraph'." << std::endl;
                    return 1;
                }
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
    std::unordered_map<int, int> id_to_index; // Map to store original IDs to contiguous indices
    auto graph = io.load_graph_from_tsv(edgelist_file, id_to_index);
    
    // Read clustering
    auto clustering = io.load_clustering_from_tsv(clustering_file);
    if (!clustering.size()) {
        std::cerr << "Error opening clustering file: " << clustering_file << std::endl;
        return 1;
    }

    // Print graph and clustering information
    if (verbose) {
        std::cout << "Graph loaded with " << graph.get_num_nodes() << " vertices and "
                  << graph.get_num_edges() << " edges." << std::endl;
        std::cout << "Clustering loaded with " << clustering.get_num_clusters() << " clusters." << std::endl;
    }

    return 0;
}
