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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <edgelist_file>" << std::endl;
        return 1;
    }

    // Initialize igraph
    igraph_t graph;
    
    // Open the input file
    FILE* instream = fopen(argv[1], "r");
    if (!instream) {
        std::cerr << "Failed to open file: " << argv[1] << std::endl;
        return 1;
    }
    
    // No predefined vertex names
    igraph_strvector_t predefnames;
    igraph_strvector_init(&predefnames, 0);
    
    // Read the graph using the NCOL format
    // Parameters:
    // - graph: output graph
    // - instream: input file
    // - predefnames: predefined vertex names (empty in this case)
    // - names: true to use symbolic names in the file
    // - weights: IGRAPH_ADD_WEIGHTS_IF_PRESENT to add weights if present
    // - directed: false for undirected graph
    igraph_error_t err = igraph_read_graph_ncol(
        &graph, 
        instream, 
        &predefnames,
        true,  // Use symbolic names from the file
        IGRAPH_ADD_WEIGHTS_IF_PRESENT,  // Add weights if present
        false  // Undirected graph
    );
    
    // Close the file
    fclose(instream);
    
    if (err != IGRAPH_SUCCESS) {
        std::cerr << "Failed to read graph. Error code: " << err << std::endl;
        igraph_strvector_destroy(&predefnames);
        return 1;
    }
    
    // Print graph info
    std::cout << "Number of vertices: " << igraph_vcount(&graph) << std::endl;
    std::cout << "Number of edges: " << igraph_ecount(&graph) << std::endl;
    
    // Clean up resources
    igraph_destroy(&graph);
    igraph_strvector_destroy(&predefnames);
    
    return 0;
}