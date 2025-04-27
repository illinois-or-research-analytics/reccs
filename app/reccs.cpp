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


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <edgelist_file>" << std::endl;
        return 1;
    }

    // Create the reader with tab as delimiter
    // Last parameter (true) enables mapping of non-continuous node IDs
    NetworKit::EdgeListReader reader('\t', 0, "#", false, true);
    
    // Read the graph
    NetworKit::Graph G = reader.read(argv[1]);
    
    // At this point, NetworKit has already created a mapping internally
    // You can access it only if needed
    // std::map<std::string, NetworKit::node> nodeMap = reader.getNodeMap();
    
    // Only create reverse mapping if you actually need it for your application
    // This step can be omitted if not needed for downstream processing
    /*
    std::vector<int> indexToOriginal(G.upperNodeIdBound());
    for (const auto& pair : nodeMap) {
        indexToOriginal[pair.second] = std::stoi(pair.first);
    }
    */
    
    std::cout << "Number of nodes: " << G.numberOfNodes() << std::endl;
    std::cout << "Number of edges: " << G.numberOfEdges() << std::endl;
    
    return 0;
}
