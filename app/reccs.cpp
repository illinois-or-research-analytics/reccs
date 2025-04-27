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


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <edgelist_file>" << std::endl;
        return 1;
    }

    // Create the reader with tab as delimiter
    NetworKit::EdgeListReader reader('\t', 0, "#", false, false);
    
    // The last parameter (true) enables mapping of non-continuous node IDs
    NetworKit::Graph G = reader.read(argv[1], true);
    
    // Get the node map (original ID → NetworKit index)
    std::map<std::string, NetworKit::node> nodeMap = reader.getNodeMap();
    
    // Create reverse mapping (NetworKit index → original ID)
    std::vector<int> indexToOriginal(G.upperNodeIdBound());
    for (const auto& pair : nodeMap) {
        indexToOriginal[pair.second] = std::stoi(pair.first);
    }
    
    std::cout << "Number of nodes: " << G.numberOfNodes() << std::endl;
    std::cout << "Number of edges: " << G.numberOfEdges() << std::endl;
    
    // Example of accessing original IDs
    G.forNodes([&](NetworKit::node u) {
        std::cout << "NetworKit node " << u << " corresponds to original ID " 
                  << indexToOriginal[u] << std::endl;
    });
    
    return 0;
}
