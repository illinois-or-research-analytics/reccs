#include <networkit/graph/Graph.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <iostream>
#include <vector>

int main() {
    // Create the reader with tab as delimiter
    NetworKit::EdgeListReader reader('\t', 0, "#", false, false);
    
    // The last parameter (true) enables mapping of non-continuous node IDs
    NetworKit::Graph G = reader.read("your_edgelist.tsv", true);
    
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
