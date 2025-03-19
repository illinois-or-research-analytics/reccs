#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <edgelist_tsv> <clustering_tsv>" << std::endl;
        return 1;
    }

    std::string edgelist_file = argv[1];
    std::string clustering_file = argv[2];

    // Read edgelist
    std::ifstream edge_in(edgelist_file);
    if (!edge_in) {
        std::cerr << "Error opening edge list file: " << edgelist_file << std::endl;
        return 1;
    }

    std::unordered_set<int> nodes;
    size_t edge_count = 0;
    std::string line;
    
    while (std::getline(edge_in, line)) {
        std::istringstream iss(line);
        int src, dst;
        if (iss >> src >> dst) {
            nodes.insert(src);
            nodes.insert(dst);
            edge_count++;
        }
    }
    edge_in.close();

    // Read clustering
    std::ifstream cluster_in(clustering_file);
    if (!cluster_in) {
        std::cerr << "Error opening clustering file: " << clustering_file << std::endl;
        return 1;
    }

    std::unordered_set<int> clusters;
    while (std::getline(cluster_in, line)) {
        std::istringstream iss(line);
        int node_id, cluster_id;
        if (iss >> node_id >> cluster_id) {
            clusters.insert(cluster_id);
        }
    }
    cluster_in.close();

    // Print statistics
    std::cout << "Number of nodes: " << nodes.size() << std::endl;
    std::cout << "Number of edges: " << edge_count << std::endl;
    std::cout << "Number of clusters: " << clusters.size() << std::endl;

    return 0;
}
