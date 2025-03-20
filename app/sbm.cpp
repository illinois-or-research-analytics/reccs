// graph_tool_sbm_integration.cpp
// Integration of graph-tool's C++ backend for SBM network generation

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <vector>
#include <iostream>
#include <random>
#include <cmath>

// You'll need to include these graph-tool headers
// Adjust paths according to your installation
#include <graph_tool/graph.hh>
#include <graph_tool/generation/generate_sbm.hh>

using namespace std;
using namespace boost;
using namespace graph_tool;

// Custom SBM network generator using graph-tool's backend
class SBMNetworkGenerator {
private:
    // Random number generator
    boost::mt19937 gen;
    
    // Graph type definition from graph-tool
    typedef graph_tool::GraphInterface::multigraph_t Graph;
    
public:
    SBMNetworkGenerator(unsigned long seed = 42) : gen(seed) {}
    
    Graph generate_sbm_network(
        const vector<vector<double>>& probability_matrix,  // Block-to-block edge probabilities
        const vector<int>& cluster_assignments,            // Node-to-block assignments
        const vector<int>& degree_sequence,                // Target degree sequence
        bool directed = false,
        bool self_loops = false) {
        
        size_t num_vertices = cluster_assignments.size();
        
        // Verify inputs
        if (degree_sequence.size() != num_vertices) {
            throw std::invalid_argument("Degree sequence length must match the number of vertices");
        }
        
        // Initialize the graph
        Graph g(directed);
        
        // Add vertices to the graph
        for (size_t i = 0; i < num_vertices; ++i) {
            add_vertex(g);
        }
        
        // Convert cluster assignments to a property map
        auto block_property = get_vertex_property<int32_t>(g, "block");
        for (size_t i = 0; i < num_vertices; ++i) {
            block_property[vertex(i, g)] = cluster_assignments[i];
        }
        
        // Convert probability matrix to a format accepted by graph-tool
        size_t num_blocks = probability_matrix.size();
        array2d<double> probs(extents[num_blocks][num_blocks]);
        
        for (size_t r = 0; r < num_blocks; ++r) {
            if (probability_matrix[r].size() != num_blocks) {
                throw std::invalid_argument("Probability matrix must be square");
            }
            
            for (size_t c = 0; c < num_blocks; ++c) {
                probs[r][c] = probability_matrix[r][c];
            }
        }
        
        // Create a property map for target degrees
        auto degree_property = get_vertex_property<int32_t>(g, "target_degree");
        for (size_t i = 0; i < num_vertices; ++i) {
            degree_property[vertex(i, g)] = degree_sequence[i];
        }
        
        // Generate the SBM network with graph-tool's algorithm
        // We use a variant of the algorithm that also accounts for degree sequence
        generate_sbm(g, block_property, probs, degree_property, self_loops, gen);
        
        return g;
    }
    
    // Utility function to export the generated graph to edge list
    void export_edges(const Graph& g, const string& filename) {
        ofstream out(filename);
        
        graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e) {
            out << source(*e, g) << " " << target(*e, g) << endl;
        }
        
        out.close();
    }
    
    // Utility function to calculate graph statistics
    void print_graph_stats(const Graph& g) {
        size_t num_vertices = num_vertices(g);
        size_t num_edges = num_edges(g);
        
        cout << "Graph Statistics:" << endl;
        cout << "  Vertices: " << num_vertices << endl;
        cout << "  Edges: " << num_edges << endl;
        
        // Calculate degree statistics
        vector<size_t> degrees(num_vertices);
        graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v) {
            degrees[*v] = out_degree(*v, g);
        }
        
        // Calculate average degree
        double avg_degree = 0;
        for (auto d : degrees) {
            avg_degree += d;
        }
        avg_degree /= num_vertices;
        
        cout << "  Average degree: " << avg_degree << endl;
    }
};

// Example usage
int main() {
    // Example inputs
    vector<vector<double>> prob_matrix = {
        {0.8, 0.2},
        {0.2, 0.6}
    };
    
    // Assign 10 nodes to block 0 and 10 nodes to block 1
    vector<int> cluster_assign(20);
    fill(cluster_assign.begin(), cluster_assign.begin() + 10, 0);
    fill(cluster_assign.begin() + 10, cluster_assign.end(), 1);
    
    // Target degree sequence (for simplicity, uniform)
    vector<int> degrees(20, 5);
    
    // Create the generator and generate the network
    SBMNetworkGenerator generator;
    auto graph = generator.generate_sbm_network(prob_matrix, cluster_assign, degrees);
    
    // Print stats and export
    generator.print_graph_stats(graph);
    generator.export_edges(graph, "sbm_network.edgelist");
    
    return 0;
}
