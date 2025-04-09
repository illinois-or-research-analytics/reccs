#include <filesystem>
#include <iostream>

/**
 * @brief Stochastic Block Model that uses Peixoto's graph-tool (gt) library.
*/
class StochasticBlockModelGT {
private:
    Graph graph; // The input graph
    Clustering clustering; // The clustering of the nodes
public:
    /**
     * @brief Constructor for StochasticBlockModelGT.
     * 
     * @param graph The input graph.
     * @param clustering The clustering of the nodes.
     */
    StochasticBlockModelGT(const Graph& graph, const Clustering& clustering)
        : graph(graph), clustering(clustering) {}


    Graph generate_graph() {
        // Make a temp directory
        std::string temp_dir = "temp";
        std::filesystem::create_directory(temp_dir);

        graph_io io;

        // Write the graph to a file
        std::string graph_file = temp_dir + "/graph.tsv";
        io.write_tsv(graph, graph_file);

        // Write the clustering to a file
        std::string clustering_file = temp_dir + "/clustering.tsv";
        io.write_tsv(clustering, clustering_file);

        // Get the absolute path to the current executable
        char resolved_path[PATH_MAX];
        std::string executable_path = realpath("/proc/self/exe", resolved_path);
        if (executable_path.empty()) {
            std::cerr << "Error: could not determine executable path" << std::endl;
            return Graph();
        }

        // Get the directory containing the executable
        std::filesystem::path exec_dir = std::filesystem::path(executable_path).parent_path();

        // Construct the path to the script (../../extlib/gen_SBM.py relative to the executable)
        std::filesystem::path script_path = exec_dir / ".." / ".." / "extlib" / "gen_SBM.py";

        // Ensure the script exists
        if (!std::filesystem::exists(script_path)) {
            std::cerr << "Error: script not found at " << script_path << std::endl;
            return Graph();
        }

        // Construct the command to run the script
        std::string command = 
            "python3 " + 
            script_path.string() + 
            " -f " + 
            graph_file + 
            " -c " + 
            clustering_file;
        std::cout << "Running command: " << command << std::endl;

        // Execute the command
        int result = system(command.c_str());
        if (result != 0) {
            std::cerr << "Error: script execution failed with code " << result << std::endl;
            return Graph();
        }

        // Read the generated graph
        std::string output_file = temp_dir + "/output_graph.tsv";
        Graph generated_graph = io.read_tsv(output_file);

        // Clean up temporary files
        std::filesystem::remove_all(temp_dir);

        return generated_graph;
    }
};