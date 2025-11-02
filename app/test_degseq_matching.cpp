/**
 * @file test_degseq_matching.cpp
 * @brief Test script for degree sequence matching algorithm
 *
 * This script tests the deg_seq_matching_pp algorithm by:
 * 1. Loading a reference graph and an SBM graph
 * 2. Computing degree deficits (reference - SBM)
 * 3. Running degree sequence matching to add edges to the SBM graph
 * 4. Outputting the augmented graph
 *
 * Usage:
 *   ./test_degseq_matching <reference.tsv> <sbm.tsv> <clustering.tsv> <output.tsv> [options]
 *
 * Options:
 *   -v, --verbose    Verbose output
 *   -t <threads>     Number of threads (default: 1 for testing)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <chrono>
#include <omp.h>

#include "../lib/data_structures/graph.h"
#include "../lib/data_structures/available_node_degrees.h"
#include "../lib/io/g_io.h"
#include "../lib/io/cluster_io.h"
#include "../lib/algorithm/deg_seq_matching_pp.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

/**
 * @brief Hash function for pairs (for unordered_set of edges)
 */
template <typename T1, typename T2>
struct pair_hash {
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

/**
 * @brief Compute degree deficits between reference and SBM graphs
 */
std::unordered_map<uint64_t, int32_t> compute_degree_deficits(
    const Graph& reference_graph,
    const Graph& sbm_graph,
    bool verbose)
{
    if (verbose) {
        std::cout << "\n=== Computing Degree Deficits ===\n";
    }

    std::unordered_map<uint64_t, int32_t> deficits;

    // Compute reference degrees
    std::unordered_map<uint64_t, uint32_t> reference_degrees;
    for (uint32_t node = 0; node < reference_graph.num_nodes; ++node) {
        uint64_t node_id = reference_graph.id_map[node];
        uint32_t degree = reference_graph.row_ptr[node + 1] - reference_graph.row_ptr[node];
        reference_degrees[node_id] = degree;
    }

    // Compute SBM degrees
    std::unordered_map<uint64_t, uint32_t> sbm_degrees;
    for (uint32_t node = 0; node < sbm_graph.num_nodes; ++node) {
        uint64_t node_id = sbm_graph.id_map[node];
        uint32_t degree = sbm_graph.row_ptr[node + 1] - sbm_graph.row_ptr[node];
        sbm_degrees[node_id] = degree;
    }

    // Compute deficits
    std::unordered_set<uint64_t> all_nodes;
    for (const auto& [node_id, _] : reference_degrees) {
        all_nodes.insert(node_id);
    }
    for (const auto& [node_id, _] : sbm_degrees) {
        all_nodes.insert(node_id);
    }

    size_t positive_deficits = 0;
    size_t negative_deficits = 0;
    size_t zero_deficits = 0;
    int64_t total_deficit = 0;

    for (uint64_t node_id : all_nodes) {
        int32_t ref_deg = reference_degrees.count(node_id) ? reference_degrees[node_id] : 0;
        int32_t sbm_deg = sbm_degrees.count(node_id) ? sbm_degrees[node_id] : 0;
        int32_t deficit = ref_deg - sbm_deg;

        if (deficit != 0) {
            deficits[node_id] = deficit;
            total_deficit += deficit;

            if (deficit > 0) {
                positive_deficits++;
            } else {
                negative_deficits++;
            }
        } else {
            zero_deficits++;
        }
    }

    if (verbose) {
        std::cout << "Total nodes: " << all_nodes.size() << "\n";
        std::cout << "Nodes with positive deficit: " << positive_deficits << "\n";
        std::cout << "Nodes with negative deficit: " << negative_deficits << "\n";
        std::cout << "Nodes with zero deficit: " << zero_deficits << "\n";
        std::cout << "Total deficit sum: " << total_deficit << "\n";
    }

    return deficits;
}

/**
 * @brief Save deficits to JSON file
 */
void save_deficits_json(
    const std::unordered_map<uint64_t, int32_t>& deficits,
    const std::string& filename)
{
    json j;
    for (const auto& [node_id, deficit] : deficits) {
        j[std::to_string(node_id)] = deficit;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
    file << j.dump(2);
    file.close();
}

/**
 * @brief Print usage information
 */
void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name
              << " <reference.tsv> <sbm.tsv> <clustering.tsv> <output.tsv> [options]\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  reference.tsv   Reference graph edgelist (TSV format)\n";
    std::cout << "  sbm.tsv         SBM graph edgelist to augment (TSV format)\n";
    std::cout << "  clustering.tsv  Node clustering file (TSV format: node_id<TAB>cluster_id)\n";
    std::cout << "  output.tsv      Output augmented graph (TSV format)\n\n";
    std::cout << "Options:\n";
    std::cout << "  -v, --verbose   Enable verbose output\n";
    std::cout << "  -t <threads>    Number of threads (default: 1)\n";
    std::cout << "  -h, --help      Show this help message\n";
    std::cout << "  --save-deficits <file>  Save computed deficits to JSON file\n";
}

int main(int argc, char* argv[]) {
    // Parse arguments
    if (argc < 5) {
        print_usage(argv[0]);
        return 1;
    }

    std::string reference_file = argv[1];
    std::string sbm_file = argv[2];
    std::string clustering_file = argv[3];
    std::string output_file = argv[4];

    bool verbose = false;
    int num_threads = 1;
    std::string deficits_output = "";

    // Parse options
    for (int i = 5; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-t" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "--save-deficits" && i + 1 < argc) {
            deficits_output = argv[++i];
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    // Set thread count
    omp_set_num_threads(num_threads);

    auto start_time = std::chrono::high_resolution_clock::now();

    if (verbose) {
        std::cout << "\n=== Degree Sequence Matching Test ===\n";
        std::cout << "Reference graph: " << reference_file << "\n";
        std::cout << "SBM graph: " << sbm_file << "\n";
        std::cout << "Clustering: " << clustering_file << "\n";
        std::cout << "Output: " << output_file << "\n";
        std::cout << "Threads: " << num_threads << "\n";
    }

    // Load reference graph
    if (verbose) {
        std::cout << "\n=== Loading Reference Graph ===\n";
    }
    Graph reference_graph = load_undirected_tsv_edgelist_parallel(reference_file, num_threads, verbose);

    if (verbose) {
        std::cout << "Reference graph loaded: " << reference_graph.num_nodes << " nodes, "
                  << reference_graph.num_edges / 2 << " edges\n";
    }

    // Load SBM graph
    if (verbose) {
        std::cout << "\n=== Loading SBM Graph ===\n";
    }
    Graph sbm_graph = load_undirected_tsv_edgelist_parallel(sbm_file, num_threads, verbose);

    if (verbose) {
        std::cout << "SBM graph loaded: " << sbm_graph.num_nodes << " nodes, "
                  << sbm_graph.num_edges / 2 << " edges\n";
    }

    // Load clustering
    if (verbose) {
        std::cout << "\n=== Loading Clustering ===\n";
    }
    Clustering clustering = load_clustering(clustering_file, sbm_graph, verbose);

    if (verbose) {
        std::cout << "Clustering loaded: " << clustering.cluster_nodes.size() << " clusters\n";
    }

    // Compute degree deficits
    auto deficits = compute_degree_deficits(reference_graph, sbm_graph, verbose);

    // Save deficits if requested
    if (!deficits_output.empty()) {
        if (verbose) {
            std::cout << "\nSaving deficits to: " << deficits_output << "\n";
        }
        save_deficits_json(deficits, deficits_output);
    }

    // Create a temporary JSON file for the degree manager
    std::string temp_deficits_file = "/tmp/test_degseq_deficits.json";
    save_deficits_json(deficits, temp_deficits_file);

    // Initialize degree manager
    if (verbose) {
        std::cout << "\n=== Initializing Degree Manager ===\n";
    }
    auto degree_manager = std::make_shared<AvailableNodeDegreesManager>(temp_deficits_file);

    // Print initial statistics
    if (verbose) {
        auto stats = degree_manager->get_stats();
        std::cout << "\n=== Initial Degree Budget Statistics ===\n";
        std::cout << "Available nodes: " << stats.total_available_nodes << "\n";
        std::cout << "Total degree budget: " << stats.total_available_degrees << "\n";
        std::cout << "Average budget per node: " << stats.avg_available_degree << "\n";
    }

    // Process the graph as a single cluster for simplicity
    // (In production, would process each cluster separately)
    if (verbose) {
        std::cout << "\n=== Running Degree Sequence Matching ===\n";
    }

    // Collect all nodes as a single "cluster"
    std::unordered_set<uint64_t> all_node_ids;
    for (uint32_t node = 0; node < sbm_graph.num_nodes; ++node) {
        all_node_ids.insert(sbm_graph.id_map[node]);
    }

    // Create a shared pointer to the SBM graph
    auto graph_ptr = std::make_shared<Graph>(sbm_graph);

    // Create a single GraphTask for all nodes
    GraphTask task(
        graph_ptr,
        TaskType::MIN_DEG_ENFORCE,
        "all_nodes",
        0,
        0,
        degree_manager,
        all_node_ids
    );

    // Run degree sequence matching
    uint64_t edges_before = graph_ptr->num_edges;
    match_degree_sequence_pp(task);
    uint64_t edges_after = graph_ptr->num_edges;

    uint64_t total_edges_added = (edges_after - edges_before) / 2;  // Divide by 2 for undirected

    // Update the original graph with the new edges
    sbm_graph = *graph_ptr;

    // Print final statistics
    if (verbose) {
        auto final_stats = degree_manager->get_stats();
        std::cout << "\n=== Final Degree Budget Statistics ===\n";
        std::cout << "Remaining available nodes: " << final_stats.total_available_nodes << "\n";
        std::cout << "Remaining degree budget: " << final_stats.total_available_degrees << "\n";
        std::cout << "Budget utilization: "
                  << (1.0 - (double)final_stats.total_available_degrees /
                      degree_manager->get_stats().total_available_degrees) * 100.0 << "%\n";
        std::cout << "\nTotal edges added: " << total_edges_added << "\n";
    }

    // Write output graph
    if (verbose) {
        std::cout << "\n=== Writing Output ===\n";
    }

    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Failed to open output file: " << output_file << "\n";
        return 1;
    }

    // Write all edges from the augmented graph
    std::unordered_set<std::pair<uint64_t, uint64_t>, pair_hash<uint64_t, uint64_t>> written_edges;

    for (uint32_t src = 0; src < sbm_graph.num_nodes; ++src) {
        uint64_t src_id = sbm_graph.id_map[src];
        for (uint32_t idx = sbm_graph.row_ptr[src]; idx < sbm_graph.row_ptr[src + 1]; ++idx) {
            uint32_t tgt = sbm_graph.col_idx[idx];
            uint64_t tgt_id = sbm_graph.id_map[tgt];

            // Only write each edge once (for undirected graph)
            auto edge_pair = std::make_pair(std::min(src_id, tgt_id), std::max(src_id, tgt_id));
            if (written_edges.find(edge_pair) == written_edges.end()) {
                out << src_id << "\t" << tgt_id << "\n";
                written_edges.insert(edge_pair);
            }
        }
    }

    out.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    if (verbose) {
        std::cout << "\n=== Processing Complete ===\n";
        std::cout << "Total execution time: " << duration.count() / 1000.0 << " seconds\n";
        std::cout << "Output written to: " << output_file << "\n";
    }

    // Clean up temp file
    std::remove(temp_deficits_file.c_str());

    return 0;
}
