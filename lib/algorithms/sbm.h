#ifndef SBM_H
#define SBM_H

#include <string>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <thread>
#include <chrono>
#include <filesystem>
#include "../data_structures/graph.h"
#include "../io/graph_io.h"

namespace fs = std::filesystem;

class SBMGenerator {
public:
    // Generate an SBM from a graph and clustering, and return the resulting graph
    static CSRGraph generate_sbm(
        const CSRGraph& input_graph,
        const std::string& clustering_filepath,
        const std::string& output_directory,
        const std::string& prefix,
        int num_threads = 1,
        bool verbose = false) {
        
        // Create output directory if it doesn't exist
        if (!fs::exists(output_directory)) {
            fs::create_directories(output_directory);
        }
        
        // Generate temporary edge list filename
        std::string edge_list_filename = output_directory + "/" + prefix + "_edges.tsv";
        
        // Save the input graph as an edge list
        save_graph_edgelist(edge_list_filename, input_graph, verbose);
        
        // Build the command to call the Python script
        std::string command = "python extlib/gen_SBM.py"
                             + std::string(" -f ") + edge_list_filename
                             + std::string(" -c ") + clustering_filepath
                             + std::string(" -o ") + output_directory;
                             
        if (verbose) {
            std::cout << "Running command: " << command << std::endl;
        }
        
        // Execute the command
        auto start_time = std::chrono::steady_clock::now();
        int ret = system(command.c_str());
        auto end_time = std::chrono::steady_clock::now();
        
        if (ret != 0) {
            std::cerr << "Error running SBM generator script" << std::endl;
            return CSRGraph(); // Return empty graph on error
        }
        
        if (verbose) {
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
            std::cout << "SBM generation completed in " << duration.count() / 1000.0 
                      << " seconds" << std::endl;
        }
        
        // Load the generated graph
        std::string output_graph_filename = output_directory + "/syn_sbm.tsv";
        
        if (verbose) {
            std::cout << "Loading generated SBM from " << output_graph_filename << std::endl;
        }
        
        // Load the generated graph
        CSRGraph sbm_graph = load_undirected_tsv_edgelist_parallel(
            output_graph_filename, num_threads, verbose);
        
        clean_graph_parallel(sbm_graph, num_threads, verbose);
        
        if (verbose) {
            std::cout << "Generated SBM has " << sbm_graph.num_nodes << " nodes and " 
                      << sbm_graph.num_edges << " edges" << std::endl;
        }
        
        return sbm_graph;
    }
    
    // Fork a process to generate an SBM and return immediately
    static pid_t fork_generate_sbm(
        const CSRGraph& input_graph,
        const std::string& clustering_filepath,
        const std::string& output_directory,
        const std::string& prefix,
        bool verbose = false) {
        
        pid_t pid = fork();
        
        if (pid < 0) {
            // Error
            std::cerr << "Error forking process for SBM generation" << std::endl;
            return -1;
        } else if (pid == 0) {
            // Child process
            
            // Create output directory if it doesn't exist
            if (!fs::exists(output_directory)) {
                fs::create_directories(output_directory);
            }
            
            // Generate temporary edge list filename
            std::string edge_list_filename = output_directory + "/" + prefix + "_edges.tsv";
            
            // Save the input graph as an edge list
            save_graph_edgelist(edge_list_filename, input_graph, verbose);
            
            // Build the command to call the Python script
            std::string command = "python extlib/gen_SBM.py"
                                + std::string(" -f ") + edge_list_filename
                                + std::string(" -c ") + clustering_filepath
                                + std::string(" -o ") + output_directory;
                                
            if (verbose) {
                std::cout << "[Child " << getpid() << "] Running command: " << command << std::endl;
            }
            
            // Execute the command
            int ret = system(command.c_str());
            
            if (ret != 0) {
                std::cerr << "[Child " << getpid() << "] Error running SBM generator script" << std::endl;
                exit(1);
            }
            
            if (verbose) {
                std::cout << "[Child " << getpid() << "] SBM generation completed" << std::endl;
            }
            
            // Exit child process
            exit(0);
        } else {
            // Parent process
            if (verbose) {
                std::cout << "Forked SBM generation process with PID " << pid << std::endl;
            }
            return pid;
        }
    }
    
    // Wait for a forked SBM generation process to complete and load the resulting graph
    static CSRGraph wait_and_load_sbm(
        pid_t pid,
        const std::string& output_directory,
        int num_threads = 1,
        bool verbose = false) {
        
        int status;
        waitpid(pid, &status, 0);
        
        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
            // Child process completed successfully
            if (verbose) {
                std::cout << "SBM generation process " << pid << " completed successfully" << std::endl;
            }
            
            // Load the generated graph
            std::string output_graph_filename = output_directory + "/syn_sbm.tsv";
            
            if (verbose) {
                std::cout << "Loading generated SBM from " << output_graph_filename << std::endl;
            }
            
            // Load the generated graph
            CSRGraph sbm_graph = load_undirected_tsv_edgelist_parallel(
                output_graph_filename, num_threads, verbose);
            
            clean_graph_parallel(sbm_graph, num_threads, verbose);
            
            if (verbose) {
                std::cout << "Generated SBM has " << sbm_graph.num_nodes << " nodes and " 
                          << sbm_graph.num_edges << " edges" << std::endl;
            }
            
            return sbm_graph;
        } else {
            // Child process failed
            std::cerr << "SBM generation process " << pid << " failed" << std::endl;
            return CSRGraph(); // Return empty graph on error
        }
    }
};

#endif // SBM_H
