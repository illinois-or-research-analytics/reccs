// lib/utils/orchestrator.h
#pragma once

#include <string>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <sys/wait.h>
#include <unistd.h>

/**
 * @class Orchestrator
 * @brief Handles orchestration of the RECCS workflow, including splitting graphs and running 
 *        various analysis scripts in parallel.
 */
class Orchestrator {
public:
    /**
     * @brief Constructor for the Orchestrator
     * 
     * @param graph_filename Path to the input graph file
     * @param cluster_filename Path to the input clustering file
     * @param temp_dir Directory for temporary files
     * @param verbose Whether to print verbose output
     */
    Orchestrator(const std::string& graph_filename, const std::string& cluster_filename, 
                const std::string& temp_dir, bool verbose)
        : graph_filename_(graph_filename), cluster_filename_(cluster_filename), 
          temp_dir_(temp_dir), verbose_(verbose) {
        // Ensure temp directory exists
        if (!std::filesystem::exists(temp_dir_)) {
            std::filesystem::create_directories(temp_dir_);
        }
    }

    /**
     * @brief Run the full RECCS workflow
     * 
     * @return int 0 if successful, non-zero error code otherwise
     */
    int run() {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Step 1: Run the splitter script
        if (!splitGraph()) {
            return 1;
        }
        
        if (verbose_) {
            auto split_end_time = std::chrono::high_resolution_clock::now();
            auto split_duration = std::chrono::duration_cast<std::chrono::seconds>(
                split_end_time - start_time).count();
            std::cout << "Graph splitting completed in " << split_duration << " seconds" << std::endl;
            std::cout << "Step 2: Processing clustered and unclustered components in parallel..." << std::endl;
        }
        
        // Step 2: Run stats and sbm in parallel
        if (!runParallelProcessing()) {
            return 1;
        }
        
        if (verbose_) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(
                end_time - start_time).count();
            std::cout << "Total execution time: " << duration << " seconds" << std::endl;
            
            // List the generated files
            std::cout << "Generated files:" << std::endl;
            std::cout << "  - " << temp_dir_ << "/clustered_stats.csv" << std::endl;
            std::cout << "  - " << temp_dir_ << "/clustered_sbm/syn_sbm.tsv" << std::endl;
            std::cout << "  - " << temp_dir_ << "/unclustered_sbm/syn_sbm.tsv" << std::endl;
        }
        
        return 0;
    }

private:
    std::string graph_filename_;
    std::string cluster_filename_;
    std::string temp_dir_;
    bool verbose_;
    
    // Paths for intermediate files
    std::string clustered_edges_;
    std::string clustered_clusters_;
    std::string unclustered_edges_;
    std::string unclustered_clusters_;

    /**
     * @brief Run the graph splitter script
     * 
     * @return bool true if successful, false otherwise
     */
    bool splitGraph() {
        if (verbose_) {
            std::cout << "Step 1: Running graph splitter..." << std::endl;
        }
        
        std::string splitter_command = "python3 extlib/splitter.py";
        splitter_command += " -f " + graph_filename_;
        splitter_command += " -c " + cluster_filename_;
        splitter_command += " -o " + temp_dir_;
        if (verbose_) {
            splitter_command += " -v";
            std::cout << "Executing: " << splitter_command << std::endl;
        }
        
        int splitter_status = system(splitter_command.c_str());
        if (splitter_status != 0) {
            std::cerr << "Error: Splitter script failed with status " << splitter_status << std::endl;
            return false;
        }
        
        // Define paths to output files
        clustered_edges_ = temp_dir_ + "/non_singleton_edges.tsv";
        clustered_clusters_ = temp_dir_ + "/non_singleton_clusters.tsv";
        unclustered_edges_ = temp_dir_ + "/singleton_edges.tsv";
        unclustered_clusters_ = temp_dir_ + "/singleton_clusters.tsv";
        
        // Check if output files exist
        if (!std::filesystem::exists(clustered_edges_) || 
            !std::filesystem::exists(clustered_clusters_) || 
            !std::filesystem::exists(unclustered_edges_) || 
            !std::filesystem::exists(unclustered_clusters_)) {
            std::cerr << "Error: Expected output files from splitter not found" << std::endl;
            return false;
        }
        
        return true;
    }
    
    /**
     * @brief Run the stats and SBM scripts in parallel
     * 
     * @return bool true if successful, false otherwise
     */
    bool runParallelProcessing() {
        pid_t pid1, pid2, pid3;
        
        // Fork for stats on clustered subgraph
        pid1 = fork();
        if (pid1 == 0) {
            // Child process 1: Run stats on the clustered subgraph
            std::string stats_command = "python3 extlib/stats.py";
            stats_command += " -i " + clustered_edges_;
            stats_command += " -e " + clustered_clusters_;
            stats_command += " -o " + temp_dir_ + "/clustered_stats.csv";
            if (verbose_) {
                stats_command += " -v";
                std::cout << "Child process 1 executing: " << stats_command << std::endl;
            }
            
            int stats_status = system(stats_command.c_str());
            exit(stats_status);
        } else if (pid1 < 0) {
            std::cerr << "Error: Failed to fork for stats process" << std::endl;
            return false;
        }
        
        // Fork for SBM on clustered subgraph
        pid2 = fork();
        if (pid2 == 0) {
            // Child process 2: Run SBM on the clustered subgraph
            std::string sbm_clustered_command = "python3 extlib/gen_SBM.py";
            sbm_clustered_command += " -f " + clustered_edges_;
            sbm_clustered_command += " -c " + clustered_clusters_;
            sbm_clustered_command += " -o " + temp_dir_ + "/clustered_sbm";
            if (verbose_) {
                std::cout << "Child process 2 executing: " << sbm_clustered_command << std::endl;
            }
            
            int sbm_clustered_status = system(sbm_clustered_command.c_str());
            exit(sbm_clustered_status);
        } else if (pid2 < 0) {
            std::cerr << "Error: Failed to fork for clustered SBM process" << std::endl;
            return false;
        }
        
        // Fork for SBM on unclustered subgraph
        pid3 = fork();
        if (pid3 == 0) {
            // Child process 3: Run SBM on the unclustered subgraph
            std::string sbm_unclustered_command = "python3 extlib/gen_SBM.py";
            sbm_unclustered_command += " -f " + unclustered_edges_;
            sbm_unclustered_command += " -c " + unclustered_clusters_;
            sbm_unclustered_command += " -o " + temp_dir_ + "/unclustered_sbm";
            if (verbose_) {
                std::cout << "Child process 3 executing: " << sbm_unclustered_command << std::endl;
            }
            
            int sbm_unclustered_status = system(sbm_unclustered_command.c_str());
            exit(sbm_unclustered_status);
        } else if (pid3 < 0) {
            std::cerr << "Error: Failed to fork for unclustered SBM process" << std::endl;
            return false;
        }
        
        // Wait for all child processes to complete
        int status1, status2, status3;
        waitpid(pid1, &status1, 0);
        waitpid(pid2, &status2, 0);
        waitpid(pid3, &status3, 0);
        
        // Check process exit statuses
        if (WEXITSTATUS(status1) != 0) {
            std::cerr << "Error: Stats process failed with status " << WEXITSTATUS(status1) << std::endl;
            return false;
        }
        if (WEXITSTATUS(status2) != 0) {
            std::cerr << "Error: Clustered SBM process failed with status " << WEXITSTATUS(status2) << std::endl;
            return false;
        }
        if (WEXITSTATUS(status3) != 0) {
            std::cerr << "Error: Unclustered SBM process failed with status " << WEXITSTATUS(status3) << std::endl;
            return false;
        }
        
        return true;
    }
};
