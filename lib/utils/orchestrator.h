// lib/utils/orchestrator.h
#pragma once

#include <string>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <fcntl.h>


/**
 * @class Orchestrator
 * @brief Handles orchestration of the RECCS workflow, optimized to run processes in parallel when possible.
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
        
        if (verbose_) {
            std::cout << "Starting RECCS workflow with optimized parallelization" << std::endl;
            std::cout << "Input graph: " << graph_filename_ << std::endl;
            std::cout << "Input clustering: " << cluster_filename_ << std::endl;
            std::cout << "Temporary directory: " << temp_dir_ << std::endl;
        }
        
        // Step 1: Run the splitter and stats in parallel
        if (!runFirstStage()) {
            return 1;
        }
        
        if (verbose_) {
            auto stage1_end_time = std::chrono::high_resolution_clock::now();
            auto stage1_duration = std::chrono::duration_cast<std::chrono::seconds>(
                stage1_end_time - start_time).count();
            std::cout << "First stage completed in " << stage1_duration << " seconds" << std::endl;
            std::cout << "Running SBM generation for both components..." << std::endl;
        }
        
        // Step 2: Run SBM generation on both components
        if (!runSecondStage()) {
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
    
    // Track process IDs
    std::vector<pid_t> pids_;

    /**
     * @brief Execute a command in a child process with proper stream handling
     * 
     * @param command The command to execute
     * @param description Description for verbose output
     * @return pid_t Process ID of the child, or -1 on error
     */
    pid_t executeCommand(const std::string& command, const std::string& description) {
        pid_t pid = fork();
        
        if (pid == 0) {
            // Child process
            
            // Close stdin to prevent the child from blocking on input
            close(STDIN_FILENO);
            
            // Redirect stdout/stderr to /dev/null if not verbose
            if (!verbose_) {
                int dev_null = open("/dev/null", O_WRONLY);
                dup2(dev_null, STDOUT_FILENO);
                dup2(dev_null, STDERR_FILENO);
                close(dev_null);
            }
            
            if (verbose_) {
                std::cout << "[" << description << "] Executing: " << command << std::endl;
            }
            
            // Execute the command
            int status = system(command.c_str());
            
            // Use _exit instead of exit to avoid flushing stdio buffers
            _exit(WEXITSTATUS(status));
        }
        
        return pid;
    }

    /**
     * @brief Run the first stage: splitter and stats in parallel
     * 
     * @return bool true if successful, false otherwise
     */
    bool runFirstStage() {
        // Prepare commands
        std::string splitter_command = "python3 extlib/splitter.py";
        splitter_command += " -f " + graph_filename_;
        splitter_command += " -c " + cluster_filename_;
        splitter_command += " -o " + temp_dir_;
        if (verbose_) {
            splitter_command += " -v";
        }
        
        std::string stats_command = "python3 extlib/stats.py";
        stats_command += " -i " + graph_filename_;
        stats_command += " -e " + cluster_filename_;
        stats_command += " -o " + temp_dir_ + "/clustered_stats.csv";
        if (verbose_) {
            stats_command += " -v";
        }
        
        // Execute commands in parallel
        pid_t splitter_pid = executeCommand(splitter_command, "Splitter");
        if (splitter_pid < 0) {
            std::cerr << "Error: Failed to fork for splitter process" << std::endl;
            return false;
        }
        
        pid_t stats_pid = executeCommand(stats_command, "Stats");
        if (stats_pid < 0) {
            std::cerr << "Error: Failed to fork for stats process" << std::endl;
            return false;
        }
        
        // Store PIDs
        pids_.push_back(stats_pid);
        
        // Wait for the splitter process to complete (we need its output for the next stage)
        int splitter_status;
        waitpid(splitter_pid, &splitter_status, 0);
        
        if (WEXITSTATUS(splitter_status) != 0) {
            std::cerr << "Error: Splitter process failed with status " << WEXITSTATUS(splitter_status) << std::endl;
            return false;
        }
        
        // Define paths to output files from the splitter
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
     * @brief Run the second stage: SBM generation for both components
     * 
     * @return bool true if successful, false otherwise
     */
    bool runSecondStage() {
        // Prepare commands
        std::string sbm_clustered_command = "python3 extlib/gen_SBM.py";
        sbm_clustered_command += " -f " + clustered_edges_;
        sbm_clustered_command += " -c " + clustered_clusters_;
        sbm_clustered_command += " -o " + temp_dir_ + "/clustered_sbm";
        
        std::string sbm_unclustered_command = "python3 extlib/gen_SBM.py";
        sbm_unclustered_command += " -f " + unclustered_edges_;
        sbm_unclustered_command += " -c " + unclustered_clusters_;
        sbm_unclustered_command += " -o " + temp_dir_ + "/unclustered_sbm";
        
        // Execute commands in parallel
        pid_t clustered_sbm_pid = executeCommand(sbm_clustered_command, "SBM-Clustered");
        if (clustered_sbm_pid < 0) {
            std::cerr << "Error: Failed to fork for clustered SBM process" << std::endl;
            return false;
        }
        
        pid_t unclustered_sbm_pid = executeCommand(sbm_unclustered_command, "SBM-Unclustered");
        if (unclustered_sbm_pid < 0) {
            std::cerr << "Error: Failed to fork for unclustered SBM process" << std::endl;
            return false;
        }
        
        // Store PIDs
        pids_.push_back(clustered_sbm_pid);
        pids_.push_back(unclustered_sbm_pid);
        
        // Wait for all remaining child processes to complete
        int status;
        for (pid_t pid : pids_) {
            waitpid(pid, &status, 0);
            if (WEXITSTATUS(status) != 0) {
                std::cerr << "Error: Process " << pid << " failed with status " << WEXITSTATUS(status) << std::endl;
                return false;
            }
        }
        
        return true;
    }
};
