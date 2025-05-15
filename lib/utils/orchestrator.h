#pragma once

#include <string>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <fcntl.h>
#include <algorithm>  // for std::find

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
        
        // Also ensure the SBM output directories exist
        std::filesystem::create_directories(temp_dir_ + "/clustered_sbm");
        std::filesystem::create_directories(temp_dir_ + "/unclustered_sbm");
    }

    /**
     * @brief Run the full RECCS workflow
     * 
     * @return int 0 if successful, non-zero error code otherwise
     */
    int run() {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (verbose_) {
            std::cout << "Starting RECCS workflow" << std::endl;
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
        
        // Wait for any remaining processes
        if (!waitForAllProcesses()) {
            std::cerr << "Error: Some processes failed to complete successfully" << std::endl;
            return 1;
        }
        
        // Verify all expected output files exist
        if (!verifyOutputFiles()) {
            std::cerr << "Error: Some expected output files are missing" << std::endl;
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
    
    // Track process IDs and their descriptions for better error reporting
    std::vector<std::pair<pid_t, std::string>> pids_;
    
    // Track process output pipes
    struct ProcessOutputPipes {
        pid_t pid;
        std::string description;
        int stdout_fd;
        int stderr_fd;
    };
    std::vector<ProcessOutputPipes> output_pipes_;

    /**
     * @brief Execute a command in a child process with proper stream handling and output capture
     * 
     * @param command The command to execute
     * @param description Description for verbose output
     * @return pid_t Process ID of the child, or -1 on error
     */
    pid_t executeCommand(const std::string& command, const std::string& description) {
        if (verbose_) {
            std::cout << "[" << description << "] Executing: " << command << std::endl;
        }
        
        // Create pipes for capturing stdout/stderr
        int stdout_pipe[2];
        int stderr_pipe[2];
        
        if (pipe(stdout_pipe) < 0 || pipe(stderr_pipe) < 0) {
            std::cerr << "Error: Failed to create pipes for " << description << std::endl;
            return -1;
        }
        
        pid_t pid = fork();
        
        if (pid == 0) {
            // Child process
            
            // Close read ends of pipes in child
            close(stdout_pipe[0]);
            close(stderr_pipe[0]);
            
            // Redirect stdout/stderr to pipes
            dup2(stdout_pipe[1], STDOUT_FILENO);
            dup2(stderr_pipe[1], STDERR_FILENO);
            
            // Close the write ends of pipes after duplication
            close(stdout_pipe[1]);
            close(stderr_pipe[1]);
            
            // Close stdin to prevent the child from blocking on input
            close(STDIN_FILENO);
            
            // Execute the command
            int status = system(command.c_str());
            
            // Use _exit instead of exit to avoid flushing stdio buffers
            _exit(WEXITSTATUS(status));
        } else if (pid > 0) {
            // Parent process
            
            // Close write ends of pipes in parent
            close(stdout_pipe[1]);
            close(stderr_pipe[1]);
            
            // Store pipes for later reading
            std::string process_log_dir = temp_dir_ + "/logs";
            std::filesystem::create_directories(process_log_dir);
            
            // Store PID, description and pipe file descriptors for later processing
            output_pipes_.push_back({pid, description, stdout_pipe[0], stderr_pipe[0]});
            
            if (verbose_) {
                std::cout << "[" << description << "] Started with PID: " << pid << std::endl;
            }
        } else {
            // Fork failed
            close(stdout_pipe[0]);
            close(stdout_pipe[1]);
            close(stderr_pipe[0]);
            close(stderr_pipe[1]);
        }
        
        return pid;
    }
    
    /**
     * @brief Process and save the output from a command
     * 
     * @param pid Process ID to process output for
     * @param description Process description
     * @param stdout_fd File descriptor for stdout
     * @param stderr_fd File descriptor for stderr
     */
    void processCommandOutput(pid_t pid, const std::string& description, int stdout_fd, int stderr_fd) {
        std::string process_log_dir = temp_dir_ + "/logs";
        std::string stdout_log = process_log_dir + "/" + description + "." + std::to_string(pid) + ".stdout.log";
        std::string stderr_log = process_log_dir + "/" + description + "." + std::to_string(pid) + ".stderr.log";
        
        // Read and save stdout
        char buffer[4096];
        ssize_t bytes_read;
        int stdout_file = open(stdout_log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        
        if (stdout_file >= 0) {
            while ((bytes_read = read(stdout_fd, buffer, sizeof(buffer))) > 0) {
                write(stdout_file, buffer, bytes_read);
                if (verbose_) {
                    write(STDOUT_FILENO, buffer, bytes_read);
                }
            }
            close(stdout_file);
        }
        
        // Read and save stderr
        int stderr_file = open(stderr_log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        
        if (stderr_file >= 0) {
            while ((bytes_read = read(stderr_fd, buffer, sizeof(buffer))) > 0) {
                write(stderr_file, buffer, bytes_read);
                if (verbose_) {
                    write(STDERR_FILENO, buffer, bytes_read);
                }
            }
            close(stderr_file);
        }
        
        // Close the pipe file descriptors
        close(stdout_fd);
        close(stderr_fd);
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
        
        // Store both PIDs with descriptions
        pids_.push_back({splitter_pid, "Splitter"});
        pids_.push_back({stats_pid, "Stats"});
        
        // Wait for the splitter process to complete (we need its output for the next stage)
        int splitter_status;
        waitpid(splitter_pid, &splitter_status, 0);
        
        if (WEXITSTATUS(splitter_status) != 0) {
            std::cerr << "Error: Splitter process failed with status " << WEXITSTATUS(splitter_status) << std::endl;
            return false;
        }
        
        // Remove splitter PID from the pending list as we already waited for it
        auto it = std::find_if(pids_.begin(), pids_.end(), 
                             [splitter_pid](const auto& pair) { return pair.first == splitter_pid; });
        if (it != pids_.end()) {
            pids_.erase(it);
        }
        
        // Define paths to output files from the splitter
        clustered_edges_ = temp_dir_ + "/non_singleton_edges.tsv";
        clustered_clusters_ = temp_dir_ + "/non_singleton_clusters.tsv";
        unclustered_edges_ = temp_dir_ + "/singleton_edges.tsv";
        unclustered_clusters_ = temp_dir_ + "/singleton_clusters.tsv";
        
        // Check if output files exist and are not empty
        if (!checkFileExistsAndNotEmpty(clustered_edges_) || 
            !checkFileExistsAndNotEmpty(clustered_clusters_) || 
            !checkFileExistsAndNotEmpty(unclustered_edges_) || 
            !checkFileExistsAndNotEmpty(unclustered_clusters_)) {
            std::cerr << "Error: Expected output files from splitter not found or empty" << std::endl;
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
        std::string sbm_script = "extlib/gen_SBM.py";
        
        std::string sbm_clustered_command = "python3 " + sbm_script;
        sbm_clustered_command += " -f " + clustered_edges_;
        sbm_clustered_command += " -c " + clustered_clusters_;
        sbm_clustered_command += " -o " + temp_dir_ + "/clustered_sbm";
        
        std::string sbm_unclustered_command = "python3 " + sbm_script;
        sbm_unclustered_command += " -f " + unclustered_edges_;
        sbm_unclustered_command += " -c " + unclustered_clusters_;
        sbm_unclustered_command += " -o " + temp_dir_ + "/unclustered_sbm";

        if (verbose_) {
            sbm_clustered_command += " -v";
            sbm_unclustered_command += " -v";
        }
        
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
        
        // Store PIDs with descriptions
        pids_.push_back({clustered_sbm_pid, "SBM-Clustered"});
        pids_.push_back({unclustered_sbm_pid, "SBM-Unclustered"});
        
        return true;
    }
    
    /**
     * @brief Wait for all remaining child processes to complete
     * 
     * @return bool true if all processes completed successfully
     */
    bool waitForAllProcesses() {
        if (verbose_) {
            std::cout << "Waiting for " << pids_.size() << " remaining processes to complete..." << std::endl;
        }
        
        int status;
        bool all_success = true;
        
        for (const auto& [pid, description] : pids_) {
            if (verbose_) {
                std::cout << "Waiting for process " << pid << " (" << description << ")..." << std::endl;
            }
            
            // Wait for the process to finish
            waitpid(pid, &status, 0);
            
            // Find and process any output pipes for this process
            auto it = std::find_if(output_pipes_.begin(), output_pipes_.end(),
                                [pid](const auto& pipes) { return pipes.pid == pid; });
            
            if (it != output_pipes_.end()) {
                processCommandOutput(it->pid, it->description, it->stdout_fd, it->stderr_fd);
                output_pipes_.erase(it);
            }
            
            if (WEXITSTATUS(status) != 0) {
                std::cerr << "Error: Process " << pid << " (" << description << ") failed with status " 
                          << WEXITSTATUS(status) << std::endl;
                all_success = false;
                // Continue waiting for other processes even if one fails
            } else {
                if (verbose_) {
                    std::cout << "Process " << pid << " (" << description << ") completed successfully" << std::endl;
                }
                
                // For SBM processes, verify output files immediately
                if (description == "SBM-Clustered") {
                    std::string output_file = temp_dir_ + "/clustered_sbm/syn_sbm.tsv";
                    if (!checkFileExistsAndNotEmpty(output_file)) {
                        std::cerr << "Error: SBM-Clustered process completed but did not generate output file: " 
                                  << output_file << std::endl;
                        all_success = false;
                    }
                } else if (description == "SBM-Unclustered") {
                    std::string output_file = temp_dir_ + "/unclustered_sbm/syn_sbm.tsv";
                    if (!checkFileExistsAndNotEmpty(output_file)) {
                        std::cerr << "Error: SBM-Unclustered process completed but did not generate output file: " 
                                  << output_file << std::endl;
                        all_success = false;
                    }
                }
            }
        }
        
        // Clear the PID list
        pids_.clear();
        
        return all_success;
    }
    
    /**
     * @brief Check if a file exists and is not empty
     * 
     * @param filepath Path to the file to check
     * @return bool true if file exists and is not empty
     */
    bool checkFileExistsAndNotEmpty(const std::string& filepath) {
        if (!std::filesystem::exists(filepath)) {
            if (verbose_) {
                std::cerr << "File does not exist: " << filepath << std::endl;
            }
            return false;
        }
        
        if (std::filesystem::file_size(filepath) == 0) {
            if (verbose_) {
                std::cerr << "File is empty: " << filepath << std::endl;
            }
            return false;
        }
        
        return true;
    }
    
    /**
     * @brief Verify that all expected output files exist
     * 
     * @return bool true if all files exist
     */
    bool verifyOutputFiles() {
        std::vector<std::string> required_files = {
            temp_dir_ + "/clustered_stats.csv",
            temp_dir_ + "/clustered_sbm/syn_sbm.tsv",
            temp_dir_ + "/unclustered_sbm/syn_sbm.tsv"
        };
        
        bool all_files_exist = true;
        
        for (const auto& file : required_files) {
            if (!checkFileExistsAndNotEmpty(file)) {
                std::cerr << "Error: Required output file not found or empty: " << file << std::endl;
                all_files_exist = false;
                // Continue checking other files
            }
        }
        
        return all_files_exist;
    }
};
