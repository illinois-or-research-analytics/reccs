# RECCS++

**Scalable REalistic Cluster Connectivity Simulator for synthetic network generation**

## 🛠️ Requirements

- **cmake** (3.26.5+)
- **gcc/g++** (11.4.1+)
- **python** (3.13.2+) *for graph-tool*
- **conda** (24.9.2+) *for graph-tool*
- **openmp** (201511+)

## 🚀 Installation

### 1. Install Dependencies

Set up and load the `graph-tool` conda environment:

```bash
source install.sh
```

### 2. Build RECCS

Compile the project:

```bash
mkdir build
cd build
cmake .. && make
```

### 3. Verify Installation

Test that RECCS built successfully:

```bash
./eval/eval_pipeline_mini.sh
```

### 4. Initiate Pre-Commit Hook

Run this command in the terminal with your conda env active:

```bash
pre-commit install
```

This will run automated tests before you commit anychanges to Python or C++ files

## 📖 Usage

RECCS++ supports two modes: **Normal mode** (runs full pipeline) and **Checkpoint mode** (skips orchestrator using pre-generated files).

### Quick Start

**Normal mode** (recommended for first-time users):

```bash
./reccs graph.tsv -c clusters.tsv -v
```

**Checkpoint mode** (for resuming from intermediate files):

```bash
./reccs --checkpoint -c clusters.tsv \
  --clustered-sbm clustered.tsv --unclustered-sbm unclustered.tsv \
  --requirements requirements.csv --degseq degseq.json -v
```

### Command Line Reference

```text
Usage: ./reccs <edgelist.tsv> [options]
       ./reccs --checkpoint [checkpoint_options]

Note: In normal mode, <edgelist.tsv> must be the first argument.
      In checkpoint mode, --checkpoint must be the first argument.

Common options:
  -c <clusters.tsv> Load clusters from TSV file (required)
  -t <num_threads>  Number of threads to use (default: hardware concurrency)
  -v                Verbose mode: print detailed progress information
  -o <output_file>  Output file (default: 'output.tsv')
  -h, --help        Show this help message and exit

Normal mode specific options:
  <edgelist.tsv>                   Input graph edgelist file

Checkpoint mode specific options:
  --checkpoint                     Enable checkpoint mode (skip orchestrator)
  --clustered-sbm <path>          Path to clustered SBM graph file
  --unclustered-sbm <path>        Path to unclustered SBM graph file
  --requirements <path>           Path to requirements CSV file
  --degseq <path>                 Path to degree sequence JSON file
```

### Modes Explained

#### Normal Mode

Runs the complete RECCS++ pipeline starting from a graph edgelist and clustering file. The orchestrator will generate all intermediate files and process the graph through the full workflow.

**Use when:**

- Running RECCS++ for the first time
- You have a raw graph and clustering that need processing
- You want the complete pipeline execution

#### Checkpoint Mode

Skips the orchestrator and loads pre-generated intermediate files directly. This is useful for:

**Use when:**

- Resuming interrupted processing
- Re-running analysis with different parameters on the same intermediate data
- You already have the required intermediate files from a previous run

**Required checkpoint files:**

- Clustered SBM graph file (`--clustered-sbm`)
- Unclustered SBM graph file (`--unclustered-sbm`)
- Requirements CSV file (`--requirements`)
- Degree sequence JSON file (`--degseq`)

### Examples

**Complete pipeline with verbose output:**

```bash
./reccs input/network.tsv -c input/clusters.tsv -t 8 -v -o results.tsv
```

**Resume from checkpoint files:**

```bash
./reccs --checkpoint -c input/clusters.tsv \
  --clustered-sbm temp123/clustered_sbm/syn_sbm.tsv \
  --unclustered-sbm temp123/unclustered_sbm/syn_sbm.tsv \
  --requirements temp123/clustered_stats.csv \
  --degseq temp123/clustered_stats_degree_sequences.json \
  -v -o final_output.tsv
```

**Multi-threaded processing:**

```bash
./reccs graph.tsv -c clusters.tsv -t 16 -o output.tsv
```
