# RECCS++

Scalable REalistic Cluster Connectivity Simulator for synthetic network generation

## Requirements

- `cmake` (3.26.5)
- `gcc/g++` (11.4.1)
- `python` (3.13.2) [for graph-tool]
- `conda` (24.9.2) [for graph-tool]

## Setting Up

### Installing Dependencies

Simply run the following to setup (if not set up already) and load the `graph-tool` conda environment and install dependencies:

```bash
source install.sh
```

### Building RECCS

Then run the following to build the project:

```bash
mkdir build
cd build
cmake .. && make
```

### Verification

To verify that RECCS built successfully, simply run

```bash
./eval/eval_pipeline_mini.sh
```

## Running RECCS

### Quick Start

Run using the following command:

```bash
./reccs -e <edgelist tsv> -c <clustering tsv> -v
```

### Arguments

```text
Usage: "./reccs <edgelist.tsv> [options]"
Options:
  -t <num_threads>  Number of threads to use (default: hardware concurrency)
  -v                Verbose mode: print detailed progress information
  -c <clusters.tsv> Load clusters from TSV file
  -o <output_dir>   Output tsv edgelist with added edges (default: 'output')
  -h, --help        Show this help message and exit
```
