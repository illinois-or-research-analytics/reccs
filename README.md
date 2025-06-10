# RECCS++

**Scalable REalistic Cluster Connectivity Simulator for synthetic network generation**

## 🛠️ Requirements

- **cmake** (3.26.5+)
- **gcc/g++** (11.4.1+)
- **python** (3.13.2+) *for graph-tool*
- **conda** (24.9.2+) *for graph-tool*

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

## 📖 Usage

### Quick Start

```bash
./reccs -e <edgelist.tsv> -c <clustering.tsv> -v
```

### Command Line Options

```
Usage: ./reccs <edgelist.tsv> [options]

Options:
  -t <num_threads>  Number of threads (default: hardware concurrency)
  -v                Enable verbose output
  -c <clusters.tsv> Load clusters from TSV file
  -o <output_dir>   Output directory (default: 'output')
  -h, --help        Show help message
```

### Example

```bash
./reccs input/network.tsv -c input/clusters.tsv -o results/ -v
```
