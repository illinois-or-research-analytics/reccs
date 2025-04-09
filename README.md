# RECCS++

Scalable REalistic Cluster Connectivity Simulator for synthetic network generation

- [RECCS++](#reccs)
  - [Requirements](#requirements)
  - [Setting Up](#setting-up)
    - [Installing Dependencies](#installing-dependencies)
    - [Setting up `graph-tool`](#setting-up-graph-tool)
  - [Running RECCS](#running-reccs)

## Requirements

- `cmake` (3.26.5)
- `gcc/g++` (11.4.1)
- `python` (3.13.2) [for graph-tool]
- `conda` (24.9.2) [for graph-tool]

## Setting Up

### Installing Dependencies

Simply run the following to install dependencies:

```bash
source install.sh
```

Then run the following to build the project:

```bash
mkdir build
cd build
cmake ..
make
```

### Setting up `graph-tool`

**Note**: This step can be skipped if you'd rather use the native iGraph Stochastic Block Model (SBM).

```bash
conda create --name gt -c conda-forge graph-tool
conda activate gt
```

## Running RECCS

Run using the following command:

```bash
./reccs -e <edgelist tsv> -c <clustering tsv>
```

Or if you'd rather use iGraph's native SBM, run:

```bash
./reccs -e <edgelist tsv> -c <clustering tsv> -sbm igraph
```
