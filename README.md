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

Simply run the following to setup (if not set up already) and load the `graph-tool` conda environment and install dependencies:

```bash
source install.sh
```

Then run the following to build the project:

```bash
mkdir build
cd build
cmake .. && make
```

## Running RECCS

### Quick Start

Run using the following command:

```bash
./reccs -e <edgelist tsv> -c <clustering tsv> -v
```
