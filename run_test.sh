#!/bin/bash
set -e

echo "Loading required modules"

module load gcc/13.3.0
module load miniconda3/24.9.2
module load python/3.13.2
module load openmpi

echo "Running install.sh"
source install.sh


echo "Building RECCS++..."
mkdir -p build
cd build
cmake .. && make
cd ..

pip install pytest

echo "Running pytest"
pytest ./test/unit_test.py

