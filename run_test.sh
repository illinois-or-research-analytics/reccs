#!/bin/bash
set -e  # Stop on any error

echo "🔧 Loading required modules..."

# Load appropriate compiler, cmake, and conda modules
module load gcc/13.3.0
module load miniconda3/24.9.2
module load python/3.13.2
module load cmake   # If needed, or install via conda if not available

echo "📦 Activating Conda environment via install.sh..."
source install.sh

echo "🏗️ Building RECCS++..."
mkdir -p build
cd build
cmake .. && make
cd ..

pip install pytest

echo "✅ Running pytest..."
pytest ./test/unit_test.py

