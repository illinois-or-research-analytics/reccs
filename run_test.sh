#!/bin/bash
set -e  # fail if any step fails

# 1. Activate your conda env with graph-tool
source install.sh

# 2. Build C++ code
mkdir -p build
cd build
cmake .. && make
cd ..

# 3. Run your test suite
pytest ./test/unit_test.py
