#!/bin/bash

# Enable conda shell support
source "$(conda info --base)/etc/profile.d/conda.sh"

# Check if gt environment exists, if not create it
if ! conda env list | grep -q "gt$"; then
    conda create --yes --name gt -c conda-forge graph-tool
fi


# Install python dependencies
conda activate gt
pip install -r requirements.txt

export CMAKE_PREFIX_PATH=$HOME/local:$CMAKE_PREFIX_PATH
