# Check if gt environment exists, if not create it
if ! conda env list | grep "gt$"; then
    conda create --name gt -c conda-forge graph-tool
fi


# Install python dependencies
conda activate gt
pip install -r requirements.txt

export CMAKE_PREFIX_PATH=$HOME/local:$CMAKE_PREFIX_PATH
export igraph_DIR=$HOME/local/lib/cmake/igraph
