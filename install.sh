# Install python dependencies
conda activate gt
pip install -r requirements.txt

export CMAKE_PREFIX_PATH=$HOME/local:$CMAKE_PREFIX_PATH
export igraph_DIR=$HOME/local/lib/cmake/igraph
