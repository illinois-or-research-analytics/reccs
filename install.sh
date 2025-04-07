# Install script for the project

# Install submodules
git submodule update --init --recursive

# Begin: Install igraph
cd extlib/igraph

# If build directory does not exist, create it
if [ ! -d "build" ]; then
    mkdir build
fi
cd build

# If igraph is not installed, install it
if [ ! -d "lib" ]; then
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local
    make -j4
    make install
else
    echo "igraph is already installed."
fi
# End: Install igraph

cd ../../../

export CMAKE_PREFIX_PATH=$HOME/local:$CMAKE_PREFIX_PATH
export igraph_DIR=$HOME/local/lib/cmake/igraph
