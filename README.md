# RECCS++
Scalable REalistic Cluster Connectivity Simulator for synthetic network generation

# Getting started
Simply run the following to install dependencies:

```
source install.sh
```

Then run the following to build the project:

```
mkdir build
cd build
cmake ..
make
```

Then run using the following command:
```
./reccs -e <edgelist tsv> -c <clustering tsv>
```
