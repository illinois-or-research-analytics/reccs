from sys import argv
import pandas as pd

reference_edgelist_df = pd.read_csv(argv[1], sep="\t", header=None, names=["source", "target"])
sbm_edgelist_df = pd.read_csv(argv[2], sep="\t", header=None, names=["source", "target"])

# Ensure integer type for source and target columns
reference_edgelist_df['source'] = reference_edgelist_df['source'].astype(int)
reference_edgelist_df['target'] = reference_edgelist_df['target'].astype(int)
sbm_edgelist_df['source'] = sbm_edgelist_df['source'].astype(int)
sbm_edgelist_df['target'] = sbm_edgelist_df['target'].astype(int)

# Get the degree counts for both reference and SBM edgelists (undirected)
reference_all_nodes = pd.concat([reference_edgelist_df['source'], reference_edgelist_df['target']])
sbm_all_nodes = pd.concat([sbm_edgelist_df['source'], sbm_edgelist_df['target']])
reference_degree_counts = reference_all_nodes.value_counts()
sbm_degree_counts = sbm_all_nodes.value_counts()

# Get the difference in degree counts for each node (set to reference degree if not present in sbm)
all_nodes = set(reference_degree_counts.index).union(set(sbm_degree_counts.index))
degree_diff = {
    node: (
        reference_degree_counts.get(node, 0) - sbm_degree_counts.get(node, 0)
        if node in sbm_degree_counts
        else reference_degree_counts.get(node, 0)
    )
    for node in all_nodes
}

# Filter out nodes with zero degree difference
degree_diff = {node: diff for node, diff in degree_diff.items() if diff != 0}

# Save as a node to degree difference mapping in a json file
pd.Series(degree_diff).to_json(argv[3], orient='index') 
