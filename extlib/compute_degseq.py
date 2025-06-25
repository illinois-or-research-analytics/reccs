from sys import argv
import pandas as pd

edgelist_df = pd.read_csv(argv[1], sep="\t", header=None, names=["source", "target"])

# Get all nodes from both source and target columns
all_nodes = pd.concat([edgelist_df['source'], edgelist_df['target']])

# Count degree for each node
degree_counts = all_nodes.value_counts()

# Sort in descending order and convert to list with int values
degree_sequence = sorted([int(x) for x in degree_counts.values], reverse=True)

# Save the degree sequence to a json file containing a list of integers
pd.Series(degree_sequence).to_json(argv[2], orient='values')
