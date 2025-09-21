import pandas as pd
from sys import argv

def normalize_edges(edges_df):
    """
    Normalize edges for undirected graph by ensuring source <= target
    This prevents counting (A,B) and (B,A) as different edges
    """
    normalized = edges_df.copy()
    # Swap source and target where source > target
    mask = normalized['source'] > normalized['target']
    normalized.loc[mask, ['source', 'target']] = normalized.loc[mask, ['target', 'source']].values
    return normalized.drop_duplicates()

# Read the edge files
if len(argv) != 4:
    print("Usage: python sanity_check.py <current_edges.tsv> <new_edges.tsv> <output_file.tsv>")
    exit(1)

current_edges = pd.read_csv(argv[1], sep='\t', header=None, names=['source', 'target'])
new_edges = pd.read_csv(argv[2], sep='\t', header=None, names=['source', 'target'])

# Normalize both edge sets (for undirected graphs)
current_edges_normalized = normalize_edges(current_edges)
new_edges_normalized = normalize_edges(new_edges)

print(f"Current edges (normalized): {current_edges_normalized.shape[0]}")
print(f"New edges (normalized): {new_edges_normalized.shape[0]}")

# Find edges that are in new_edges but not in current_edges
current_edge_set = set(zip(current_edges_normalized['source'], current_edges_normalized['target']))
new_edge_set = set(zip(new_edges_normalized['source'], new_edges_normalized['target']))

# Find truly new edges
truly_new_edges = new_edge_set - current_edge_set
print(f"\nTruly new edges: {len(truly_new_edges)}")
if truly_new_edges:
    print(f"Adding {len(truly_new_edges)} new edges.")
    
    # Append truly new edges to current edges
    new_edges_to_add = pd.DataFrame(list(truly_new_edges), columns=['source', 'target'])

    # Concatenate with current edges
    updated_edges = pd.concat([current_edges_normalized, new_edges_to_add], ignore_index=True)

    updated_edges.to_csv(argv[3], sep='\t', header=False, index=False)
else:
    print("All edges in new edges are present in current edges.")
