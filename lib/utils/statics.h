#ifndef STATICS_H
#define STATICS_H

#include <vector>
#include <memory>
#include <unordered_set>
#include <functional>
#include "../data_structures/graph.h"

namespace statics {
    // Empty degree sequence to use when no sequence is available
    inline std::shared_ptr<const std::vector<uint32_t>> empty_sequence = 
        std::make_shared<const std::vector<uint32_t>>();

    // Edge lookup utilities
    
    /**
     * Compute hash set of existing edges for O(1) lookup.
     * Each edge is stored once as (min_node << 32) | max_node.
     * 
     * @param g The graph to process
     * @return Unordered set containing encoded edges
     */
    inline std::unordered_set<uint64_t> compute_existing_edges(const Graph& g) {
        std::unordered_set<uint64_t> existing_edges;
        
        for (uint32_t u = 0; u < g.num_nodes; ++u) {
            for (uint32_t idx = g.row_ptr[u]; idx < g.row_ptr[u + 1]; ++idx) {
                uint32_t v = g.col_idx[idx];
                if (u < v) { // Store each edge once
                    existing_edges.insert((static_cast<uint64_t>(u) << 32) | v);
                }
            }
        }
        
        return existing_edges;
    }
    
    /**
     * Create a fast edge existence checker using precomputed edge set.
     * 
     * @param existing_edges Precomputed set from compute_existing_edges()
     * @return Lambda function for O(1) edge existence checking
     */
    inline auto create_edge_exists_checker(const std::unordered_set<uint64_t>& existing_edges) {
        return [&existing_edges](uint32_t u, uint32_t v) -> bool {
            if (u > v) std::swap(u, v);
            return existing_edges.count((static_cast<uint64_t>(u) << 32) | v) > 0;
        };
    }
    
    /**
     * Encode an edge as a 64-bit integer for use in hash sets.
     * Always ensures smaller node ID is in upper 32 bits.
     * 
     * @param u First node ID
     * @param v Second node ID
     * @return 64-bit encoded edge
     */
    inline uint64_t encode_edge(uint32_t u, uint32_t v) {
        if (u > v) std::swap(u, v);
        return (static_cast<uint64_t>(u) << 32) | v;
    }
    
    /**
     * Decode a 64-bit edge back to node pair.
     * 
     * @param encoded_edge 64-bit encoded edge
     * @return Pair of (smaller_node, larger_node)
     */
    inline std::pair<uint32_t, uint32_t> decode_edge(uint64_t encoded_edge) {
        uint32_t u = static_cast<uint32_t>(encoded_edge >> 32);
        uint32_t v = static_cast<uint32_t>(encoded_edge & 0xFFFFFFFF);
        return {u, v};
    }
}

#endif // STATICS_H
