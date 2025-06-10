#ifndef ENFORCE_MINCUT_H
#define ENFORCE_MINCUT_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

void enforce_mincut(const Graph& g, uint32_t min_cut_size) {
    // Get minimum cut using the configured algorithm
    MincutResult result = compute_mincut(g);
    if (result.get_light_partition().empty() || result.get_heavy_partition().empty()) {
        std::cerr << "Error: Mincut algorithm returned empty partitions." << std::endl;
        return;
    }
    
    if (result.get_cut_size() < min_cut_size) {
        // Output the partitions
        std::cerr << "Light partition: ";
        for (uint32_t node : result.get_light_partition()) {
            std::cerr << node << " ";
        }
        std::cerr << "\nHeavy partition: ";
        for (uint32_t node : result.get_heavy_partition()) {
            std::cerr << node << " ";
        }
        std::cerr << "\nCut size: " << result.get_cut_size() << std::endl;
        std::cerr << "Required minimum cut size: " << min_cut_size << std::endl;
    }
}

#endif // ENFORCE_MINCUT_H
