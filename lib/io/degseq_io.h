#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using json = nlohmann::json;

json load_degseq_json(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return json();
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        std::cerr << "Error: Could not parse JSON from file " << filename << " - " << e.what() << std::endl;
        return json();
    }

    file.close();
    return j;
}
