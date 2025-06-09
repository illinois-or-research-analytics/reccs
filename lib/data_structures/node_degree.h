#ifndef NODE_DEGREE_H
#define NODE_DEGREE_H

// Node with degree for heap
struct NodeDegree {
    uint32_t node;
    uint32_t degree;
    
    // For min heap
    bool operator>(const NodeDegree& other) const {
        return degree > other.degree;
    }

    // For max heap
    bool operator<(const NodeDegree& other) const {
        return degree < other.degree;
    }
};

#endif // NODE_DEGREE_H
