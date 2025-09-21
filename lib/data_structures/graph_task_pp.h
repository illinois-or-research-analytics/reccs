#ifndef GRAPH_TASK_PP_H
#define GRAPH_TASK_PP_H

#include "graph_task_base.h"

// TaskType enum for reccspp workflow (merged stages)
enum class TaskTypePP {
    CONNECTIVITY_ENFORCE = 0,  // Combined degree + connectivity enforcement
    WCC_STITCHING = 1,
    DEG_SEQ_MATCHING = 2
};

// GraphTask struct for reccspp workflow
struct GraphTaskPP : public GraphTaskBase {
    // Current task to perform
    TaskTypePP task_type;

    // Constructor
    GraphTaskPP(std::shared_ptr<Graph> g, TaskTypePP t, const std::string& cid, 
                uint32_t cidx, uint32_t min_deg = 1,
                std::shared_ptr<const std::vector<uint32_t>> deg_seq = nullptr)
        : GraphTaskBase(g, cid, cidx, min_deg, deg_seq), task_type(t) {}
};

// Helper functions for TaskTypePP
namespace TaskTypePPHelper {
    inline std::string get_task_name(TaskTypePP type) {
        switch (type) {
            case TaskTypePP::CONNECTIVITY_ENFORCE: return "connectivity_enforce";
            case TaskTypePP::WCC_STITCHING: return "wcc_stitching";
            case TaskTypePP::DEG_SEQ_MATCHING: return "deg_seq_matching";
            default: return "unknown";
        }
    }
}

#endif // GRAPH_TASK_PP_H
