#ifndef GRAPH_TASK_H
#define GRAPH_TASK_H

#include "graph_task_base.h"

// TaskType enum for original reccs workflow (3 separate stages)
enum class TaskType {
    MIN_DEG_ENFORCE = 0,
    CC_STITCHING = 1,
    WCC_STITCHING = 2
};

// GraphTask struct for original reccs workflow
struct GraphTask : public GraphTaskBase {
    // Current task to perform
    TaskType task_type;

    // Constructor
    GraphTask(std::shared_ptr<Graph> g, TaskType t, const std::string& cid, 
              uint32_t cidx, uint32_t min_deg = 1,
              std::shared_ptr<const std::vector<uint32_t>> deg_seq = nullptr)
        : GraphTaskBase(g, cid, cidx, min_deg, deg_seq), task_type(t) {}
};

// Helper functions for TaskType
namespace TaskTypeHelper {
    inline std::string get_task_name(TaskType type) {
        switch (type) {
            case TaskType::MIN_DEG_ENFORCE: return "min_deg_enforce";
            case TaskType::CC_STITCHING: return "cc_stitching";
            case TaskType::WCC_STITCHING: return "wcc_stitching";
            default: return "unknown";
        }
    }
}

#endif // GRAPH_TASK_H
