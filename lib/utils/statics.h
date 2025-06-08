#ifndef STATICS_H
#define STATICS_H

#include <vector>
#include <memory>

namespace statics {
    // Empty degree sequence to use when no sequence is available
    inline std::shared_ptr<const std::vector<uint32_t>> empty_sequence = 
        std::make_shared<const std::vector<uint32_t>>();
}

#endif // STATICS_H
