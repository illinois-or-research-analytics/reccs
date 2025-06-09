#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

// Viecut includes
#include "io/graph_io.h"
#include "algorithms/global_mincut/algorithms.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"