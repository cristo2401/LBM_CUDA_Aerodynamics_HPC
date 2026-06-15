#pragma once

#include "field.hpp"
#include "geometry.hpp"
#include "params.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace lbm::case_builder {

struct Case {
    Field field;
    geometry::GeometryConfig geometry;
    std::size_t solid_cells = 0;
    std::size_t fluid_cells = 0;
    std::vector<std::string> warnings;
};

[[nodiscard]] Case build(const SimulationParams& params);

} // namespace lbm::case_builder
