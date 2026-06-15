#pragma once

#include "field.hpp"
#include "params.hpp"

#include <cstddef>

namespace lbm::forces {

struct ForceSample {
    std::size_t step = 0;
    double time = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double cd = 0.0;
    double cl = 0.0;
};

[[nodiscard]] ForceSample compute(const Field& field,
                                  const SimulationParams& params,
                                  std::size_t step);

} // namespace lbm::forces
