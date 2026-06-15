#pragma once

#include "field.hpp"
#include "params.hpp"

#include <cstddef>

namespace lbm::observables {

struct Sample {
    std::size_t step = 0;
    double time = 0.0;
    double total_mass = 0.0;
    double fluid_mass = 0.0;
    double mean_rho = 0.0;
    double min_rho = 0.0;
    double max_rho = 0.0;
    double max_speed = 0.0;
    double kinetic_energy = 0.0;
    std::size_t solid_cells = 0;
    std::size_t fluid_cells = 0;
    bool finite = true;
};

[[nodiscard]] Sample compute(const Field& field,
                             const SimulationParams& params,
                             std::size_t step);

} // namespace lbm::observables
