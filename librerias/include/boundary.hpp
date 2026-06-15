#pragma once

#include "field.hpp"
#include "params.hpp"

namespace lbm::boundary {

void apply_external_boundaries(Field& field,
                               const SimulationParams& params,
                               double inlet_velocity);
void set_equilibrium_cell(Field& field,
                          std::size_t x,
                          std::size_t y,
                          double rho,
                          double ux,
                          double uy);

} // namespace lbm::boundary
