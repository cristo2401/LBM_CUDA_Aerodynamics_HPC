#pragma once

#include "field.hpp"

#include <vector>

namespace lbm::lbm_utils {

[[nodiscard]] double equilibrium(int direction, double rho, double ux, double uy);
[[nodiscard]] double pressure(double rho);
[[nodiscard]] double speed(double ux, double uy);

void compute_vorticity(const Field& field, std::vector<float>& out);
void compute_pressure(const Field& field, std::vector<float>& out);

} // namespace lbm::lbm_utils
