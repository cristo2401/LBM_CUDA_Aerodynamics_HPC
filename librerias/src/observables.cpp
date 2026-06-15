#include "observables.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace lbm::observables {

Sample compute(const Field& field, const SimulationParams& params, std::size_t step)
{
    Sample sample;
    sample.step = step;
    sample.time = static_cast<double>(step) * params.dt;
    sample.min_rho = std::numeric_limits<double>::infinity();
    sample.max_rho = -std::numeric_limits<double>::infinity();

    for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
        const double rho = field.rho[cell];
        const double ux = field.ux[cell];
        const double uy = field.uy[cell];

        if (!std::isfinite(rho) || !std::isfinite(ux) || !std::isfinite(uy)) {
            sample.finite = false;
        }

        sample.total_mass += rho;
        if (field.solid[cell]) {
            ++sample.solid_cells;
            continue;
        }

        ++sample.fluid_cells;
        sample.fluid_mass += rho;
        sample.min_rho = std::min(sample.min_rho, rho);
        sample.max_rho = std::max(sample.max_rho, rho);

        const double u2 = ux * ux + uy * uy;
        sample.max_speed = std::max(sample.max_speed, std::sqrt(u2));
        sample.kinetic_energy += 0.5 * rho * u2;
    }

    if (sample.fluid_cells > 0) {
        sample.mean_rho = sample.fluid_mass / static_cast<double>(sample.fluid_cells);
    } else {
        sample.min_rho = 0.0;
        sample.max_rho = 0.0;
        sample.finite = false;
    }

    return sample;
}

} // namespace lbm::observables
