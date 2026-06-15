#include "forces.hpp"

#include "lattice.hpp"

#include <cmath>

namespace lbm::forces {

ForceSample compute(const Field& field, const SimulationParams& params, std::size_t step)
{
    ForceSample sample;
    sample.step = step;
    sample.time = static_cast<double>(step) * params.dt;

    for (std::size_t y = 0; y < field.ny(); ++y) {
        for (std::size_t x = 0; x < field.nx(); ++x) {
            const std::size_t cell = field.index(x, y);
            if (field.solid[cell]) {
                continue;
            }

            for (int q = 1; q < D2Q9::q; ++q) {
                const int sx = static_cast<int>(x) + D2Q9::cx[q];
                const int sy = static_cast<int>(y) + D2Q9::cy[q];
                if (!field.in_bounds(sx, sy)) {
                    continue;
                }

                const std::size_t neighbor = field.index(static_cast<std::size_t>(sx),
                                                         static_cast<std::size_t>(sy));
                if (!field.solid[neighbor]) {
                    continue;
                }

                const double fi = field.distribution(field.f, q, cell);
                sample.fx += 2.0 * fi * static_cast<double>(D2Q9::cx[q]);
                sample.fy += 2.0 * fi * static_cast<double>(D2Q9::cy[q]);
            }
        }
    }

    const double denom = 0.5 * params.rho0 * params.u_in * params.u_in * params.reference_length();
    if (std::abs(denom) > 1.0e-14) {
        sample.cd = sample.fx / denom;
        sample.cl = sample.fy / denom;
    }

    return sample;
}

} // namespace lbm::forces
