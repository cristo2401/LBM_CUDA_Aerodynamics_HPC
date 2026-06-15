#include "lbm_utils.hpp"

#include "lattice.hpp"

#include <algorithm>
#include <cmath>

namespace lbm::lbm_utils {

double equilibrium(int direction, double rho, double ux, double uy)
{
    const double cu = static_cast<double>(D2Q9::cx[direction]) * ux +
                      static_cast<double>(D2Q9::cy[direction]) * uy;
    const double u2 = ux * ux + uy * uy;
    return D2Q9::w[direction] * rho *
           (1.0 + D2Q9::inv_cs2 * cu + D2Q9::inv_2cs4 * cu * cu - 1.5 * u2);
}

double pressure(double rho)
{
    return D2Q9::cs2 * rho;
}

double speed(double ux, double uy)
{
    return std::sqrt(ux * ux + uy * uy);
}

void compute_vorticity(const Field& field, std::vector<float>& out)
{
    out.assign(field.cell_count(), 0.0F);
    if (field.nx() == 0 || field.ny() == 0) {
        return;
    }

    for (std::size_t y = 0; y < field.ny(); ++y) {
        const std::size_t y0 = y == 0 ? y : y - 1;
        const std::size_t y1 = std::min(y + 1, field.ny() - 1);
        const double inv_dy = y1 == y0 ? 0.0 : 1.0 / static_cast<double>(y1 - y0);

        for (std::size_t x = 0; x < field.nx(); ++x) {
            const std::size_t x0 = x == 0 ? x : x - 1;
            const std::size_t x1 = std::min(x + 1, field.nx() - 1);
            const double inv_dx = x1 == x0 ? 0.0 : 1.0 / static_cast<double>(x1 - x0);

            const std::size_t cell = field.index(x, y);
            const double duy_dx = (field.uy[field.index(x1, y)] - field.uy[field.index(x0, y)]) * inv_dx;
            const double dux_dy = (field.ux[field.index(x, y1)] - field.ux[field.index(x, y0)]) * inv_dy;
            out[cell] = static_cast<float>(duy_dx - dux_dy);
        }
    }
}

void compute_pressure(const Field& field, std::vector<float>& out)
{
    out.resize(field.cell_count());
    for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
        out[cell] = static_cast<float>(pressure(field.rho[cell]));
    }
}

} // namespace lbm::lbm_utils
