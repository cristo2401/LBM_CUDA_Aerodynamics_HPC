#include "lbm_solver.hpp"

#include "boundary.hpp"
#include "lattice.hpp"
#include "lbm_utils.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace lbm {

LBMSolver::LBMSolver(SimulationParams params, Field field)
    : params_(std::move(params))
    , field_(std::move(field))
{
    if (field_.nx() != params_.nx || field_.ny() != params_.ny) {
        throw std::invalid_argument("field dimensions do not match simulation parameters.");
    }
}

void LBMSolver::initialize()
{
    for (std::size_t y = 0; y < field_.ny(); ++y) {
        for (std::size_t x = 0; x < field_.nx(); ++x) {
            const std::size_t cell = field_.index(x, y);
            const double ux = field_.solid[cell] ? 0.0 : params_.initial_velocity();
            const double uy = 0.0;

            field_.rho[cell] = params_.rho0;
            field_.ux[cell] = ux;
            field_.uy[cell] = uy;

            for (int q = 0; q < D2Q9::q; ++q) {
                const double feq = lbm_utils::equilibrium(q, params_.rho0, ux, uy);
                field_.distribution(field_.f, q, cell) = feq;
                field_.distribution(field_.f_next, q, cell) = feq;
            }
        }
    }

    boundary::apply_external_boundaries(field_, params_, params_.inlet_velocity(0));
    compute_macroscopic();
}

void LBMSolver::compute_macroscopic()
{
    const int nx = static_cast<int>(field_.nx());
    const int ny = static_cast<int>(field_.ny());

#ifdef LBM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const std::size_t cell = field_.index(static_cast<std::size_t>(x),
                                                 static_cast<std::size_t>(y));
            if (field_.solid[cell]) {
                field_.rho[cell] = params_.rho0;
                field_.ux[cell] = 0.0;
                field_.uy[cell] = 0.0;
                continue;
            }

            double rho = 0.0;
            double momentum_x = 0.0;
            double momentum_y = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                const double fi = field_.distribution(field_.f, q, cell);
                rho += fi;
                momentum_x += fi * static_cast<double>(D2Q9::cx[q]);
                momentum_y += fi * static_cast<double>(D2Q9::cy[q]);
            }

            if (rho <= 1.0e-14) {
                rho = params_.rho0;
                momentum_x = 0.0;
                momentum_y = 0.0;
            }

            field_.rho[cell] = rho;
            field_.ux[cell] = momentum_x / rho;
            field_.uy[cell] = momentum_y / rho;
        }
    }
}

double LBMSolver::collide_value(int direction, std::size_t cell) const
{
    const double f_old = field_.distribution(field_.f, direction, cell);
    const double feq = lbm_utils::equilibrium(direction,
                                             field_.rho[cell],
                                             field_.ux[cell],
                                             field_.uy[cell]);
    return f_old - params_.omega() * (f_old - feq);
}

void LBMSolver::collide_and_stream()
{
    std::fill(field_.f_next.begin(), field_.f_next.end(), 0.0);

    const int nx = static_cast<int>(field_.nx());
    const int ny = static_cast<int>(field_.ny());

#ifdef LBM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const std::size_t cell = field_.index(static_cast<std::size_t>(x),
                                                 static_cast<std::size_t>(y));

            if (field_.solid[cell]) {
                for (int q = 0; q < D2Q9::q; ++q) {
                    field_.distribution(field_.f_next, q, cell) =
                        lbm_utils::equilibrium(q, params_.rho0, 0.0, 0.0);
                }
                continue;
            }

            for (int q = 0; q < D2Q9::q; ++q) {
                const int src_x = x - D2Q9::cx[q];
                const int src_y = y - D2Q9::cy[q];
                double streamed = field_.distribution(field_.f, q, cell);

                if (field_.in_bounds(src_x, src_y)) {
                    const std::size_t src = field_.index(static_cast<std::size_t>(src_x),
                                                        static_cast<std::size_t>(src_y));
                    if (field_.solid[src]) {
                        streamed = collide_value(D2Q9::opposite[q], cell);
                    } else {
                        streamed = collide_value(q, src);
                    }
                }

                field_.distribution(field_.f_next, q, cell) = streamed;
            }
        }
    }
}

void LBMSolver::step()
{
    compute_macroscopic();
    collide_and_stream();
    field_.f.swap(field_.f_next);
    boundary::apply_external_boundaries(field_, params_, params_.inlet_velocity(current_step_ + 1));
    ++current_step_;
    compute_macroscopic();
}

} // namespace lbm
