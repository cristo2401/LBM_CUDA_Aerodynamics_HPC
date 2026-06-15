#include "boundary.hpp"

#include "lattice.hpp"
#include "lbm_utils.hpp"

namespace lbm::boundary {

namespace {

void copy_cell(Field& field,
               std::size_t dst_x,
               std::size_t dst_y,
               std::size_t src_x,
               std::size_t src_y)
{
    const std::size_t dst = field.index(dst_x, dst_y);
    const std::size_t src = field.index(src_x, src_y);
    for (int q = 0; q < D2Q9::q; ++q) {
        field.distribution(field.f, q, dst) = field.distribution(field.f, q, src);
    }
    field.rho[dst] = field.rho[src];
    field.ux[dst] = field.ux[src];
    field.uy[dst] = field.uy[src];
}

void apply_channel_walls(Field& field, const SimulationParams& params)
{
    const std::size_t top_y = field.ny() - 1;
    for (std::size_t x = 0; x < field.nx(); ++x) {
        const std::size_t bottom = field.index(x, 0);
        const std::size_t top = field.index(x, top_y);

        if (!field.solid[bottom]) {
            field.distribution(field.f, 2, bottom) = field.distribution(field.f, 4, bottom);
            field.distribution(field.f, 5, bottom) = field.distribution(field.f, 7, bottom);
            field.distribution(field.f, 6, bottom) = field.distribution(field.f, 8, bottom);
            field.rho[bottom] = params.rho0;
            field.ux[bottom] = 0.0;
            field.uy[bottom] = 0.0;
        }

        if (!field.solid[top]) {
            field.distribution(field.f, 4, top) = field.distribution(field.f, 2, top);
            field.distribution(field.f, 7, top) = field.distribution(field.f, 5, top);
            field.distribution(field.f, 8, top) = field.distribution(field.f, 6, top);
            field.rho[top] = params.rho0;
            field.ux[top] = 0.0;
            field.uy[top] = 0.0;
        }
    }
}

void apply_open_farfield(Field& field)
{
    const std::size_t top_y = field.ny() - 1;
    for (std::size_t x = 0; x < field.nx(); ++x) {
        if (!field.solid[field.index(x, 0)] && !field.solid[field.index(x, 1)]) {
            copy_cell(field, x, 0, x, 1);
        }
        if (!field.solid[field.index(x, top_y)] && !field.solid[field.index(x, top_y - 1)]) {
            copy_cell(field, x, top_y, x, top_y - 1);
        }
    }
}

} // namespace

void set_equilibrium_cell(Field& field,
                          std::size_t x,
                          std::size_t y,
                          double rho,
                          double ux,
                          double uy)
{
    const std::size_t cell = field.index(x, y);
    for (int q = 0; q < D2Q9::q; ++q) {
        field.distribution(field.f, q, cell) = lbm_utils::equilibrium(q, rho, ux, uy);
    }
    field.rho[cell] = rho;
    field.ux[cell] = ux;
    field.uy[cell] = uy;
}

void apply_external_boundaries(Field& field,
                               const SimulationParams& params,
                               double inlet_velocity)
{
    if (field.nx() < 2 || field.ny() < 2) {
        return;
    }

    if (params.boundary_mode == BoundaryMode::Channel) {
        apply_channel_walls(field, params);
    } else {
        apply_open_farfield(field);
    }

    const std::size_t right_x = field.nx() - 1;
    for (std::size_t y = 0; y < field.ny(); ++y) {
        if (!field.solid[field.index(0, y)]) {
            set_equilibrium_cell(field, 0, y, params.rho0, inlet_velocity, 0.0);
        }

        if (!field.solid[field.index(right_x, y)] && !field.solid[field.index(right_x - 1, y)]) {
            copy_cell(field, right_x, y, right_x - 1, y);
        } else if (!field.solid[field.index(right_x, y)]) {
            set_equilibrium_cell(field, right_x, y, params.rho0, inlet_velocity, 0.0);
        }
    }
}

} // namespace lbm::boundary
