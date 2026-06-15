#include "case_builder.hpp"

#include <sstream>

namespace lbm::case_builder {

namespace {

void inspect_mask(Case& result)
{
    const Field& field = result.field;
    result.solid_cells = 0;
    result.fluid_cells = 0;

    for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
        if (field.solid[cell]) {
            ++result.solid_cells;
        } else {
            ++result.fluid_cells;
        }
    }
}

bool has_solid_on_boundary(const Field& field)
{
    if (field.nx() == 0 || field.ny() == 0) {
        return false;
    }

    const std::size_t right = field.nx() - 1;
    const std::size_t top = field.ny() - 1;

    for (std::size_t x = 0; x < field.nx(); ++x) {
        if (field.solid[field.index(x, 0)] || field.solid[field.index(x, top)]) {
            return true;
        }
    }

    for (std::size_t y = 0; y < field.ny(); ++y) {
        if (field.solid[field.index(0, y)] || field.solid[field.index(right, y)]) {
            return true;
        }
    }

    return false;
}

void add_warnings(const SimulationParams& params, Case& result)
{
    if (result.solid_cells == 0) {
        result.warnings.emplace_back("geometry produced zero solid cells; the run is just open flow.");
    }

    if (result.fluid_cells == 0) {
        result.warnings.emplace_back("geometry filled the whole domain; no fluid cells remain.");
    }

    if (has_solid_on_boundary(result.field)) {
        result.warnings.emplace_back("obstacle touches an external boundary; bounce-back and inlet/outlet will interact.");
    }

    const double diameter = params.reference_length();
    const double inlet_clearance = params.obstacle_center_x() - 0.5 * diameter;
    const bool slender_body = params.obstacle.type == ObstacleType::Airfoil ||
                              params.obstacle.type == ObstacleType::ImportedSdf;
    const double min_clearance = slender_body
        ? 0.55 * diameter
        : diameter;
    if (inlet_clearance < min_clearance) {
        std::ostringstream msg;
        msg << "obstacle is close to the inlet; clearance is " << inlet_clearance
            << " lattice units.";
        result.warnings.push_back(msg.str());
    }
}

} // namespace

Case build(const SimulationParams& params)
{
    Case result;
    result.field.resize(params.nx, params.ny);
    result.geometry = geometry::from_params(params);
    geometry::build(result.geometry, result.field);
    inspect_mask(result);
    add_warnings(params, result);
    return result;
}

} // namespace lbm::case_builder
