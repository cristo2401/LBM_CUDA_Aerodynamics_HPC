#pragma once

#include "field.hpp"
#include "params.hpp"
#include "sdf.hpp"

namespace lbm::geometry {

struct GeometryConfig {
    ObstacleType type = ObstacleType::Circle;
    double center_x = 0.0;
    double center_y = 0.0;
    double radius = 20.0;
    double half_width = 30.0;
    double half_height = 12.0;
    double angle_radians = 0.0;
    double airfoil_chord = 80.0;
    double airfoil_max_camber = 0.02;
    double airfoil_camber_position = 0.40;
    double airfoil_max_thickness = 0.12;
    std::size_t airfoil_samples = 180;
    std::string sdf_file;
    double sdf_scale = 120.0;
};

[[nodiscard]] GeometryConfig from_params(const SimulationParams& params);
[[nodiscard]] sdf::Function make_sdf(const GeometryConfig& config);
void build(const GeometryConfig& config, Field& field);

} // namespace lbm::geometry
