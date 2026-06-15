#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace lbm {

enum class BoundaryMode {
    OpenFlow,
    Channel
};

enum class ObstacleType {
    Circle,
    Box,
    Airfoil,
    ImportedSdf
};

struct ObstacleParams {
    ObstacleType type = ObstacleType::Circle;
    double center_x = -1.0;
    double center_y = -1.0;
    double radius = 20.0;
    double half_width = 30.0;
    double half_height = 12.0;
    double angle_degrees = 0.0;
    double airfoil_chord = 80.0;
    double airfoil_max_camber = 0.02;
    double airfoil_camber_position = 0.40;
    double airfoil_max_thickness = 0.12;
    std::size_t airfoil_samples = 180;
    std::string sdf_file;
    double sdf_scale = 120.0;
};

struct SimulationParams {
    std::size_t nx = 420;
    std::size_t ny = 160;
    double dx = 1.0;
    double dt = 1.0;

    double rho0 = 1.0;
    double u_in = 0.06;
    double tau = 0.60;

    std::size_t steps = 2000;
    std::size_t frame_stride = 100;
    std::size_t force_stride = 10;
    std::size_t log_stride = 100;
    std::size_t inlet_ramp_steps = 100;

    BoundaryMode boundary_mode = BoundaryMode::OpenFlow;
    ObstacleParams obstacle;

    std::string case_name = "cylinder_d2q9";
    std::string output_dir = "resultados/default";
    bool write_initial_frame = true;
    bool start_from_rest = true;

    [[nodiscard]] double omega() const
    {
        return 1.0 / tau;
    }

    [[nodiscard]] double viscosity() const
    {
        return (tau - 0.5) / 3.0;
    }

    [[nodiscard]] double obstacle_center_x() const
    {
        if ((obstacle.type == ObstacleType::Airfoil ||
             obstacle.type == ObstacleType::ImportedSdf) &&
            obstacle.center_x < 0.0) {
            return 0.42 * static_cast<double>(nx);
        }
        return obstacle.center_x >= 0.0 ? obstacle.center_x : 0.25 * static_cast<double>(nx);
    }

    [[nodiscard]] double obstacle_center_y() const
    {
        return obstacle.center_y >= 0.0 ? obstacle.center_y : 0.5 * static_cast<double>(ny);
    }

    [[nodiscard]] double reference_length() const
    {
        if (obstacle.type == ObstacleType::Circle) {
            return 2.0 * obstacle.radius;
        }
        if (obstacle.type == ObstacleType::Airfoil) {
            return obstacle.airfoil_chord;
        }
        if (obstacle.type == ObstacleType::ImportedSdf) {
            return obstacle.sdf_scale;
        }
        return 2.0 * std::max(obstacle.half_width, obstacle.half_height);
    }

    [[nodiscard]] double reynolds() const
    {
        const double nu = viscosity();
        return nu > 0.0 ? u_in * reference_length() / nu : 0.0;
    }

    [[nodiscard]] double initial_velocity() const
    {
        return start_from_rest ? 0.0 : u_in;
    }

    [[nodiscard]] double inlet_velocity(std::size_t step) const
    {
        if (inlet_ramp_steps == 0) {
            return u_in;
        }

        const double t_raw = static_cast<double>(step) / static_cast<double>(inlet_ramp_steps);
        const double t = std::clamp(t_raw, 0.0, 1.0);
        const double smooth = t * t * (3.0 - 2.0 * t);
        return u_in * smooth;
    }

    void validate() const
    {
        if (nx < 16 || ny < 16) {
            throw std::invalid_argument("nx and ny must be at least 16.");
        }
        if (tau <= 0.5) {
            throw std::invalid_argument("tau must be greater than 0.5 for positive viscosity.");
        }
        if (rho0 <= 0.0) {
            throw std::invalid_argument("rho0 must be positive.");
        }
        if (u_in < 0.0) {
            throw std::invalid_argument("u_in must be non-negative.");
        }
        if (frame_stride == 0 || force_stride == 0 || log_stride == 0) {
            throw std::invalid_argument("strides must be greater than zero.");
        }
        if (obstacle.type == ObstacleType::Circle && obstacle.radius <= 0.0) {
            throw std::invalid_argument("circle radius must be positive.");
        }
        if (obstacle.type == ObstacleType::Box &&
            (obstacle.half_width <= 0.0 || obstacle.half_height <= 0.0)) {
            throw std::invalid_argument("box half extents must be positive.");
        }
        if (obstacle.type == ObstacleType::Airfoil) {
            if (obstacle.airfoil_chord <= 0.0) {
                throw std::invalid_argument("airfoil chord must be positive.");
            }
            if (obstacle.airfoil_max_camber < 0.0 || obstacle.airfoil_max_camber > 0.20) {
                throw std::invalid_argument("airfoil max camber must be in [0, 0.20].");
            }
            if (obstacle.airfoil_camber_position <= 0.0 ||
                obstacle.airfoil_camber_position >= 1.0) {
                throw std::invalid_argument("airfoil camber position must be in (0, 1).");
            }
            if (obstacle.airfoil_max_thickness <= 0.0 || obstacle.airfoil_max_thickness > 0.40) {
                throw std::invalid_argument("airfoil max thickness must be in (0, 0.40].");
            }
            if (obstacle.airfoil_samples < 16) {
                throw std::invalid_argument("airfoil samples must be at least 16.");
            }
        }
        if (obstacle.type == ObstacleType::ImportedSdf) {
            if (obstacle.sdf_file.empty()) {
                throw std::invalid_argument("imported SDF requires --sdf-file.");
            }
            if (obstacle.sdf_scale <= 0.0) {
                throw std::invalid_argument("sdf scale must be positive.");
            }
        }
    }
};

} // namespace lbm
