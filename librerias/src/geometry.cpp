#include "geometry.hpp"

#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lbm::geometry {

namespace {

Vec2 rotate_about_center(double local_x,
                         double local_y,
                         double center_x,
                         double center_y,
                         double angle_radians)
{
    const double c = std::cos(angle_radians);
    const double s = std::sin(angle_radians);
    return {
        center_x + c * local_x - s * local_y,
        center_y + s * local_x + c * local_y
    };
}

void camber_line(double x_fraction,
                 double max_camber,
                 double camber_position,
                 double& yc,
                 double& dyc_dx)
{
    if (max_camber <= 0.0) {
        yc = 0.0;
        dyc_dx = 0.0;
        return;
    }

    const double x = std::clamp(x_fraction, 0.0, 1.0);
    const double p = camber_position;
    const double m = max_camber;

    if (x < p) {
        yc = m / (p * p) * (2.0 * p * x - x * x);
        dyc_dx = 2.0 * m / (p * p) * (p - x);
    } else {
        const double one_minus_p = 1.0 - p;
        yc = m / (one_minus_p * one_minus_p) *
             ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
        dyc_dx = 2.0 * m / (one_minus_p * one_minus_p) * (p - x);
    }
}

double thickness_distribution(double x_fraction, double max_thickness)
{
    const double x = std::clamp(x_fraction, 0.0, 1.0);
    return 5.0 * max_thickness *
           (0.2969 * std::sqrt(x) -
            0.1260 * x -
            0.3516 * x * x +
            0.2843 * x * x * x -
            0.1036 * x * x * x * x);
}

std::vector<Vec2> make_airfoil_polygon(const GeometryConfig& config)
{
    const std::size_t n = std::max<std::size_t>(config.airfoil_samples, 16);
    std::vector<Vec2> upper;
    std::vector<Vec2> lower;
    upper.reserve(n);
    lower.reserve(n);

    // For airfoils, positive angle_degrees means positive angle of attack:
    // nose up relative to a left-to-right free stream. In screen/world
    // coordinates that is a clockwise rotation of the chord line.
    const double angle = -config.angle_radians;
    const double pi = std::acos(-1.0);

    for (std::size_t i = 0; i < n; ++i) {
        const double beta = pi * static_cast<double>(i) / static_cast<double>(n - 1);
        const double x_frac = 0.5 * (1.0 - std::cos(beta));

        double yc = 0.0;
        double dyc_dx = 0.0;
        camber_line(x_frac,
                    config.airfoil_max_camber,
                    config.airfoil_camber_position,
                    yc,
                    dyc_dx);

        const double yt = thickness_distribution(x_frac, config.airfoil_max_thickness);
        const double theta = std::atan(dyc_dx);
        const double chord = config.airfoil_chord;

        const double x_base = (x_frac - 0.5) * chord;
        const double y_base = yc * chord;
        const double thickness = yt * chord;

        const double xu = x_base - thickness * std::sin(theta);
        const double yu = y_base + thickness * std::cos(theta);
        const double xl = x_base + thickness * std::sin(theta);
        const double yl = y_base - thickness * std::cos(theta);

        upper.push_back(rotate_about_center(xu,
                                            yu,
                                            config.center_x,
                                            config.center_y,
                                            angle));
        lower.push_back(rotate_about_center(xl,
                                            yl,
                                            config.center_x,
                                            config.center_y,
                                            angle));
    }

    std::vector<Vec2> polygon;
    polygon.reserve(upper.size() + lower.size());
    polygon.insert(polygon.end(), upper.begin(), upper.end());
    for (std::size_t i = lower.size() - 1; i > 0; --i) {
        polygon.push_back(lower[i]);
    }
    return polygon;
}

std::vector<Vec2> load_sdf_polygon(const GeometryConfig& config)
{
    std::ifstream in(config.sdf_file);
    if (!in) {
        throw std::runtime_error("could not open SDF polygon file: " + config.sdf_file);
    }

    std::vector<Vec2> vertices;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::istringstream stream(line);
        std::string tag_or_x;
        if (!(stream >> tag_or_x)) {
            continue;
        }

        double local_x = 0.0;
        double local_y = 0.0;
        if (tag_or_x == "v" || tag_or_x == "V") {
            if (!(stream >> local_x >> local_y)) {
                throw std::runtime_error("bad vertex in " + config.sdf_file +
                                         " at line " + std::to_string(line_number));
            }
        } else {
            std::size_t consumed = 0;
            local_x = std::stod(tag_or_x, &consumed);
            if (consumed != tag_or_x.size() || !(stream >> local_y)) {
                throw std::runtime_error("bad vertex in " + config.sdf_file +
                                         " at line " + std::to_string(line_number));
            }
        }

        const double angle = -config.angle_radians;
        vertices.push_back(rotate_about_center(local_x * config.sdf_scale,
                                               local_y * config.sdf_scale,
                                               config.center_x,
                                               config.center_y,
                                               angle));
    }

    if (vertices.size() < 3) {
        throw std::runtime_error("SDF polygon file needs at least three vertices: " +
                                 config.sdf_file);
    }
    return vertices;
}

} // namespace

GeometryConfig from_params(const SimulationParams& params)
{
    GeometryConfig config;
    config.type = params.obstacle.type;
    config.center_x = params.obstacle_center_x();
    config.center_y = params.obstacle_center_y();
    config.radius = params.obstacle.radius;
    config.half_width = params.obstacle.half_width;
    config.half_height = params.obstacle.half_height;
    config.angle_radians = utils::deg_to_rad(params.obstacle.angle_degrees);
    config.airfoil_chord = params.obstacle.airfoil_chord;
    config.airfoil_max_camber = params.obstacle.airfoil_max_camber;
    config.airfoil_camber_position = params.obstacle.airfoil_camber_position;
    config.airfoil_max_thickness = params.obstacle.airfoil_max_thickness;
    config.airfoil_samples = params.obstacle.airfoil_samples;
    config.sdf_file = params.obstacle.sdf_file;
    config.sdf_scale = params.obstacle.sdf_scale;
    return config;
}

sdf::Function make_sdf(const GeometryConfig& config)
{
    if (config.type == ObstacleType::ImportedSdf) {
        const std::vector<Vec2> vertices = load_sdf_polygon(config);
        return [vertices](double x, double y) {
            return sdf::polygon(x, y, vertices);
        };
    }

    if (config.type == ObstacleType::Airfoil) {
        const std::vector<Vec2> vertices = make_airfoil_polygon(config);
        return [vertices](double x, double y) {
            return sdf::polygon(x, y, vertices);
        };
    }

    if (config.type == ObstacleType::Box) {
        return [config](double x, double y) {
            return sdf::box(x,
                            y,
                            config.center_x,
                            config.center_y,
                            config.half_width,
                            config.half_height,
                            config.angle_radians);
        };
    }

    return [config](double x, double y) {
        return sdf::circle(x, y, config.center_x, config.center_y, config.radius);
    };
}

void build(const GeometryConfig& config, Field& field)
{
    const auto phi_fn = make_sdf(config);
    for (std::size_t y = 0; y < field.ny(); ++y) {
        for (std::size_t x = 0; x < field.nx(); ++x) {
            const std::size_t cell = field.index(x, y);
            const double px = static_cast<double>(x);
            const double py = static_cast<double>(y);
            const double value = phi_fn(px, py);
            const Vec2 n = sdf::normal(phi_fn, px, py);

            field.phi[cell] = value;
            field.solid[cell] = value <= 0.0 ? 1U : 0U;
            field.normal_x[cell] = n.x;
            field.normal_y[cell] = n.y;
        }
    }
}

} // namespace lbm::geometry
