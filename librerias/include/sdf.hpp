#pragma once

#include <functional>
#include <vector>

namespace lbm {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

namespace sdf {

using Function = std::function<double(double, double)>;

[[nodiscard]] double circle(double x, double y, double cx, double cy, double radius);
[[nodiscard]] double box(double x,
                         double y,
                         double cx,
                         double cy,
                         double half_width,
                         double half_height,
                         double angle_radians = 0.0);
[[nodiscard]] double polygon(double x, double y, const std::vector<Vec2>& vertices);

[[nodiscard]] double unite(double a, double b);
[[nodiscard]] double intersect(double a, double b);
[[nodiscard]] double subtract(double a, double b);

[[nodiscard]] Vec2 gradient(const Function& phi, double x, double y, double h = 1.0);
[[nodiscard]] Vec2 normal(const Function& phi, double x, double y, double h = 1.0);

} // namespace sdf
} // namespace lbm
