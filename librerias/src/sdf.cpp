#include "sdf.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace lbm::sdf {

double circle(double x, double y, double cx, double cy, double radius)
{
    const double dx = x - cx;
    const double dy = y - cy;
    return std::sqrt(dx * dx + dy * dy) - radius;
}

double box(double x,
           double y,
           double cx,
           double cy,
           double half_width,
           double half_height,
           double angle_radians)
{
    const double dx = x - cx;
    const double dy = y - cy;
    const double c = std::cos(angle_radians);
    const double s = std::sin(angle_radians);

    const double local_x = c * dx + s * dy;
    const double local_y = -s * dx + c * dy;

    const double qx = std::abs(local_x) - half_width;
    const double qy = std::abs(local_y) - half_height;
    const double ox = std::max(qx, 0.0);
    const double oy = std::max(qy, 0.0);
    const double outside = std::sqrt(ox * ox + oy * oy);
    const double inside = std::min(std::max(qx, qy), 0.0);
    return outside + inside;
}

double polygon(double x, double y, const std::vector<Vec2>& vertices)
{
    if (vertices.size() < 3) {
        return std::numeric_limits<double>::infinity();
    }

    double min_distance2 = std::numeric_limits<double>::infinity();
    bool inside = false;

    for (std::size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++) {
        const Vec2 a = vertices[j];
        const Vec2 b = vertices[i];
        const double abx = b.x - a.x;
        const double aby = b.y - a.y;
        const double apx = x - a.x;
        const double apy = y - a.y;
        const double ab2 = abx * abx + aby * aby;
        const double t = ab2 > 0.0
            ? std::clamp((apx * abx + apy * aby) / ab2, 0.0, 1.0)
            : 0.0;
        const double closest_x = a.x + t * abx;
        const double closest_y = a.y + t * aby;
        const double dx = x - closest_x;
        const double dy = y - closest_y;
        min_distance2 = std::min(min_distance2, dx * dx + dy * dy);

        const bool crosses_y = (a.y > y) != (b.y > y);
        if (crosses_y) {
            const double x_intersection = (b.x - a.x) * (y - a.y) / (b.y - a.y) + a.x;
            if (x < x_intersection) {
                inside = !inside;
            }
        }
    }

    const double distance = std::sqrt(min_distance2);
    return inside ? -distance : distance;
}

double unite(double a, double b)
{
    return std::min(a, b);
}

double intersect(double a, double b)
{
    return std::max(a, b);
}

double subtract(double a, double b)
{
    return std::max(a, -b);
}

Vec2 gradient(const Function& phi, double x, double y, double h)
{
    const double inv = 0.5 / h;
    return {
        (phi(x + h, y) - phi(x - h, y)) * inv,
        (phi(x, y + h) - phi(x, y - h)) * inv
    };
}

Vec2 normal(const Function& phi, double x, double y, double h)
{
    Vec2 g = gradient(phi, x, y, h);
    const double n = std::sqrt(g.x * g.x + g.y * g.y);
    if (n <= 1.0e-14) {
        return {};
    }
    return {g.x / n, g.y / n};
}

} // namespace lbm::sdf
