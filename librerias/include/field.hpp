#pragma once

#include "lattice.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace lbm {

class Field {
public:
    Field() = default;

    Field(std::size_t nx, std::size_t ny)
    {
        resize(nx, ny);
    }

    void resize(std::size_t nx, std::size_t ny)
    {
        nx_ = nx;
        ny_ = ny;
        const std::size_t cells = cell_count();
        f.assign(static_cast<std::size_t>(D2Q9::q) * cells, 0.0);
        f_next.assign(static_cast<std::size_t>(D2Q9::q) * cells, 0.0);
        rho.assign(cells, 0.0);
        ux.assign(cells, 0.0);
        uy.assign(cells, 0.0);
        solid.assign(cells, 0);
        phi.assign(cells, 0.0);
        normal_x.assign(cells, 0.0);
        normal_y.assign(cells, 0.0);
    }

    [[nodiscard]] std::size_t nx() const { return nx_; }
    [[nodiscard]] std::size_t ny() const { return ny_; }
    [[nodiscard]] std::size_t cell_count() const { return nx_ * ny_; }

    [[nodiscard]] bool in_bounds(int x, int y) const
    {
        return x >= 0 && y >= 0 &&
               x < static_cast<int>(nx_) &&
               y < static_cast<int>(ny_);
    }

    [[nodiscard]] std::size_t index(std::size_t x, std::size_t y) const
    {
        return x + nx_ * y;
    }

    [[nodiscard]] std::size_t dist_index(int direction, std::size_t cell) const
    {
        return static_cast<std::size_t>(direction) * cell_count() + cell;
    }

    [[nodiscard]] double distribution(const std::vector<double>& buffer,
                                      int direction,
                                      std::size_t cell) const
    {
        return buffer[dist_index(direction, cell)];
    }

    double& distribution(std::vector<double>& buffer, int direction, std::size_t cell)
    {
        return buffer[dist_index(direction, cell)];
    }

    std::vector<double> f;
    std::vector<double> f_next;
    std::vector<double> rho;
    std::vector<double> ux;
    std::vector<double> uy;
    std::vector<std::uint8_t> solid;
    std::vector<double> phi;
    std::vector<double> normal_x;
    std::vector<double> normal_y;

private:
    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
};

} // namespace lbm
