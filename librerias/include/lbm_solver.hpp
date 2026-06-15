#pragma once

#include "field.hpp"
#include "params.hpp"

#include <cstddef>

namespace lbm {

class LBMSolver {
public:
    LBMSolver(SimulationParams params, Field field);

    void initialize();
    void step();
    void compute_macroscopic();

    [[nodiscard]] const Field& field() const { return field_; }
    [[nodiscard]] Field& field() { return field_; }
    [[nodiscard]] const SimulationParams& params() const { return params_; }
    [[nodiscard]] std::size_t current_step() const { return current_step_; }

private:
    [[nodiscard]] double collide_value(int direction, std::size_t cell) const;
    void collide_and_stream();

    SimulationParams params_;
    Field field_;
    std::size_t current_step_ = 0;
};

} // namespace lbm
