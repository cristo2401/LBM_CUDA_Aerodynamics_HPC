#pragma once

#include "field.hpp"
#include "forces.hpp"
#include "observables.hpp"
#include "params.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace lbm {

class Exporter {
public:
    explicit Exporter(SimulationParams params);

    void prepare(const Field& field);
    void write_frame(const Field& field, std::size_t step);
    void append_force(const forces::ForceSample& sample);
    void append_diagnostics(const observables::Sample& sample);
    void write_metadata(const Field& field) const;

    [[nodiscard]] std::size_t frame_count() const { return frame_count_; }

private:
    [[nodiscard]] std::filesystem::path frame_path(std::size_t step) const;
    [[nodiscard]] std::vector<std::string> field_names() const;

    SimulationParams params_;
    std::filesystem::path output_dir_;
    std::filesystem::path frames_dir_;
    std::size_t frame_count_ = 0;
};

} // namespace lbm
