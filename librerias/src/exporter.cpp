#include "exporter.hpp"

#include "lbm_utils.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace lbm {

namespace {

std::string obstacle_type_name(ObstacleType type)
{
    if (type == ObstacleType::Box) {
        return "box";
    }
    if (type == ObstacleType::Airfoil) {
        return "airfoil";
    }
    if (type == ObstacleType::ImportedSdf) {
        return "sdf";
    }
    return "circle";
}

std::string boundary_mode_name(BoundaryMode mode)
{
    return mode == BoundaryMode::Channel ? "channel" : "open_flow";
}

void write_raw(std::ofstream& out, const std::vector<float>& values)
{
    out.write(reinterpret_cast<const char*>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(float)));
    if (!out) {
        throw std::runtime_error("failed while writing binary frame.");
    }
}

void write_double_as_float(std::ofstream& out, const std::vector<double>& values)
{
    std::vector<float> tmp(values.size());
    std::transform(values.begin(), values.end(), tmp.begin(), [](double value) {
        return static_cast<float>(value);
    });
    write_raw(out, tmp);
}

void write_solid_as_float(std::ofstream& out, const std::vector<std::uint8_t>& values)
{
    std::vector<float> tmp(values.size());
    std::transform(values.begin(), values.end(), tmp.begin(), [](std::uint8_t value) {
        return value ? 1.0F : 0.0F;
    });
    write_raw(out, tmp);
}

} // namespace

Exporter::Exporter(SimulationParams params)
    : params_(std::move(params))
    , output_dir_(params_.output_dir)
    , frames_dir_(output_dir_ / "frames")
{
}

void Exporter::prepare(const Field& field)
{
    utils::ensure_directory(output_dir_);
    utils::ensure_directory(frames_dir_);

    for (const auto& entry : std::filesystem::directory_iterator(frames_dir_)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin") {
            std::filesystem::remove(entry.path());
        }
    }

    std::ofstream forces_file(output_dir_ / "forces.csv", std::ios::trunc);
    if (!forces_file) {
        throw std::runtime_error("could not open forces.csv for writing.");
    }
    forces_file << "step,time,Fx,Fy,Cd,Cl\n";

    std::ofstream diagnostics_file(output_dir_ / "diagnostics.csv", std::ios::trunc);
    if (!diagnostics_file) {
        throw std::runtime_error("could not open diagnostics.csv for writing.");
    }
    diagnostics_file
        << "step,time,total_mass,fluid_mass,mean_rho,min_rho,max_rho,"
        << "max_speed,kinetic_energy,solid_cells,fluid_cells,finite\n";

    frame_count_ = 0;
    write_metadata(field);
}

std::filesystem::path Exporter::frame_path(std::size_t step) const
{
    return frames_dir_ / ("frame_" + utils::format_step(step) + ".bin");
}

std::vector<std::string> Exporter::field_names() const
{
    return {"rho", "ux", "uy", "pressure", "vorticity", "solid", "phi"};
}

void Exporter::write_frame(const Field& field, std::size_t step)
{
    std::ofstream out(frame_path(step), std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("could not open frame file for writing.");
    }

    std::vector<float> pressure;
    std::vector<float> vorticity;
    lbm_utils::compute_pressure(field, pressure);
    lbm_utils::compute_vorticity(field, vorticity);

    write_double_as_float(out, field.rho);
    write_double_as_float(out, field.ux);
    write_double_as_float(out, field.uy);
    write_raw(out, pressure);
    write_raw(out, vorticity);
    write_solid_as_float(out, field.solid);
    write_double_as_float(out, field.phi);

    ++frame_count_;
    write_metadata(field);
}

void Exporter::append_force(const forces::ForceSample& sample)
{
    std::ofstream out(output_dir_ / "forces.csv", std::ios::app);
    if (!out) {
        throw std::runtime_error("could not open forces.csv for appending.");
    }

    out << std::setprecision(17)
        << sample.step << ','
        << sample.time << ','
        << sample.fx << ','
        << sample.fy << ','
        << sample.cd << ','
        << sample.cl << '\n';
}

void Exporter::append_diagnostics(const observables::Sample& sample)
{
    std::ofstream out(output_dir_ / "diagnostics.csv", std::ios::app);
    if (!out) {
        throw std::runtime_error("could not open diagnostics.csv for appending.");
    }

    out << std::setprecision(17)
        << sample.step << ','
        << sample.time << ','
        << sample.total_mass << ','
        << sample.fluid_mass << ','
        << sample.mean_rho << ','
        << sample.min_rho << ','
        << sample.max_rho << ','
        << sample.max_speed << ','
        << sample.kinetic_energy << ','
        << sample.solid_cells << ','
        << sample.fluid_cells << ','
        << (sample.finite ? 1 : 0) << '\n';
}

void Exporter::write_metadata(const Field& field) const
{
    std::ofstream out(output_dir_ / "metadata.json", std::ios::trunc);
    if (!out) {
        throw std::runtime_error("could not open metadata.json for writing.");
    }

    const auto names = field_names();
    out << std::setprecision(17);
    out << "{\n";
    out << "  \"case_name\": \"" << utils::escape_json(params_.case_name) << "\",\n";
    out << "  \"nx\": " << field.nx() << ",\n";
    out << "  \"ny\": " << field.ny() << ",\n";
    out << "  \"dx\": " << params_.dx << ",\n";
    out << "  \"dt\": " << params_.dt << ",\n";
    out << "  \"dtype\": \"float32\",\n";
    out << "  \"layout\": \"row-major\",\n";
    out << "  \"index_formula\": \"index = x + nx * y\",\n";
    out << "  \"fields\": [";
    for (std::size_t i = 0; i < names.size(); ++i) {
        out << (i == 0 ? "" : ", ") << '"' << names[i] << '"';
    }
    out << "],\n";
    out << "  \"frame_stride\": " << params_.frame_stride << ",\n";
    out << "  \"force_stride\": " << params_.force_stride << ",\n";
    out << "  \"diagnostics_stride\": " << params_.force_stride << ",\n";
    out << "  \"steps\": " << params_.steps << ",\n";
    out << "  \"frames_written\": " << frame_count_ << ",\n";
    out << "  \"rho0\": " << params_.rho0 << ",\n";
    out << "  \"u_in\": " << params_.u_in << ",\n";
    out << "  \"tau\": " << params_.tau << ",\n";
    out << "  \"omega\": " << params_.omega() << ",\n";
    out << "  \"viscosity\": " << params_.viscosity() << ",\n";
    out << "  \"reynolds\": " << params_.reynolds() << ",\n";
    out << "  \"start_from_rest\": " << (params_.start_from_rest ? "true" : "false") << ",\n";
    out << "  \"inlet_ramp_steps\": " << params_.inlet_ramp_steps << ",\n";
    out << "  \"boundary_mode\": \"" << boundary_mode_name(params_.boundary_mode) << "\",\n";
    out << "  \"obstacle\": {\n";
    out << "    \"type\": \"" << obstacle_type_name(params_.obstacle.type) << "\",\n";
    out << "    \"center_x\": " << params_.obstacle_center_x() << ",\n";
    out << "    \"center_y\": " << params_.obstacle_center_y() << ",\n";
    out << "    \"radius\": " << params_.obstacle.radius << ",\n";
    out << "    \"half_width\": " << params_.obstacle.half_width << ",\n";
    out << "    \"half_height\": " << params_.obstacle.half_height << ",\n";
    out << "    \"angle_degrees\": " << params_.obstacle.angle_degrees << ",\n";
    out << "    \"airfoil_chord\": " << params_.obstacle.airfoil_chord << ",\n";
    out << "    \"airfoil_max_camber\": " << params_.obstacle.airfoil_max_camber << ",\n";
    out << "    \"airfoil_camber_position\": " << params_.obstacle.airfoil_camber_position << ",\n";
    out << "    \"airfoil_max_thickness\": " << params_.obstacle.airfoil_max_thickness << ",\n";
    out << "    \"airfoil_samples\": " << params_.obstacle.airfoil_samples << ",\n";
    out << "    \"sdf_file\": \"" << utils::escape_json(params_.obstacle.sdf_file) << "\",\n";
    out << "    \"sdf_scale\": " << params_.obstacle.sdf_scale << "\n";
    out << "  }\n";
    out << "}\n";
}

} // namespace lbm
