#pragma once

#include <chrono>
#include <filesystem>
#include <string>

namespace lbm::utils {

template <class T>
constexpr T clamp(T value, T low, T high)
{
    return value < low ? low : (value > high ? high : value);
}

template <class T>
constexpr T lerp(T a, T b, T t)
{
    return a + (b - a) * t;
}

[[nodiscard]] double deg_to_rad(double degrees);
[[nodiscard]] std::string format_step(std::size_t step, int width = 6);
[[nodiscard]] std::string escape_json(const std::string& value);

void ensure_directory(const std::filesystem::path& path);

class Timer {
public:
    Timer();
    void reset();
    [[nodiscard]] double elapsed_seconds() const;

private:
    std::chrono::steady_clock::time_point start_;
};

} // namespace lbm::utils
