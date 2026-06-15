#include "utils.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace lbm::utils {

double deg_to_rad(double degrees)
{
    return degrees * std::acos(-1.0) / 180.0;
}

std::string format_step(std::size_t step, int width)
{
    std::ostringstream out;
    out << std::setw(width) << std::setfill('0') << step;
    return out.str();
}

std::string escape_json(const std::string& value)
{
    std::ostringstream out;
    for (const char ch : value) {
        switch (ch) {
        case '\\':
            out << "\\\\";
            break;
        case '"':
            out << "\\\"";
            break;
        case '\n':
            out << "\\n";
            break;
        case '\r':
            out << "\\r";
            break;
        case '\t':
            out << "\\t";
            break;
        default:
            out << ch;
            break;
        }
    }
    return out.str();
}

void ensure_directory(const std::filesystem::path& path)
{
    std::error_code ec;
    std::filesystem::create_directories(path, ec);
    if (ec) {
        throw std::runtime_error("could not create directory '" + path.string() + "': " + ec.message());
    }
}

Timer::Timer()
{
    reset();
}

void Timer::reset()
{
    start_ = std::chrono::steady_clock::now();
}

double Timer::elapsed_seconds() const
{
    const auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_).count();
}

} // namespace lbm::utils
