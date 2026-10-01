#include "dsm_coreg/logging.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dsm_coreg {
namespace {

std::string timestamp_now() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);

    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &tt);
#else
    localtime_r(&tt, &tm);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

}  // namespace

Logger::Logger(const std::filesystem::path& path)
    : stream_(path, std::ios::out | std::ios::trunc) {
    if (!stream_) {
        throw std::runtime_error("Failed to open log file: " + path.string());
    }
}

void Logger::info(const std::string& message) { write("INFO", message); }
void Logger::warning(const std::string& message) { write("WARN", message); }
void Logger::error(const std::string& message) { write("ERROR", message); }
void Logger::debug(const std::string& message) { write("DEBUG", message); }

void Logger::write(const char* level, const std::string& message) {
    std::lock_guard<std::mutex> lock(mutex_);
    stream_ << timestamp_now() << " [" << level << "] " << message << '\n';
    stream_.flush();
}

}  // namespace dsm_coreg
