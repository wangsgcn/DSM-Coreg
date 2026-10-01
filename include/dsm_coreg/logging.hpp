#pragma once

#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>

namespace dsm_coreg {

// Thread-safe file logger.  Bootstrap fits may run in parallel, so even though
// the current implementation keeps per-replicate logging intentionally minimal,
// making the logger thread-safe avoids surprises as diagnostics evolve.
class Logger {
public:
    explicit Logger(const std::filesystem::path& path);

    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);
    void debug(const std::string& message);

private:
    void write(const char* level, const std::string& message);

    std::ofstream stream_;
    std::mutex mutex_;
};

}  // namespace dsm_coreg
