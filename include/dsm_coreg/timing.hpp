#pragma once

#include <chrono>
#include <map>
#include <mutex>
#include <string>

namespace dsm_coreg {

// Central timing registry used to report wall-clock time for the expensive
// stages of the workflow.  The registry stores seconds rather than formatted
// strings so the results JSON remains machine-readable.
class TimingRegistry {
public:
    void add_seconds(const std::string& name, double seconds) {
        std::lock_guard<std::mutex> lock(mutex_);
        seconds_[name] += seconds;
    }

    [[nodiscard]] std::map<std::string, double> snapshot() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return seconds_;
    }

private:
    mutable std::mutex mutex_;
    std::map<std::string, double> seconds_;
};

// RAII timer: construct at the start of a scope, and the elapsed wall-clock
// duration is automatically accumulated when the object leaves scope.
class ScopedTimer {
public:
    ScopedTimer(TimingRegistry& registry, std::string name)
        : registry_(registry), name_(std::move(name)), start_(Clock::now()) {}

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

    ~ScopedTimer() {
        const auto stop = Clock::now();
        const std::chrono::duration<double> dt = stop - start_;
        registry_.add_seconds(name_, dt.count());
    }

private:
    using Clock = std::chrono::steady_clock;

    TimingRegistry& registry_;
    std::string name_;
    Clock::time_point start_;
};

}  // namespace dsm_coreg
