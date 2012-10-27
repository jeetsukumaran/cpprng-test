#include <iostream>
#include <map>
#include <chrono>
#include <random>
#include <string>
#include "rng.hpp"

#define DEFAULT_NREPS   10000

class RunClock {

    typedef std::chrono::time_point<std::chrono::high_resolution_clock> TimePointType;

    public:
        RunClock() {}
        void start() {
            this->begin_ = std::chrono::system_clock::now();
        }
        void stop() {
            this->end_ = std::chrono::system_clock::now();
        }
        long get_elapsed_microseconds() {
            return std::chrono::duration_cast<std::chrono::microseconds>(this->end_-this->begin_).count();
        }
        double get_elapsed_seconds() {
            return static_cast<double>(this->get_elapsed_microseconds())/1e6;
        }

    private:
        TimePointType begin_;
        TimePointType end_;

}; // RunClock

class TimeLogger {

    public:
        RunClock * new_timer(const std::string& title) {
            this->logs_[title] = RunClock();
            return &this->logs_[title];
        }

    private:
        std::map<std::string, RunClock>     logs_;


}; // TimeLogger


void run_gsl_rng_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::cerr << "-- Starting GSL-based RNG tests\n";
    RunClock * clock = nullptr;
    RandomNumberGenerator rng;
    std::vector<double> params;

    std::cerr << "GSL: generating " << nreps << " uniform real random variates [0,1)\n";
    params.reserve(nreps);
    clock = time_logger.new_timer("GSL: Unif[0,1]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        params.push_back(rng.uniform_real());
    }
    clock->stop();
    std::cerr << "GSL: uniform real random variates [0,1) x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    std::cerr << "GSL: generating " << nreps << " exponential random variates with fixed parameter\n";
    clock = time_logger.new_timer("GSL: Exp[0.02]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.exponential(0.02);
    }
    clock->stop();
    std::cerr << "GSL: exponential with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    std::cerr << "GSL: generating " << nreps << " exponential random variates with varying parameters\n";
    clock = time_logger.new_timer("GS: Exp[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.exponential(params[rep]);
    }
    clock->stop();
    std::cerr << "GSL: exponential with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    std::cerr << "GSL: generating " << nreps << " poisson random variates with fixed parameter\n";
    clock = time_logger.new_timer("GSL: Exp[0.02]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.poisson(0.02);
    }
    clock->stop();
    std::cerr << "GSL: poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    std::cerr << "GSL: generating " << nreps << " poisson random variates with varying parameters\n";
    clock = time_logger.new_timer("GS: Exp[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.poisson(params[rep]);
    }
    clock->stop();
    std::cerr << "GSL: poisson with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;




}

int main() {
    unsigned int    nreps = DEFAULT_NREPS;
    TimeLogger      time_logger;
    run_gsl_rng_tests(time_logger, nreps);
}
