#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <gsl/gsl_rng.h>
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

bool cmp_results(std::pair<std::string, double> a, std::pair<std::string, double> b) {
    if (a.second < b.second) {
        return true;
    } else {
        return false;
    }
}

class TimeLogger {

    public:
        RunClock * new_timer(const std::string& implementation, const std::string& subtype, const std::string& operation) {
            std::ostringstream o;
            o << std::setw(6) << std::left << implementation << "  " << std::setw(20) << std::left << subtype << "  " << std::setw(15) << std::left << operation;
            std::string title = o.str();
            this->logs_[title] = RunClock();
            return &this->logs_[title];
        }
        void summarize(std::ostream& out) {
            std::vector<std::pair<std::string, double> > results;
            results.reserve(this->logs_.size());
            for (auto log : this->logs_) {
                const std::string& title = log.first;
                double seconds = log.second.get_elapsed_seconds();
                results.push_back(std::pair<std::string, double>(title, seconds));
            }
            std::sort(results.begin(), results.end(), &cmp_results);
            for (auto result : results) {
                out << std::setw(50) << std::left << result.first << "  " << result.second << std::endl;
            }
        }

    private:
        std::map<std::string, RunClock>     logs_;

}; // TimeLogger

template <class Generator>
void run_c11_tests(const std::string& rng_name, Generator& rng, TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::cerr << "-- Starting " << rng_name << "-based RNG tests\n";
    RunClock * clock = nullptr;
    std::vector<double> params;

    params.reserve(nreps);
    clock = time_logger.new_timer("C++11", rng_name, "Unif[0,1]");
    std::uniform_real_distribution<> u01(0, 1);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        params.push_back(u01(rng));
    }
    clock->stop();
    std::cerr << rng_name << ": uniform real random variates [0,1) x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.new_timer("C++11", rng_name, "Exp[0.02]");
    std::exponential_distribution<> e1(0.02);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        e1(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": exponential with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.new_timer("C++11", rng_name, "Exp[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        std::exponential_distribution<> e2(params[rep]);
        e2(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": exponential with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.new_timer("C++11", rng_name, "Poisson[0.02]");
    std::poisson_distribution<> p1(0.02);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        p1(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.new_timer("C++11", rng_name, "Poisson[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        std::poisson_distribution<> p2(params[rep]);
        p2(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": poisson with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;
    std::cerr << std::endl;
}

void run_c11_mt19937_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::mt19937 rng(rd());
    run_c11_tests("MT19937", rng, time_logger, nreps);
}

void run_c11_ranlux24_base_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::ranlux24_base rng(rd());
    run_c11_tests("RANLUX24_BASE", rng, time_logger, nreps);
}


void run_gsl_rng_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::string implementation = "GSL";
    std::map<std::string, const gsl_rng_type *> rng_types;
    rng_types["taus"] = gsl_rng_taus;
    rng_types["mt19937"] = gsl_rng_mt19937;
    for (auto rng_type : rng_types) {
        const std::string& gen_alg = rng_type.first;
        RandomNumberGenerator rng(rng_type.second);
        RunClock * clock = nullptr;
        std::vector<double> params;

        std::cerr << "-- Starting " << implementation << " " << gen_alg << " RNG tests\n";

        params.reserve(nreps);
        clock = time_logger.new_timer(implementation, gen_alg, "Unif[0,1]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            params.push_back(rng.uniform_real());
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "uniform real random variates [0,1) x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.new_timer(implementation, gen_alg, "Exp[0.02]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.exponential(0.02);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "exponential with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.new_timer(implementation, gen_alg, "Exp[p]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.exponential(params[rep]);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "exponential with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.new_timer(implementation, gen_alg, "Poisson[0.02]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.poisson(0.02);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.new_timer(implementation, gen_alg, "Poisson[p]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.poisson(params[rep]);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "poisson with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;
        std::cerr << std::endl;
    }
}

int main() {
    unsigned int    nreps = DEFAULT_NREPS;
    TimeLogger      time_logger;
    run_gsl_rng_tests(time_logger, nreps);
    run_c11_mt19937_tests(time_logger, nreps);
    run_c11_ranlux24_base_tests(time_logger, nreps);
    std::cerr << "\n\n---\nResults:\n---\n\n";
    std::cerr << std::flush;
    time_logger.summarize(std::cout);
}
