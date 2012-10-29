#include <utility>
#include <algorithm>
#include <exception>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <gsl/gsl_rng.h>
#include "rng.hpp"
#include "rng11.hpp"

#define DEFAULT_REPEAT  100
#define DEFAULT_NREPS   100000

class RunClock {

    typedef std::chrono::time_point<std::chrono::high_resolution_clock> TimePointType;

    public:
        RunClock(const std::string& implementation,
            const std::string& algorithm,
            const std::string& operation) :
            implementation_(implementation),
            algorithm_(algorithm),
            operation_(operation),
            elapsed_seconds_(-1) {
        }
        void start() {
            this->elapsed_seconds_ = -1;
            this->begin_ = std::move(std::chrono::system_clock::now());
        }
        void stop() {
            TimePointType end = std::chrono::system_clock::now();
            unsigned long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - this->begin_).count();
            this->elapsed_microseconds_.push_back(ms);
        }
        double get_elapsed_seconds() {
            if (this->elapsed_seconds_ < 0) {
                unsigned long total_microseconds = 0;
                for (auto ms : this->elapsed_microseconds_) {
                    total_microseconds += ms;
                }
                this->elapsed_seconds_ = static_cast<double>(total_microseconds) / (static_cast<double>(this->elapsed_microseconds_.size()) * 1e6);
            }
            return this->elapsed_seconds_;
        }
        void print(std::ostream& out) {
            out << std::setw(6) << std::left << this->implementation_ << " ";
            out << std::setw(25) << std::left << this->algorithm_ << " ";
            out << std::setw(25) << std::left << this->operation_ << " ";
            out << this->get_elapsed_seconds() << std::endl;
        }

    private:
        std::string                 implementation_;
        std::string                 algorithm_;
        std::string                 operation_;
        TimePointType               begin_;
        std::vector<unsigned long>         elapsed_microseconds_;
        double                      elapsed_seconds_;

}; // RunClock

bool cmp_results(RunClock * a, RunClock * b) {
    if (a->get_elapsed_seconds() < b->get_elapsed_seconds()) {
        return true;
    } else {
        return false;
    }
}

class TimeLogger {

    public:
        ~TimeLogger() {
            for (auto log : this->logs_) {
                delete log;
            }
        }
        RunClock * get_timer(const std::string& implementation, const std::string& algorithm, const std::string& operation) {
            RunClock * rc = nullptr;
            std::string rc_hash = implementation + "::" + algorithm + "::" + operation;
            std::map<std::string, RunClock *>::iterator rci = this->run_clocks_.find(rc_hash);
            if (rci != this->run_clocks_.end()) {
                rc = rci->second;
            } else {
                rc = new RunClock(implementation, algorithm, operation);
                this->run_clocks_[rc_hash] = rc;
                this->logs_.push_back(rc);
                this->logs_by_operation_[operation].push_back(rc);
            }
            // rc = new RunClock(implementation, algorithm, operation);
            // this->run_clocks_[rc_hash] = rc;
            // this->logs_.push_back(rc);
            // this->logs_by_operation_[operation].push_back(rc);
            return rc;
        }
        void summarize(std::ostream& out) {
            std::sort(this->logs_.begin(), this->logs_.end(), &cmp_results);
            for (auto rc : this->logs_) {
                rc->print(out);
            }
        }
        void summarize_by_operation(std::ostream& out) {
            for (auto log_by_operation : this->logs_by_operation_) {
                std::string operation = log_by_operation.first;
                out << "\n\n### " << operation << " ###\n\n";
                std::vector<RunClock *> logs = log_by_operation.second;
                std::sort(logs.begin(), logs.end(), &cmp_results);
                for (auto rc : logs) {
                    rc->print(out);
                }
            }
        }
        void summarize_best_by_operation(std::ostream& out) {
            for (auto log_by_operation : this->logs_by_operation_) {
                std::string operation = log_by_operation.first;
                std::vector<RunClock *> logs = log_by_operation.second;
                std::sort(logs.begin(), logs.end(), &cmp_results);
                (*logs.begin())->print(out);
            }
        }

    private:
        std::vector<RunClock *>                        logs_;
        std::map<std::string, std::vector<RunClock *>> logs_by_operation_;
        std::map<std::string, RunClock *>              run_clocks_;

}; // TimeLogger

template <class Generator>
void run_rng_tests(
        Generator& rng,
        const std::string& implementation,
        const std::string& gen_alg,
        TimeLogger& time_logger,
        unsigned int nreps=DEFAULT_NREPS) {

    RunClock * clock = nullptr;
    std::vector<double> params;
    std::vector<double> params2;

    std::cerr << "-- Starting " << implementation << " " << gen_alg << " RNG tests\n";

    params.reserve(nreps);
    params2.reserve(nreps);
    clock = time_logger.get_timer(implementation, gen_alg, "Unif[0,1]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        params.push_back(rng.uniform_real());
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " uniform real random variates [0,1) x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    for (auto p : params) {
        params2.push_back(1.0/p);
    }

    clock = time_logger.get_timer(implementation, gen_alg, "Exp[0.2]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.exponential(0.2);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " exponential with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Exp[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.exponential(params[rep]);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " exponential with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Poisson[2.0]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.poisson(2.0);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Poisson[20.0]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.poisson(20.0);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Poisson[200]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.poisson(200);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Poisson[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.poisson(params2[rep]);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " poisson with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Geometric[0.2]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.geometric(0.2);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " geometric with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer(implementation, gen_alg, "Geometric[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        rng.geometric(params[rep]);
    }
    clock->stop();
    std::cerr << implementation << " " << gen_alg << " geometric with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;
    std::cerr << std::endl;
}


void run_gsl_rng_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::string implementation = "GSL";
    std::map<std::string, const gsl_rng_type *> rng_types;
    rng_types["taus"] = gsl_rng_taus2;
    rng_types["mt19937"] = gsl_rng_mt19937;
    rng_types["ranlux"] = gsl_rng_ranlux;
    rng_types["gfsr4"] = gsl_rng_gfsr4;
    rng_types["mrg"] = gsl_rng_mrg;
    rng_types["cmrg"] = gsl_rng_cmrg;
    rng_types["ranlxd1"] = gsl_rng_ranlxd1;
    rng_types["ranlxs0"] = gsl_rng_ranlxs0;
    for (auto rng_type : rng_types) {
        const std::string& gen_alg = rng_type.first;
        RandomNumberGenerator rng(rng_type.second);
        run_rng_tests(rng, implementation, gen_alg, time_logger, nreps);
    }
}

// int main() {
//     std::vector<double>     probs{0.01, 0.1, 1.0, 10, 100, 1000, 1e6, -0.01, -0.1, -1.0, -1000.0};
//     std::cout << "\n-- GSL\n" << std::endl;
//     RandomNumberGenerator gsl;
//     for (auto p : probs) {
//         auto r = gsl.poisson(p);
//         std::cout << std::setw(20) << std::left << p << r << std::endl;
//     }
//     std::cout << "\n\n--C++11\n" << std::endl;
//     std::random_device rd;
//     std::mt19937 rng(rd());
//     for (auto p : probs) {
//         std::poisson_distribution<> pois(p);
//         auto r = pois(rng);
//         std::cout << std::setw(20) << std::left << p << r << std::endl;
//     }
// }

void run_cpp11_tests( TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::string implementation = "C++11";
    RandomNumberGenerator11<std::minstd_rand0>      rng_minstd_rand0;
    RandomNumberGenerator11<std::minstd_rand>       rng_minstd_rand;
    RandomNumberGenerator11<std::mt19937>           rng_mt19937;
    RandomNumberGenerator11<std::mt19937_64>        rng_mt19937_64;
    RandomNumberGenerator11<std::ranlux24_base>     rng_ranlux24_base;
    RandomNumberGenerator11<std::ranlux48_base>     rng_ranlux48_base;
    RandomNumberGenerator11<std::ranlux24>          rng_ranlux24;
    RandomNumberGenerator11<std::ranlux48>          rng_ranlux48;
    RandomNumberGenerator11<std::knuth_b>           rng_knuth_b;
    run_rng_tests(rng_minstd_rand0, implementation, "minstd_rand0", time_logger, nreps);
    run_rng_tests(rng_minstd_rand, implementation, "minstd_rand", time_logger, nreps);
    run_rng_tests(rng_mt19937, implementation, "mt19937", time_logger, nreps);
    run_rng_tests(rng_mt19937_64, implementation, "mt19937_64", time_logger, nreps);
    run_rng_tests(rng_ranlux24_base, implementation, "ranlux24_base", time_logger, nreps);
    run_rng_tests(rng_ranlux48_base, implementation, "ranlux48_base", time_logger, nreps);
    run_rng_tests(rng_ranlux24, implementation, "ranlux24", time_logger, nreps);
    run_rng_tests(rng_ranlux48, implementation, "ranlux48", time_logger, nreps);
    run_rng_tests(rng_knuth_b, implementation, "knuth_b", time_logger, nreps);
}

int main() {
    unsigned int    nreps = DEFAULT_NREPS;
    TimeLogger      time_logger;
    for (unsigned int i = 0; i < DEFAULT_REPEAT; ++i) {
        run_gsl_rng_tests(time_logger, nreps);
        run_cpp11_tests(time_logger, nreps);
    }
    std::cerr << std::flush;
    std::ofstream flat("results.flat.txt");
    time_logger.summarize(flat);
    std::ofstream grouped("results.grouped.txt");
    time_logger.summarize_by_operation(grouped);
    std::ofstream best("results.best.txt");
    time_logger.summarize_best_by_operation(best);
}
