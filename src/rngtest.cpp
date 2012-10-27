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
void run_c11_tests(const std::string& rng_name, Generator& rng, TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::cerr << "-- Starting C++11 " << rng_name << " RNG tests\n";
    RunClock * clock = nullptr;
    std::vector<double> params;

    params.reserve(nreps);
    clock = time_logger.get_timer("C++11", rng_name, "Unif[0,1]");
    std::uniform_real_distribution<> u01(0, 1);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        params.push_back(u01(rng));
    }
    clock->stop();
    std::cerr << rng_name << ": uniform real random variates [0,1) x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Exp[0.2]");
    std::exponential_distribution<> e1(0.2);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        e1(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": exponential with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Exp[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        std::exponential_distribution<> e2(params[rep]);
        e2(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": exponential with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Poisson[0.2]");
    std::poisson_distribution<> p1(0.2);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        p1(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Poisson[0.02]");
    std::poisson_distribution<> p2(0.02);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        p2(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Poisson[0.002]");
    std::poisson_distribution<> p3(0.002);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        p3(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Poisson[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        std::poisson_distribution<> p4(params[rep]);
        p4(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": poisson with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Geometric[0.2]");
    std::geometric_distribution<> g1(0.2);
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        g1(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": geometric with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

    clock = time_logger.get_timer("C++11", rng_name, "Geometric[p]");
    clock->start();
    for (unsigned int rep = 0; rep < nreps; ++rep) {
        std::geometric_distribution<> g2(params[rep]);
        g2(rng);
    }
    clock->stop();
    std::cerr << rng_name << ": geometric with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;
    std::cerr << std::endl;
}

void run_c11_mt19937_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::mt19937 rng(rd());
    run_c11_tests("MT19937", rng, time_logger, nreps);
}

void run_c11_mt19937_64_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    run_c11_tests("MT19937_64", rng, time_logger, nreps);
}

void run_c11_minstd_rand_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::minstd_rand rng(rd());
    run_c11_tests("MINSTD_RAND", rng, time_logger, nreps);
}

void run_c11_minstd_rand0_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::minstd_rand0 rng(rd());
    run_c11_tests("MINSTD_RAND0", rng, time_logger, nreps);
}

void run_c11_ranlux24_base_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::ranlux24_base rng(rd());
    run_c11_tests("RANLUX24_BASE", rng, time_logger, nreps);
}

void run_c11_ranlux24_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::ranlux24 rng(rd());
    run_c11_tests("RANLUX24", rng, time_logger, nreps);
}

void run_c11_ranlux48_base_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::ranlux48_base rng(rd());
    run_c11_tests("RANLUX48_BASE", rng, time_logger, nreps);
}

void run_c11_ranlux48_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::ranlux48 rng(rd());
    run_c11_tests("RANLUX48", rng, time_logger, nreps);
}

void run_c11_knuth_b_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::random_device rd;
    std::knuth_b rng(rd());
    run_c11_tests("KNUTHB", rng, time_logger, nreps);
}


void run_gsl_rng_tests(TimeLogger& time_logger, unsigned int nreps=DEFAULT_NREPS) {
    std::string implementation = "GSL";
    std::map<std::string, const gsl_rng_type *> rng_types;
    rng_types["taus"] = gsl_rng_taus2;
    rng_types["mt19937"] = gsl_rng_mt19937;
    rng_types["ranlux"] = gsl_rng_ranlux;
    // rng_types["gfsr4"] = gsl_rng_gfsr4;
    // rng_types["mrg"] = gsl_rng_mrg;
    // rng_types["cmrg"] = gsl_rng_cmrg;
    // rng_types["ranlxd1"] = gsl_rng_ranlxd1;
    // rng_types["ranlxs0"] = gsl_rng_ranlxs0;
    for (auto rng_type : rng_types) {
        const std::string& gen_alg = rng_type.first;
        RandomNumberGenerator rng(rng_type.second);
        RunClock * clock = nullptr;
        std::vector<double> params;

        std::cerr << "-- Starting " << implementation << " " << gen_alg << " RNG tests\n";

        params.reserve(nreps);
        clock = time_logger.get_timer(implementation, gen_alg, "Unif[0,1]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            params.push_back(rng.uniform_real());
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "uniform real random variates [0,1) x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Exp[0.2]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.exponential(0.2);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "exponential with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Exp[p]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.exponential(params[rep]);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "exponential with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Poisson[0.2]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.poisson(0.2);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Poisson[0.02]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.poisson(0.02);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Poisson[0.002]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.poisson(0.002);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "poisson with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Poisson[p]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.poisson(params[rep]);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "poisson with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Geometric[0.2]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.geometric(0.2);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "geometric with fixed parameter x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;

        clock = time_logger.get_timer(implementation, gen_alg, "Geometric[p]");
        clock->start();
        for (unsigned int rep = 0; rep < nreps; ++rep) {
            rng.geometric(params[rep]);
        }
        clock->stop();
        std::cerr << implementation << gen_alg << "geometric with varying parameters x " << nreps << ":\t" << clock->get_elapsed_seconds() << std::endl;
        std::cerr << std::endl;
    }
}

int main() {
    unsigned int    nreps = DEFAULT_NREPS;
    TimeLogger      time_logger;
    for (unsigned int i = 0; i < DEFAULT_REPEAT; ++i) {
        run_gsl_rng_tests(time_logger, nreps);
        run_c11_mt19937_tests(time_logger, nreps);
        run_c11_mt19937_64_tests(time_logger, nreps);
        run_c11_ranlux24_base_tests(time_logger, nreps);
        run_c11_ranlux24_tests(time_logger, nreps);
        run_c11_ranlux48_base_tests(time_logger, nreps);
        run_c11_ranlux48_tests(time_logger, nreps);
        run_c11_minstd_rand_tests(time_logger, nreps);
        run_c11_minstd_rand0_tests(time_logger, nreps);
        run_c11_knuth_b_tests(time_logger, nreps);
    }
    std::cerr << std::flush;
    std::ofstream flat("results.flat.txt");
    time_logger.summarize(flat);
    std::ofstream grouped("results.grouped.txt");
    time_logger.summarize_by_operation(grouped);
    std::ofstream best("results.best.txt");
    time_logger.summarize_best_by_operation(best);
}
