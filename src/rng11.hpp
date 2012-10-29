#ifndef RNG11_HPP_
#define RNG11_HPP_

#include <random>

template <typename Engine>
class RandomNumberGenerator11 {

    public:

        RandomNumberGenerator11() {
            this->set_seed_from_time();
        }
        RandomNumberGenerator11(RandomSeedType rng_seed) {
            this->set_seed_(rng_seed);
        }
        ~RandomNumberGenerator11() {};

        RandomSeedType get_seed() const {
            return this->engine_;
        }

        void set_seed(typename Engine::result_type seed) {
            this->seed_ = seed;
            this->engine_.seed(seed);
        }

        void set_seed_from_time() {
            this->set_seed(static_cast<typename Engine::result_type>(std::time(NULL)));
        }

        void set_seed_from_from_device() {
            std::random_device rd;
            this->set_seed(rd());
        }

        // returns real value uniformly distributed in [0, 1]
        inline double uniform_real() {
            return this->uniform_real_rng_(this->engine_,
                    typename decltype(this->uniform_real_rng_)::param_type (0.0, 1.0));
        }

        // // returns real value uniformly distributed in [a, b]
        inline double uniform_real(double a, double b) {
            return this->uniform_real_rng_(this->engine_,
                    typename decltype(this->uniform_real_rng_)::param_type (a, b));
        }

        // returns an integer in {0, 1} sampled from a Bernoulli distribution with probability of 1
        // given by p
        inline bool bernoulli(double p=0.5) {
            return this->bernoulli_rng_(this->engine_,
                    typename decltype(this->bernoulli_rng_)::param_type (p));
        }

        // returns a binomial random variate with probability of success in n trials given by p
        inline unsigned long binomial(double p, unsigned long n) {
            return this->binomial_rng_(this->engine_,
                    typename decltype(this->binomial_rng_)::param_type (p, n));
        }

        // returns a random variate from an exponential distribution with rate parameter p
        inline double exponential(double p) {
            if (p == 0) {
                return 0.0;
            }
            return this->exponential_rng_(this->engine_,
                    typename decltype(this->exponential_rng_)::param_type (p));
        }

        // returns a random variate from normal distribution with given mean and standard deviation
        inline double normal(double mean, double stddev) {
            return this->normal_rng_(this->engine_,
                    typename decltype(this->normal_rng_)::param_type (mean, stddev));
        }

        // returns an integer sampled from a Geometric distribution with rate parameter p,
        // with support in {0, 1, 2, ...}
        inline unsigned long geometric(double p) {
            return this->geometric_rng_(this->engine_,
                    typename decltype(this->geometric_rng_)::param_type (p));
        }

        // returns an integer value sampled from a Poisson distribution with mean mu
        inline unsigned long poisson(double mu) {
            if (mu <= 0) {
                return 0;
            }
            return this->poisson_rng_(this->engine_,
                    typename decltype(this->poisson_rng_)::param_type (mu));
        }

    private:
        typename Engine::result_type                seed_;
        Engine                                      engine_;
        std::uniform_real_distribution<double>      uniform_real_rng_;
        std::uniform_int_distribution<long>         uniform_int_rng_;
        std::bernoulli_distribution                 bernoulli_rng_;
        std::binomial_distribution<unsigned long>   binomial_rng_;
        std::exponential_distribution<double>       exponential_rng_;
        std::normal_distribution<double>            normal_rng_;
        std::poisson_distribution<unsigned long>    poisson_rng_;
        std::geometric_distribution<unsigned long>  geometric_rng_;

    private:
        RandomNumberGenerator11(const RandomNumberGenerator11 &);
        RandomNumberGenerator11 & operator=(const RandomNumberGenerator11 &);

}; // RandomNumberGenerator11


#endif

