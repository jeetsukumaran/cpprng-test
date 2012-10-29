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
            typename decltype(this->uniform_real_rng_)::param_type param(0.0, 1.0);
            return this->uniform_real_rng_(this->engine_, param);
        }

        // // returns real value uniformly distributed in [a, b]
        inline double uniform_real(double a, double b) {
            typename decltype(this->uniform_real_rng_)::param_type param(a, b);
            return this->uniform_real_rng_(this->engine_, param);
        }

        // // returns an integer in {0, 1} sampled from a Bernoulli distribution with probability of 1
        // // given by p
        // inline unsigned int bernoulli(double p=0.5) {
        //     return gsl_ran_bernoulli(this->rgen_, p);
        // }

        // // returns a binomial random variate with probability of success in n trials given by p
        // inline unsigned long binomial(double p, unsigned long n) {
        //     return gsl_ran_binomial(this->rgen_, p, static_cast<unsigned int>(n));
        // }

        // returns a random variate from an exponential distribution with rate parameter p
        inline double exponential(double p) {
            if (p == 0) {
                return 0.0;
            }
            typename decltype(this->exponential_rng_)::param_type param(p);
            return this->exponential_rng_(this->engine_, param);
        }

        // // returns a random variate from normal distribution with mean 0 and standard deviation sd
        // inline double gaussian(double sd) {
        //     return gsl_ran_gaussian_ziggurat(this->rgen_, sd);
        // }

        // returns an integer sampled from a Geometric distribution with rate parameter p,
        // with support in {0, 1, 2, ...}
        inline unsigned long geometric(double p) {
            typename decltype(this->geometric_rng_)::param_type param(p);
            return this->geometric_rng_(this->engine_, param);
        }

        // // returns a set of numbers sampled from the multinomial distribution, where:
        // //  - K is the number of categories
        // //  - N is the total number of elements in each category
        // //  - p[] is a K-length array of weights (need not sum to 1)
        // //  - n[] is a K-length array that will be populated with the sampled numbers
        // unsigned int * multinomial(unsigned long K, unsigned int N, const double * p, unsigned int n[]) {
        //     gsl_ran_multinomial(this->rgen_, K, N, p, n);
        //     return n;
        // }

        // returns an integer value sampled from a Poisson distribution with mean mu
        inline unsigned long poisson(double mu) {
            if (mu <= 0) {
                return 0;
            }
            typename decltype(this->poisson_rng_)::param_type param(mu);
            return this->poisson_rng_(this->engine_, param);
        }

    private:
        typename Engine::result_type                seed_;
        Engine                                      engine_;
        std::uniform_real_distribution<double>      uniform_real_rng_;
        std::uniform_int_distribution<long>         uniform_int_rng_;
        std::exponential_distribution<double>       exponential_rng_;
        std::poisson_distribution<unsigned long>    poisson_rng_;
        std::geometric_distribution<unsigned long>  geometric_rng_;

    private:
        RandomNumberGenerator11(const RandomNumberGenerator11 &);
        RandomNumberGenerator11 & operator=(const RandomNumberGenerator11 &);

}; // RandomNumberGenerator11


#endif

