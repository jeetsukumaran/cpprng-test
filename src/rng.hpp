#ifndef RNG_HPP_
#define RNG_HPP_

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

// inline calls to gsl
#define HAVE_INLINE 1

typedef unsigned long RandomSeedType;

class RandomNumberGenerator {

    public:

        RandomNumberGenerator();
        RandomNumberGenerator(const gsl_rng_type * T);
        RandomNumberGenerator(RandomSeedType rng_seed);
        ~RandomNumberGenerator();

        void init(const gsl_rng_type * T = nullptr);
        RandomSeedType get_seed() const;
        void set_seed(RandomSeedType seed);
        void set_seed_from_time();
        inline unsigned long max_int() const {
            return this->max_int_;
        }
        inline unsigned long min_int() const {
            return this->min_int_;
        }

        // returns next number from random sequence
        unsigned long get() {
            return gsl_rng_get(this->rgen_);
        }

        // returns real value uniformly distributed in [0, 1)
        inline double uniform_real() {
            return gsl_rng_uniform(this->rgen_);
        }

        // returns real value uniformly distributed in [a, b]
        inline double uniform_real(double a, double b) {
            return gsl_ran_flat(this->rgen_, a, b);
        }

        // returns real value uniformly distributed in (0, 1)
        inline double uniform_real_pos() {
            return gsl_rng_uniform(this->rgen_);
        }

        // returns an integer value in [0, n), i.e. 0 to n-1 inclusive
        inline unsigned long uniform_0n(unsigned long n) {
            return gsl_rng_uniform_int(this->rgen_, n);
        }

        // returns integer value uniformly distributed in [a, b]
        inline long uniform_int(long a, long b) {
            return (static_cast<long>(gsl_rng_uniform_int(this->rgen_, b-a+1))) + a;
        }

        // returns an integer in {0, 1} sampled from a Bernoulli distribution with probability of 1
        // given by p
        inline unsigned int bernoulli(double p=0.5) {
            return gsl_ran_bernoulli(this->rgen_, p);
        }

        // returns a binomial random variate with probability of success in n trials given by p
        inline unsigned long binomial(double p, unsigned long n) {
            return gsl_ran_binomial(this->rgen_, p, static_cast<unsigned int>(n));
        }

        // returns a random variate from an exponential distribution with rate parameter p
        inline double exponential(double p) {
            if (p == 0) {
                return 0;
            }
            return gsl_ran_exponential(this->rgen_, 1.0/p);
        }

        // returns a random variate from normal distribution with mean 0 and standard deviation sd
        inline double gaussian(double sd) {
            return gsl_ran_gaussian_ziggurat(this->rgen_, sd);
        }

        // returns an integer sampled from a Geometric distribution with rate parameter p,
        // with support in {1, 2, ...}
        inline unsigned long geometric(double p) {
            return gsl_ran_geometric(this->rgen_, p);
        }

        // returns a set of numbers sampled from the multinomial distribution, where:
        //  - K is the number of categories
        //  - N is the total number of elements in each category
        //  - p[] is a K-length array of weights (need not sum to 1)
        //  - n[] is a K-length array that will be populated with the sampled numbers
        unsigned int * multinomial(unsigned long K, unsigned int N, const double * p, unsigned int n[]) {
            gsl_ran_multinomial(this->rgen_, K, N, p, n);
            return n;
        }

        // returns a real value sampled from a Poisson distribution with rate parameter p
        inline double poisson(double p) {
            if (p == 0) {
                return 0;
            }
            return gsl_ran_poisson(this->rgen_, 1.0/p);
        }

    private:
        RandomSeedType          seed_;
        gsl_rng *               rgen_;
        unsigned long           max_int_;
        unsigned long           min_int_;

    private:
        RandomNumberGenerator(const RandomNumberGenerator &);
        RandomNumberGenerator & operator=(const RandomNumberGenerator &);

}; // RandomNumberGenerator


#endif

