#include <ctime>
#include "rng.hpp"

RandomNumberGenerator::RandomNumberGenerator() {
    this->init();
    this->set_seed_from_time();
}

RandomNumberGenerator::RandomNumberGenerator(RandomSeedType rng_seed) {
    this->init();
    this->set_seed(rng_seed);
}

RandomNumberGenerator::~RandomNumberGenerator() {
    gsl_rng_free(this->rgen_);
}

void RandomNumberGenerator::init() {
    this->rgen_ = gsl_rng_alloc(gsl_rng_taus);
    this->max_int_ = gsl_rng_max(this->rgen_);
    this->min_int_ = gsl_rng_min(this->rgen_);
}

RandomSeedType RandomNumberGenerator::get_seed() const {
    return this->seed_;
}

void RandomNumberGenerator::set_seed(RandomSeedType seed) {
    this->seed_ = seed;
    gsl_rng_set(this->rgen_, static_cast<unsigned long int>(this->seed_));
}

void RandomNumberGenerator::set_seed_from_time() {
    this->set_seed(std::time(NULL));
}
