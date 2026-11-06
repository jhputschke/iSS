// Copyright @ Chun Shen 2018

#include "Random.h"

#include <gsl/gsl_cdf.h>

#include <cmath>
#include <math.h>   // lgamma_r (glibc)

namespace RandomUtil {

Random::Random(long seed_in, double min, double max) :
    rand_uniform_dist(min, max) {
    seed = seed_in;
    if (seed == -1) {
        seed = ran_dev();
    }
    ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));
}


namespace {

//! lgamma without the global signgam: glibc's lgamma() stores the sign there,
//! and threads calling it at once (correlated sampling) all write that one
//! cache line. lgamma_r is the same computation with the sign in a local.
inline double lgamma_no_signgam(const double x) {
#if defined(__GLIBC__)
    int sign;
    return ::lgamma_r(x, &sign);
#else
    return std::lgamma(x);
#endif
}

}  // namespace

long poisson_inverse(const double u, const double lambda) {
// The smallest k with P(k; lambda) >= u: monotonic in u and continuous in
// lambda, so the same u gives the same k for nearly equal lambdas.
    if (lambda <= 0.) return 0;
    const long mode = static_cast<long>(lambda);
    double p_k = std::exp(-lambda + mode*std::log(lambda)
                          - lgamma_no_signgam(mode + 1.));
    double cdf_k = (mode == 0) ? p_k : gsl_cdf_poisson_P(mode, lambda);
    long k = mode;
    if (u <= cdf_k) {
        // walk down: P(k-1) = P(k) - p(k)
        while (k > 0) {
            const double cdf_below = cdf_k - p_k;
            if (u > cdf_below) break;
            p_k *= k/lambda;
            cdf_k = cdf_below;
            k--;
        }
    } else {
        while (u > cdf_k) {
            k++;
            p_k *= lambda/k;
            cdf_k += p_k;
            if (p_k < 1e-300 && k > lambda) break;  // u rounds above the CDF
        }
    }
    return k;
}

}  // namespace RandomUtil
