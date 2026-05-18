#ifndef PGPHASE_FISHER_EXACT_HPP
#define PGPHASE_FISHER_EXACT_HPP

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace pgphase_collect {

inline double log_hypergeom(int a, int b, int c, int d) {
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N = n1 + n2;
    if (N <= 0) return -std::numeric_limits<double>::infinity();
    if (n1 > n2) return log_hypergeom(c, d, a, b);
    if (m1 > m2) return log_hypergeom(b, a, d, c);
    return std::lgamma(static_cast<double>(n1 + 1)) + std::lgamma(static_cast<double>(n2 + 1)) +
           std::lgamma(static_cast<double>(m1 + 1)) + std::lgamma(static_cast<double>(m2 + 1)) -
           (std::lgamma(static_cast<double>(a + 1)) + std::lgamma(static_cast<double>(b + 1)) +
            std::lgamma(static_cast<double>(c + 1)) + std::lgamma(static_cast<double>(d + 1)) +
            std::lgamma(static_cast<double>(N + 1)));
}

/**
 * Two-sided Fisher exact test on a 2x2 contingency table.
 * @param a Alt-forward count.
 * @param b Alt-reverse count.
 * @param c,d Expected counts (symmetric layout).
 * @return Two-sided p-value in [0,1].
 */
inline double fisher_exact_two_tail(int a, int b, int c, int d) {
    if (a + b + c + d <= 0) return 1.0;
    const double log_p_observed = log_hypergeom(a, b, c, d);
    const double p_observed = std::exp(log_p_observed);
    double total_p = 0.0;
    int min_a = (0 > (a + c) - (a + b + c + d)) ? 0 : (a + c) - (b + d);
    const int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    const int denom = a + b + c + d;
    const int mode_a =
        denom > 0
            ? static_cast<int>(static_cast<double>(a + b) * static_cast<double>(a + c) /
                               static_cast<double>(denom))
            : 0;
    for (int delta = 0; delta <= max_a - min_a; ++delta) {
        int current_a = mode_a + delta;
        if (current_a <= max_a) {
            int current_b = (a + b) - current_a;
            int current_c = (a + c) - current_a;
            int current_d = (b + d) - current_b;
            if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                const double p = std::exp(log_hypergeom(current_a, current_b, current_c, current_d));
                if (p <= p_observed + DBL_EPSILON) total_p += p;
            }
        }
        if (delta > 0) {
            current_a = mode_a - delta;
            if (current_a >= min_a) {
                int current_b = (a + b) - current_a;
                int current_c = (a + c) - current_a;
                int current_d = (b + d) - current_b;
                if (current_b >= 0 && current_c >= 0 && current_d >= 0) {
                    const double p =
                        std::exp(log_hypergeom(current_a, current_b, current_c, current_d));
                    if (p <= p_observed + DBL_EPSILON) total_p += p;
                }
            }
        }
    }
    return total_p;
}

} // namespace pgphase_collect

#endif
