#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "minpack.h"

#define check(x, ...) \
    _Generic((x), \
        int: check_int, \
     double: check_double \
            )(x, __VA_ARGS__)

static inline bool
check_int(int expected, int actual, const char *msg)
{
    if (expected == actual) {
        return true;
    }
    fprintf(stderr, "FAIL: %s: expected %d, got %d\n", msg, expected, actual);
    return false;
}

static inline bool
check_double(double expected, double actual, double tol, const char *msg)
{
    if (fabs(expected - actual) < tol) {
        return true;
    }
    fprintf(stderr, "FAIL: %s: expected %g, got %g\n", msg, expected, actual);
    return false;
}

static inline bool
is_close(double a, double b, double tol)
{
    return fabs(a - b) < tol;
}

static double
enorm(const int n, const double* x)
{
    double norm = 0.0;
    for(int i = 0; i < n; i++) {
        norm += x[i]*x[i];
    }
    return sqrt(norm);
}

void
trial_hybrd_fcn(int n, const double* x, double* fvec, int* iflag, void* udata) {
    assert(!udata);  // we always pass NULL in the test
    if (*iflag == 0) return;
    for(int k = 0; k < n; k++) {
        double temp = (3.0 - 2.0*x[k])*x[k];
        double temp1 = k != 0 ? x[k - 1] : 0.0;
        double temp2 = k != n-1 ? x[k + 1] : 0.0;
        fvec[k] = temp - temp1 - 2.0*temp2 + 1.0;
    }
}

static int
test_hybrd1 (void)
{
    int n = 9;
    int info = 0;
    int lwa = 180;
    double x[9] = {-1.0};
    double fvec[9] = {0.0};
    double wa[180] = {0.0};
    double tol = sqrt(minpack_dpmpar(1));
    double reference[9] = {
        -0.5706545, -0.6816283, -0.7017325,
        -0.7042129, -0.7013690, -0.6918656,
        -0.6657920, -0.5960342, -0.4164121};

    minpack_hybrd1(trial_hybrd_fcn, n, x, fvec, tol, &info, wa, lwa, NULL);

    if (!check(info, 1, "Unexpected info value")) return 1;

    if (!check(enorm(n, fvec), 0.0, tol, "Unexpected residual")) return 1;

    for(int i = 0; i < 9; i++) {
        if (!check(x[i], reference[i], 10*tol, "Unexpected solution")) return 1;
    }

    return 0;
}

int
main (void) {
    int stat = 0;

    stat += test_hybrd1 ();

    if (stat > 0) {
        fprintf(stderr, "[FAIL] %d tests failed\n", stat);
    } else {
        fprintf(stderr, "[PASS] all tests passed\n");
    }
    return stat == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
