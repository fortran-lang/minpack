#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "testsuite.h"
#include "minpack.h"

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

int
test_hybrd1 (void)
{
    int n = 9;
    int info = 0;
    int lwa = 180;
    double x[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    double fvec[9];
    double wa[180];
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
test_hybrd (void)
{
    int n = 9;
    int info = 0, nfev = 0;
    double x[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    double diag[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double fvec[9];
    double fjac[9*9];
    double r[45], qtf[9], wa1[9], wa2[9], wa3[9], wa4[9];
    double tol = sqrt(minpack_dpmpar(1));
    double reference[9] = {
        -0.5706545, -0.6816283, -0.7017325,
        -0.7042129, -0.7013690, -0.6918656,
        -0.6657920, -0.5960342, -0.4164121};

    minpack_hybrd(trial_hybrd_fcn, n, x, fvec, tol, 2000, 1, 1, 0.0, diag, 2, 100.0, 0,
            &info, &nfev, fjac, 9, r, 45, qtf, wa1, wa2, wa3, wa4, NULL);

    if (!check(info, 1, "Unexpected info value")) return 1;
    if (!check(nfev, 14, "Unexpected number of function calls")) return 1;

    if (!check(enorm(n, fvec), 0.0, tol, "Unexpected residual")) return 1;

    for(int i = 0; i < 9; i++) {
        if (!check(x[i], reference[i], 10*tol, "Unexpected solution")) return 1;
    }

    return 0;
}

void
trial_lmder_fcn(int m, int n, const double* x, double* fvec, double* fjac,
                int ldfjac, int* iflag, void* data) {
    assert(!!data);
    double* y = (double*)data;
    assert(m == 15);
    assert(n == 3);

    if (*iflag == 1) {
        for (int i = 0; i < m; i++) {
            double tmp1 = i + 1;
            double tmp2 = 16 - i - 1;
            double tmp3 = i >= 8 ? tmp2 : tmp1;
            fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
        }
    } else if (*iflag == 2) {
        assert(ldfjac == m);
        for (int i = 0; i < m; i++) {
            double tmp1 = i + 1;
            double tmp2 = 16 - i - 1;
            double tmp3 = i >= 8 ? tmp2 : tmp1;
            double tmp4 = pow(x[1]*tmp2 + x[2]*tmp3, 2);
            fjac[i] = -1.0;
            fjac[i + ldfjac] = tmp1*tmp2/tmp4;
            fjac[i + 2*ldfjac] = tmp1*tmp3/tmp4;
        }
    }
}

int
test_lmder1 (void)
{
    const double y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1, 3.9e-1,
        3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34e0, 2.1e0, 4.39e0};
    const int m = 15, n = 3;
    int info = 0;
    double x[3] = {1.0, 1.0, 1.0}, xp[n];
    double fvec[m], fvecp[m], err[m];
    double fjac[m*n];
    int ipvt[n];
    double wa[5*n+m];
    double tol = sqrt(minpack_dpmpar(1));

    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 1, err);
    info = 1;
    trial_lmder_fcn(m, n, x, fvec, fjac, m, &info, y);
    info = 2;
    trial_lmder_fcn(m, n, x, fvec, fjac, m, &info, y);
    info = 1;
    trial_lmder_fcn(m, n, xp, fvecp, fjac, m, &info, y);
    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 2, err);

    for (int i = 0; i < 15; i++) {
        if (!check(err[i], 1.0, tol, "Unexpected derivatives")) return 1;
    }

    minpack_lmder1(trial_lmder_fcn, m, n, x, fvec, fjac, m, tol, &info, ipvt, wa, 30, y);
    if (!check(info, 1, "Unexpected info value")) return 1;
    if (!check(x[0], 0.8241058e-1, 100*tol, "Unexpected x[0]")) return 1;
    if (!check(x[1], 0.1133037e+1, 100*tol, "Unexpected x[1]")) return 1;
    if (!check(x[2], 0.2343695e+1, 100*tol, "Unexpected x[2]")) return 1;
    if (!check(enorm(m, fvec), 0.9063596e-1, tol, "Unexpected residual")) return 1;

    return 0;
}

void
trial_lmdif_fcn(int m, int n, const double* x, double* fvec, int* iflag, void* data) {
    assert(!!data);
    double* y = (double*)data;
    assert(m == 15);
    assert(n == 3);
    if (*iflag == 0) return;

    for (int i = 0; i < m; i++) {
        double tmp1 = i + 1;
        double tmp2 = 16 - i - 1;
        double tmp3 = i >= 8 ? tmp2 : tmp1;
        fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
    }
}

int
test_lmdif1 (void)
{
    const int m = 15, n = 3;
    double y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1, 3.9e-1,
        3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34e0, 2.1e0, 4.39e0};
    double x[3] = {1.0, 1.0, 1.0}, fvec[15];
    int info = 0;
    double tol = sqrt(minpack_dpmpar(1));
    int ipvt[n];
    int lwa = m*n + 5*n + m;
    double wa[lwa];

    minpack_lmdif1(trial_lmdif_fcn, 15, 3, x, fvec, tol, &info, ipvt, wa, lwa, y);
    if (!check(info, 1, "Unexpected info value")) return 1;
    if (!check(x[0], 0.8241058e-1, 100*tol, "Unexpected x[0]")) return 1;
    if (!check(x[1], 0.1133037e+1, 100*tol, "Unexpected x[1]")) return 1;
    if (!check(x[2], 0.2343695e+1, 100*tol, "Unexpected x[2]")) return 1;
    if (!check(enorm(m, fvec), 0.9063596e-1, tol, "Unexpected residual")) return 1;

    return 0;
}

int
main (void) {
    int stat = 0;

    stat += run("hybrd1", test_hybrd1);
    stat += run("hybrd ", test_hybrd);
    stat += run("lmder1", test_lmder1);
    stat += run("lmdif1", test_lmdif1);

    if (stat > 0) {
        fprintf(stderr, "[FAIL] %d test(s) failed\n", stat);
    } else {
        fprintf(stderr, "[PASS] all tests passed\n");
    }
    return stat == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
