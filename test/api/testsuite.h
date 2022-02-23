#pragma once

#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#define check(x, ...) \
    _Generic((x), \
        int: check_int, \
     double: check_double \
            )(x, __VA_ARGS__)

static inline bool
check_int(int actual, int expected, const char *msg)
{
    if (expected == actual) {
        return true;
    }
    fprintf(stderr, "FAIL: %s: expected %d, got %d\n", msg, expected, actual);
    return false;
}

static inline bool
check_double(double actual, double expected, double tol, const char *msg)
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

static inline int
run(const char* name, int (*test)(void)) {
    printf("Testing %s ...", name);
    int stat = test();
    if (stat) {
        printf(" FAILED\n");
    } else {
        printf(" OK\n");
    }
    return stat;
}
