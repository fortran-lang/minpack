#pragma once

#ifdef __cplusplus
#define MINPACK_EXTERN extern "C"
#else
#define MINPACK_EXTERN extern
#endif

// In case we have to add some attributes to all functions (e.g. for DLL exports)
#define MINPACK_CALL

MINPACK_EXTERN double MINPACK_CALL
minpack_dpmpar(int /* i */);

typedef void (*minpack_func)(
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    int* /* iflag */,
    void* /* udata */);

#ifdef MINPACK_CFFI
extern "Python" void MINPACK_CALL
func(
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    int* /* iflag */,
    void* /* udata */);
#endif

/*
 * the purpose of hybrd is to find a zero of a system of
 * n nonlinear functions in n variables by a modification
 * of the powell hybrid method. the user must provide a
 * subroutine which calculates the functions. the jacobian is
 * then calculated by a forward-difference approximation.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_hybrd(
    minpack_func /* fcn */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double /* xtol */,
    int /* maxfev */,
    int /* ml */,
    int /* mu */,
    double /* epsfcn */,
    double* /* diag */,
    int /* mode */,
    double /* factor */,
    int /* nprint */,
    int* /* info */,
    int* /* nfev */,
    double* /* fjac */,
    int /* ldfjac */,
    double* /* r */,
    int /* lr */,
    double* /* qtf */,
    double* /* wa1 */,
    double* /* wa2 */,
    double* /* wa3 */,
    double* /* wa4 */,
    void* /* udata */);

/*
 * the purpose of hybrd1 is to find a zero of a system of
 * n nonlinear functions in n variables by a modification
 * of the powell hybrid method. this is done by using the
 * more general nonlinear equation solver hybrd. the user
 * must provide a subroutine which calculates the functions.
 * the jacobian is then calculated by a forward-difference
 * approximation.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_hybrd1(
    minpack_func /* fcn */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double /* tol */,
    int* /* info */,
    double* /* wa */,
    int /* lwa */,
    void* /* udata */);

typedef void (*minpack_fcn_hybrj)(
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    int* /* iflag */,
    void* /* udata */);

#ifdef MINPACK_CFFI
extern "Python" void MINPACK_CALL
fcn_hybrj(
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    int* /* iflag */,
    void* /* udata */);
#endif

/*
 * the purpose of hybrj is to find a zero of a system of
 * n nonlinear functions in n variables by a modification
 * of the powell hybrid method. the user must provide a
 * subroutine which calculates the functions and the jacobian.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_hybrj(
    minpack_fcn_hybrj /* fcn */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    double /* xtol */,
    int /* maxfev */,
    double* /* diag */,
    int /* mode */,
    double /* factor */,
    int /* nprint */,
    int* /* info */,
    int* /* nfev */,
    int* /* njev */,
    double* /* r */,
    int /* lr */,
    double* /* qtf */,
    double* /* wa1 */,
    double* /* wa2 */,
    double* /* wa3 */,
    double* /* wa4 */,
    void* /* udata */);

/*
 * The purpose of hybrj1 is to find a zero of a system of
 * n nonlinear functions in n variables by a modification
 * of the powell hybrid method. this is done by using the
 * more general nonlinear equation solver hybrj. the user
 * must provide a subroutine which calculates the functions
 * and the jacobian.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_hybrj1(
    minpack_fcn_hybrj /* fcn */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    double /* tol */,
    int* /* info */,
    double* /* wa */,
    int /* lwa */,
    void* /* udata */);

typedef void (*minpack_fcn_lmder)(
    int /* m */,
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    int* /* iflag */,
    void* /* udata */);

#ifdef MINPACK_CFFI
extern "Python" void MINPACK_CALL
fcn_lmder(
    int /* m */,
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    int* /* iflag */,
    void* /* udata */);
#endif

/*
 * the purpose of lmder is to minimize the sum of the squares of
 * m nonlinear functions in n variables by a modification of
 * the levenberg-marquardt algorithm. the user must provide a
 * subroutine which calculates the functions and the jacobian.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_lmder(
    minpack_fcn_lmder /* fcn */,
    int /* m */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    double /* ftol */,
    double /* xtol */,
    double /* gtol */,
    int /* maxfev */,
    double* /* diag */,
    int /* mode */,
    double /* factor */,
    int /* nprint */,
    int* /* info */,
    int* /* nfev */,
    int* /* njev */,
    int* /* ipvt */,
    double* /* qtf */,
    double* /* wa1 */,
    double* /* wa2 */,
    double* /* wa3 */,
    double* /* wa4 */,
    void* /* udata */);

/*
 * the purpose of lmder1 is to minimize the sum of the squares of
 * m nonlinear functions in n variables by a modification of the
 * levenberg-marquardt algorithm. this is done by using the more
 * general least-squares solver lmder. the user must provide a
 * subroutine which calculates the functions and the jacobian.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_lmder1(
    minpack_fcn_lmder /* fcn */,
    int /* m */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    double /* tol */,
    int* /* info */,
    int* /* ipvt */,
    double* /* wa */,
    int /* lwa */,
    void* /* udata */);

typedef void (*minpack_func2)(
    int /* m */,
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    int* /* iflag */,
    void* /* udata */);

#ifdef MINPACK_CFFI
extern "Python" void MINPACK_CALL
func2(
    int /* m */,
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    int* /* iflag */,
    void* /* udata */);
#endif

/*
 * the purpose of lmdif is to minimize the sum of the squares of
 * m nonlinear functions in n variables by a modification of
 * the levenberg-marquardt algorithm. the user must provide a
 * subroutine which calculates the functions. the jacobian is
 * then calculated by a forward-difference approximation.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_lmdif(
    minpack_func2 /* fcn */,
    int /* m */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double /* ftol */,
    double /* xtol */,
    double /* gtol */,
    int /* maxfev */,
    double /* epsfcn */,
    double* /* diag */,
    int /* mode */,
    double /* factor */,
    int /* nprint */,
    int* /* info */,
    int* /* nfev */,
    double* /* fjac */,
    int /* ldfjac */,
    int* /* ipvt */,
    double* /* qtf */,
    double* /* wa1 */,
    double* /* wa2 */,
    double* /* wa3 */,
    double* /* wa4 */,
    void* /* udata */);

/*
 * the purpose of lmdif1 is to minimize the sum of the squares of
 * m nonlinear functions in n variables by a modification of the
 * levenberg-marquardt algorithm. this is done by using the more
 * general least-squares solver lmdif. the user must provide a
 * subroutine which calculates the functions. the jacobian is
 * then calculated by a forward-difference approximation.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_lmdif1(
    minpack_func2 /* fcn */,
    int /* m */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double /* tol */,
    int* /* info */,
    int* /* iwa */,
    double* /* wa */,
    int /* lwa */,
    void* /* udata */);

typedef void (*minpack_fcn_lmstr)(
    int /* m */,
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    double* /* fjrow */,
    int* /* iflag */,
    void* /* udata */);

#ifdef MINPACK_CFFI
extern "Python" void MINPACK_CALL
fcn_lmstr(
    int /* m */,
    int /* n */,
    const double* /* x */,
    double* /* fvec */,
    double* /* fjrow */,
    int* /* iflag */,
    void* /* udata */);
#endif

/*
 * the purpose of lmstr is to minimize the sum of the squares of
 * m nonlinear functions in n variables by a modification of
 * the levenberg-marquardt algorithm which uses minimal storage.
 * the user must provide a subroutine which calculates the
 * functions and the rows of the jacobian.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_lmstr(
    minpack_fcn_lmstr /* fcn */,
    int /* m */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    double /* ftol */,
    double /* xtol */,
    double /* gtol */,
    int /* maxfev */,
    double* /* diag */,
    int /* mode */,
    double /* factor */,
    int /* nprint */,
    int* /* info */,
    int* /* nfev */,
    int* /* njev */,
    int* /* ipvt */,
    double* /* qtf */,
    double* /* wa1 */,
    double* /* wa2 */,
    double* /* wa3 */,
    double* /* wa4 */,
    void* /* udata */);

/*
 * the purpose of lmstr1 is to minimize the sum of the squares of
 * m nonlinear functions in n variables by a modification of
 * the levenberg-marquardt algorithm which uses minimal storage.
 * this is done by using the more general least-squares solver
 * lmstr. the user must provide a subroutine which calculates
 * the functions and the rows of the jacobian.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_lmstr1(
    minpack_fcn_lmstr /* fcn */,
    int /* m */,
    int /* n */,
    double* /* x */,
    double* /* fvec */,
    double* /* fjac */,
    int /* ldfjac */,
    double /* tol */,
    int* /* info */,
    int* /* ipvt */,
    double* /* wa */,
    int /* lwa */,
    void* /* udata */);

/*
 * this subroutine checks the gradients of m nonlinear functions
 * in n variables, evaluated at a point x, for consistency with
 * the functions themselves.
 *
 * the subroutine does not perform reliably if cancellation or
 * rounding errors cause a severe loss of significance in the
 * evaluation of a function. therefore, none of the components
 * of x should be unusually small (in particular, zero) or any
 * other value which may cause loss of significance.
 */
MINPACK_EXTERN void MINPACK_CALL
minpack_chkder(
    int /* m */,
    int /* n */,
    const double* /* x */,
    const double* /* fvec */,
    const double* /* fjac */,
    int /* ldfjac */,
    double* /* xp */,
    const double* /* fvecp */,
    int /* mode */,
    double* /* err */);
