!*****************************************************************************************
!>
!  Modernized Minpack
!
!### Authors
!  * argonne national laboratory. minpack project. march 1980.
!    burton s. garbow, kenneth e. hillstrom, jorge j. more.
!  * Jacob Williams, Sept 2021, updated to modern standards.

module minpack_legacy
    use minpack_module, only : &
        & hybrd_ => hybrd, hybrd1_ => hybrd1, func_ => func, &
        & hybrj_ => hybrj, hybrj1_ => hybrj1, fcn_hybrj_ => fcn_hybrj, &
        & lmdif_ => lmdif, lmdif1_ => lmdif1, func2_ => func2, &
        & lmder_ => lmder, lmder1_ => lmder1, fcn_lmder_ => fcn_lmder, &
        & lmstr_ => lmstr, lmstr1_ => lmstr1, fcn_lmstr_ => fcn_lmstr, &
        & chkder_ => chkder, wp
    implicit none
    private

    public :: &
        & hybrd, hybrd1, func, &
        & hybrj, hybrj1, fcn_hybrj, &
        & lmdif, lmdif1, func2, &
        & lmder, lmder1, fcn_lmder, &
        & lmstr, lmstr1, fcn_lmstr, &
        & chkder

    abstract interface
        subroutine func(n, x, fvec, iflag)
            !! user-supplied subroutine for [[hybrd]], [[hybrd1]], and [[fdjac1]]
            import :: wp
            implicit none
            integer, intent(in) :: n !! the number of variables.
            real(wp), intent(in) :: x(n) !! independent variable vector
            real(wp), intent(out) :: fvec(n) !! value of function at `x`
            integer, intent(inout) :: iflag !! set to <0 to terminate execution
        end subroutine func

        subroutine func2(m, n, x, fvec, iflag)
            !! user-supplied subroutine for [[fdjac2]], [[lmdif]], and [[lmdif1]]
            import :: wp
            implicit none
            integer, intent(in) :: m !! the number of functions.
            integer, intent(in) :: n !! the number of variables.
            real(wp), intent(in) :: x(n) !! independent variable vector
            real(wp), intent(out) :: fvec(m) !! value of function at `x`
            integer, intent(inout) :: iflag !! the value of iflag should not be changed unless
                                           !! the user wants to terminate execution of lmdif.
                                           !! in this case set iflag to a negative integer.
        end subroutine func2

        subroutine fcn_hybrj(n, x, fvec, fjac, ldfjac, iflag)
            !! user-supplied subroutine for [[hybrj]] and [[hybrj1]]
            import :: wp
            implicit none
            integer, intent(in) :: n !! the number of variables.
            real(wp), dimension(n), intent(in) :: x !! independent variable vector
            integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
            real(wp), dimension(n), intent(inout) :: fvec !! value of function at `x`
            real(wp), dimension(ldfjac, n), intent(inout) :: fjac !! jacobian matrix at `x`
            integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                            !! return this vector in fvec. do not alter fjac.
                                            !! if iflag = 2 calculate the jacobian at x and
                                            !! return this matrix in fjac. do not alter fvec.
                                            !!
                                            !! the value of iflag should not be changed by fcn unless
                                            !! the user wants to terminate execution of hybrj.
                                            !! in this case set iflag to a negative integer.
        end subroutine fcn_hybrj

        subroutine fcn_lmder(m, n, x, fvec, fjac, ldfjac, iflag)
            !! user-supplied subroutine for [[lmder]] and [[lmder1]]
            import :: wp
            implicit none
            integer, intent(in) :: m !! the number of functions.
            integer, intent(in) :: n !! the number of variables.
            integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
            integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                           !! return this vector in fvec. do not alter fjac.
                                           !! if iflag = 2 calculate the jacobian at x and
                                           !! return this matrix in fjac. do not alter fvec.
                                           !!
                                           !! the value of iflag should not be changed by fcn unless
                                           !! the user wants to terminate execution of lmder.
                                           !! in this case set iflag to a negative integer.
            real(wp), intent(in) :: x(n) !! independent variable vector
            real(wp), intent(inout) :: fvec(m) !! value of function at `x`
            real(wp), intent(inout) :: fjac(ldfjac, n) !! jacobian matrix at `x`
        end subroutine fcn_lmder

        subroutine fcn_lmstr(m, n, x, fvec, fjrow, iflag)
            import :: wp
            implicit none
            integer, intent(in) :: m !! the number of functions.
            integer, intent(in) :: n !! the number of variables.
            integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                       !! return this vector in fvec.
                                       !! if iflag = i calculate the (i-1)-st row of the
                                       !! jacobian at x and return this vector in fjrow.
                                       !!
                                       !! the value of iflag should not be changed by fcn unless
                                       !! the user wants to terminate execution of lmstr.
                                       !! in this case set iflag to a negative integer.
            real(wp), intent(in) :: x(n) !! independent variable vector
            real(wp), intent(inout) :: fvec(m) !! value of function at `x`
            real(wp), intent(inout) :: fjrow(n) !! jacobian row
        end subroutine fcn_lmstr

    end interface


    type :: hybrd_data
        procedure(func), pointer, nopass :: fcn => null()
    end type hybrd_data

    type :: hybrj_data
        procedure(fcn_hybrj), pointer, nopass :: fcn => null()
    end type hybrj_data

    type :: lmdif_data
        procedure(func2), pointer, nopass :: fcn => null()
    end type lmdif_data

    type :: lmder_data
        procedure(fcn_lmder), pointer, nopass :: fcn => null()
    end type lmder_data

    type :: lmstr_data
        procedure(fcn_lmstr), pointer, nopass :: fcn => null()
    end type lmstr_data

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  this subroutine checks the gradients of m nonlinear functions
!  in n variables, evaluated at a point x, for consistency with
!  the functions themselves.
!
!  the subroutine does not perform reliably if cancellation or
!  rounding errors cause a severe loss of significance in the
!  evaluation of a function. therefore, none of the components
!  of x should be unusually small (in particular, zero) or any
!  other value which may cause loss of significance.

    subroutine chkder(m, n, x, Fvec, Fjac, Ldfjac, Xp, Fvecp, Mode, Err)

        implicit none

        integer, intent(in) :: m !! a positive integer input variable set to the number
                            !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of variables.
        integer, intent(in) :: Ldfjac !! a positive integer input parameter not less than m
                                    !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: Mode !! an integer input variable set to 1 on the first call
                                !! and 2 on the second. other values of mode are equivalent
                                !! to mode = 1.
                                !!
                                !! the user must call chkder twice,
                                !! first with mode = 1 and then with mode = 2.
                                !!
                                !!  * mode = 1. **on input**, x must contain the point of evaluation.
                                !!    **on output**, xp is set to a neighboring point.
                                !!
                                !!  * mode = 2. **on input**, fvec must contain the functions and the
                                !!    rows of fjac must contain the gradients
                                !!    of the respective functions each evaluated
                                !!    at x, and fvecp must contain the functions
                                !!    evaluated at xp.
                                !!    **on output**, err contains measures of correctness of
                                !!    the respective gradients.
        real(wp), intent(in) :: x(n) !! input array
        real(wp), intent(in) :: Fvec(m) !! an array of length m. on input when mode = 2,
                                    !! fvec must contain the functions evaluated at x.
        real(wp), intent(in) :: Fjac(Ldfjac, n) !! an m by n array. on input when mode = 2,
                                            !! the rows of fjac must contain the gradients of
                                            !! the respective functions evaluated at x.
        real(wp), intent(out) :: Xp(n) !! an array of length n. on output when mode = 1,
                                    !! xp is set to a neighboring point of x.
        real(wp), intent(in) :: Fvecp(m) !! an array of length m. on input when mode = 2,
                                    !! fvecp must contain the functions evaluated at xp.
        real(wp), intent(out) :: Err(m) !! an array of length m. on output when mode = 2,
                                    !! err contains measures of correctness of the respective
                                    !! gradients. if there is no severe loss of significance,
                                    !! then if err(i) is 1.0 the i-th gradient is correct,
                                    !! while if err(i) is 0.0 the i-th gradient is incorrect.
                                    !! for values of err between 0.0 and 1.0, the categorization
                                    !! is less certain. in general, a value of err(i) greater
                                    !! than 0.5 indicates that the i-th gradient is probably
                                    !! correct, while a value of err(i) less than 0.5 indicates
                                    !! that the i-th gradient is probably incorrect.

        call chkder_(m, n, x, Fvec, Fjac, Ldfjac, Xp, Fvecp, Mode, Err)
    end subroutine chkder
!*****************************************************************************************


!*****************************************************************************************
!>
!  the purpose of hybrd is to find a zero of a system of
!  n nonlinear functions in n variables by a modification
!  of the powell hybrid method. the user must provide a
!  subroutine which calculates the functions. the jacobian is
!  then calculated by a forward-difference approximation.

    subroutine hybrd(fcn, n, x, Fvec, Xtol, Maxfev, Ml, Mu, Epsfcn, Diag, Mode, &
                     Factor, Nprint, Info, Nfev, Fjac, Ldfjac, r, Lr, Qtf, Wa1, &
                     Wa2, Wa3, Wa4)
        procedure(func) :: fcn                  !! user-supplied subroutine which calculates the functions
        integer, intent(in) :: n                 !! a positive integer input variable set to the number
                                            !! of functions and variables.
        integer, intent(in) :: maxfev            !! a positive integer input variable. termination
                                            !! occurs when the number of calls to `fcn` is at least `maxfev`
                                            !! by the end of an iteration.
        integer, intent(in) :: ml                !! a nonnegative integer input variable which specifies
                                            !! the number of subdiagonals within the band of the
                                            !! jacobian matrix. if the jacobian is not banded, set
                                            !! `ml` to at least `n - 1`.
        integer, intent(in) :: mu                !! a nonnegative integer input variable which specifies
                                            !! the number of superdiagonals within the band of the
                                            !! jacobian matrix. if the jacobian is not banded, set
                                            !! `mu` to at least` n - 1`.
        integer, intent(in) :: mode              !! if `mode = 1`, the
                                            !! variables will be scaled internally. if `mode = 2`,
                                            !! the scaling is specified by the input `diag`. other
                                            !! values of `mode` are equivalent to `mode = 1`.
        integer, intent(in)  :: nprint           !! an integer input variable that enables controlled
                                            !! printing of iterates if it is positive. in this case,
                                            !! `fcn` is called with `iflag = 0` at the beginning of the first
                                            !! iteration and every `nprint` iterations thereafter and
                                            !! immediately prior to return, with `x` and `fvec` available
                                            !! for printing. if `nprint` is not positive, no special calls
                                            !! of `fcn` with `iflag = 0` are made.
        integer, intent(out) :: info             !! an integer output variable. if the user has
                                            !! terminated execution, `info` is set to the (negative)
                                            !! value of `iflag`. see description of `fcn`. otherwise,
                                            !! `info` is set as follows:
                                            !!
                                            !!  * ***info = 0*** improper input parameters.
                                            !!  * ***info = 1*** relative error between two consecutive iterates
                                            !!    is at most `xtol`.
                                            !!  * ***info = 2*** number of calls to `fcn` has reached or exceeded
                                            !!    `maxfev`.
                                            !!  * ***info = 3*** `xtol` is too small. no further improvement in
                                            !!    the approximate solution `x` is possible.
                                            !!  * ***info = 4*** iteration is not making good progress, as
                                            !!    measured by the improvement from the last
                                            !!    five jacobian evaluations.
                                            !!  * ***info = 5*** iteration is not making good progress, as
                                            !!    measured by the improvement from the last
                                            !!    ten iterations.
        integer, intent(out) :: nfev             !! output variable set to the number of calls to `fcn`.
        integer, intent(in):: ldfjac             !! a positive integer input variable not less than `n`
                                            !! which specifies the leading dimension of the array `fjac`.
        integer, intent(in) :: lr                !! a positive integer input variable not less than `(n*(n+1))/2`.
        real(wp), intent(in) :: xtol             !! a nonnegative input variable. termination
                                            !! occurs when the relative error between two consecutive
                                            !! iterates is at most `xtol`.
        real(wp), intent(in) :: epsfcn           !! an input variable used in determining a suitable
                                            !! step length for the forward-difference approximation. this
                                            !! approximation assumes that the relative errors in the
                                            !! functions are of the order of `epsfcn`. if `epsfcn` is less
                                            !! than the machine precision, it is assumed that the relative
                                            !! errors in the functions are of the order of the machine
                                            !! precision.
        real(wp), intent(in) :: factor           !! a positive input variable used in determining the
                                            !! initial step bound. this bound is set to the product of
                                            !! `factor` and the euclidean norm of `diag*x` if nonzero, or else
                                            !! to `factor` itself. in most cases factor should lie in the
                                            !! interval (.1,100.). 100. is a generally recommended value.
        real(wp), intent(inout) :: x(n)          !! array of length n. on input `x` must contain
                                            !! an initial estimate of the solution vector. on output `x`
                                            !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: fvec(n)         !! an output array of length `n` which contains
                                            !! the functions evaluated at the output `x`.
        real(wp), intent(inout) :: diag(n)       !! an array of length `n`. if `mode = 1` (see
                                            !! below), `diag` is internally set. if `mode = 2`, `diag`
                                            !! must contain positive entries that serve as
                                            !! multiplicative scale factors for the variables.
        real(wp), intent(out) :: fjac(ldfjac, n)  !! array which contains the
                                            !! orthogonal matrix `q` produced by the QR factorization
                                            !! of the final approximate jacobian.
        real(wp), intent(out) :: r(lr)           !! an output array which contains the
                                            !! upper triangular matrix produced by the QR factorization
                                            !! of the final approximate jacobian, stored rowwise.
        real(wp), intent(out) :: qtf(n)          !! an output array of length `n` which contains
                                            !! the vector `(q transpose)*fvec`.
        real(wp), intent(inout) :: wa1(n)  !! work array
        real(wp), intent(inout) :: wa2(n)  !! work array
        real(wp), intent(inout) :: wa3(n)  !! work array
        real(wp), intent(inout) :: wa4(n)  !! work array

        type(hybrd_data) :: wrapper
        wrapper%fcn => fcn

        call hybrd_(wrap_func, n, x, Fvec, Xtol, Maxfev, Ml, Mu, Epsfcn, Diag, Mode, &
                    Factor, Nprint, Info, Nfev, Fjac, Ldfjac, r, Lr, Qtf, Wa1, &
                    Wa2, Wa3, Wa4, wrapper)
    end subroutine hybrd
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of hybrd1 is to find a zero of a system of
!  n nonlinear functions in n variables by a modification
!  of the powell hybrid method. this is done by using the
!  more general nonlinear equation solver hybrd. the user
!  must provide a subroutine which calculates the functions.
!  the jacobian is then calculated by a forward-difference
!  approximation.

    subroutine hybrd1(fcn, n, x, Fvec, Tol, Info, Wa, Lwa)

        implicit none

        procedure(func)                     :: fcn      !! user-supplied subroutine which calculates the functions
        integer, intent(in)                  :: n        !! a positive integer input variable set to the number
                                                    !! of functions and variables.
        integer, intent(out)                 :: info     !! an integer output variable. if the user has
                                                    !! terminated execution, info is set to the (negative)
                                                    !! value of `iflag`. see description of `fcn`. otherwise,
                                                    !! `info` is set as follows:
                                                    !!
                                                    !!  * ***info = 0*** improper input parameters.
                                                    !!  * ***info = 1*** algorithm estimates that the relative error
                                                    !!  between `x` and the solution is at most `tol`.
                                                    !!  * ***info = 2*** number of calls to `fcn` has reached or exceeded
                                                    !!  `200*(n+1)`.
                                                    !!  * ***info = 3*** `tol` is too small. no further improvement in
                                                    !!  the approximate solution `x` is possible.
                                                    !!  * ***info = 4*** iteration is not making good progress.
        real(wp), intent(in)                 :: tol      !! a nonnegative input variable. termination occurs
                                                    !! when the algorithm estimates that the relative error
                                                    !! between `x` and the solution is at most `tol`.
        real(wp), dimension(n), intent(inout) :: x        !! an array of length `n`. on input `x` must contain
                                                    !! an initial estimate of the solution vector. on output `x`
                                                    !! contains the final estimate of the solution vector.
        real(wp), dimension(n), intent(out)   :: fvec     !! an output array of length `n` which contains
                                                    !! the functions evaluated at the output `x`.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than
                              !! (n*(3*n+13))/2.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

        type(hybrd_data) :: wrapper
        wrapper%fcn => fcn

        call hybrd1_(wrap_func, n, x, Fvec, Tol, Info, Wa, Lwa, wrapper)
    end subroutine hybrd1
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of hybrj is to find a zero of a system of
!  n nonlinear functions in n variables by a modification
!  of the powell hybrid method. the user must provide a
!  subroutine which calculates the functions and the jacobian.

    subroutine hybrj(fcn, n, x, Fvec, Fjac, Ldfjac, Xtol, Maxfev, Diag, Mode, &
                     Factor, Nprint, Info, Nfev, Njev, r, Lr, Qtf, Wa1, Wa2, &
                     Wa3, Wa4)

        implicit none

        procedure(fcn_hybrj) :: fcn !! the user-supplied subroutine which
                                !! calculates the functions and the jacobian
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of functions and variables.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
                                    !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
                                    !! occurs when the number of calls to fcn with iflag = 1
                                    !! has reached maxfev.
        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
                                !! variables will be scaled internally. if mode = 2,
                                !! the scaling is specified by the input diag. other
                                !! values of mode are equivalent to mode = 1.
        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
                                    !! printing of iterates if it is positive. in this case,
                                    !! fcn is called with iflag = 0 at the beginning of the first
                                    !! iteration and every nprint iterations thereafter and
                                    !! immediately prior to return, with x and fvec available
                                    !! for printing. fvec and fjac should not be altered.
                                    !! if nprint is not positive, no special calls of fcn
                                    !! with iflag = 0 are made.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                !! terminated execution, info is set to the (negative)
                                !! value of iflag. see description of fcn. otherwise,
                                !! info is set as follows:
                                !!
                                !!  * ***info = 0***   improper input parameters.
                                !!  * ***info = 1***   relative error between two consecutive iterates
                                !!    is at most xtol.
                                !!  * ***info = 2***   number of calls to fcn with iflag = 1 has
                                !!    reached maxfev.
                                !!  * ***info = 3***   xtol is too small. no further improvement in
                                !!    the approximate solution x is possible.
                                !!  * ***info = 4***   iteration is not making good progress, as
                                !!    measured by the improvement from the last
                                !!    five jacobian evaluations.
                                !!  * ***info = 5***   iteration is not making good progress, as
                                !!    measured by the improvement from the last
                                !!    ten iterations.
        integer, intent(out) :: Nfev !! an integer output variable set to the number of
                                !! calls to fcn with iflag = 1.
        integer, intent(out) :: Njev !! an integer output variable set to the number of
                                !! calls to fcn with iflag = 2.
        integer, intent(in) :: Lr !! a positive integer input variable not less than
                                !! (n*(n+1))/2.
        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
                                !! occurs when the relative error between two consecutive
                                !! iterates is at most xtol.
        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
                                    !! initial step bound. this bound is set to the product of
                                    !! factor and the euclidean norm of diag*x if nonzero, or else
                                    !! to factor itself. in most cases factor should lie in the
                                    !! interval (.1,100.). 100. is a generally recommended value.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                    !! an initial estimate of the solution vector. on output x
                                    !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(n) !! an output array of length n which contains
                                    !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array which contains the
                                            !! orthogonal matrix q produced by the qr factorization
                                            !! of the final approximate jacobian.
        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
                                        !! below), diag is internally set. if mode = 2, diag
                                        !! must contain positive entries that serve as
                                        !! multiplicative scale factors for the variables.
        real(wp), intent(out) :: r(Lr) !! an output array of length lr which contains the
                                    !! upper triangular matrix produced by the qr factorization
                                    !! of the final approximate jacobian, stored rowwise.
        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
                                    !! the vector (q transpose)*fvec.
        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
        real(wp), intent(inout) :: Wa4(n) !! work array of length n.

        type(hybrj_data) :: wrapper
        wrapper%fcn => fcn

        call hybrj_(wrap_fcn_hybrj, n, x, Fvec, Fjac, Ldfjac, Xtol, Maxfev, Diag, Mode, &
                    Factor, Nprint, Info, Nfev, Njev, r, Lr, Qtf, Wa1, Wa2, &
                    Wa3, Wa4, wrapper)
    end subroutine hybrj
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of hybrj1 is to find a zero of a system of
!  n nonlinear functions in n variables by a modification
!  of the powell hybrid method. this is done by using the
!  more general nonlinear equation solver hybrj. the user
!  must provide a subroutine which calculates the functions
!  and the jacobian.

    subroutine hybrj1(fcn, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Wa, Lwa)

        implicit none

        procedure(fcn_hybrj) :: fcn !! the user-supplied subroutine which
                                !! calculates the functions and the jacobian
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of functions and variables.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
                                 !! which specifies the leading dimension of the array fjac.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                !! terminated execution, info is set to the (negative)
                                !! value of iflag. see description of fcn. otherwise,
                                !! info is set as follows:
                                !!
                                !!  * ***info = 0***   improper input parameters.
                                !!  * ***info = 1***   algorithm estimates that the relative error
                                !!    between x and the solution is at most tol.
                                !!  * ***info = 2***   number of calls to fcn with iflag = 1 has
                                !!    reached 100*(n+1).
                                !!  * ***info = 3***   tol is too small. no further improvement in
                                !!    the approximate solution x is possible.
                                !!  * ***info = 4***   iteration is not making good progress.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than
                              !! (n*(n+13))/2.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                !! when the algorithm estimates that the relative error
                                !! between x and the solution is at most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                    !! an initial estimate of the solution vector. on output x
                                    !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(n) !! an output array of length n which contains
                                    !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array which contains the
                                            !! orthogonal matrix q produced by the qr factorization
                                            !! of the final approximate jacobian.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

        type(hybrj_data) :: wrapper
        wrapper%fcn => fcn

        call hybrj1_(wrap_fcn_hybrj, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Wa, Lwa, wrapper)
    end subroutine hybrj1
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of lmder is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm. the user must provide a
!  subroutine which calculates the functions and the jacobian.

    subroutine lmder(fcn, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
                     Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
                     Wa1, Wa2, Wa3, Wa4)

        implicit none

        procedure(fcn_lmder) :: fcn !! the user-supplied subroutine which
                                !! calculates the functions and the jacobian
        integer, intent(in) :: m !! a positive integer input variable set to the number
                            !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                 !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
                                 !! occurs when the number of calls to fcn with iflag = 1
                                 !! has reached maxfev.
        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
                                !! variables will be scaled internally. if mode = 2,
                                !! the scaling is specified by the input diag. other
                                !! values of mode are equivalent to mode = 1.
        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
                                 !! printing of iterates if it is positive. in this case,
                                 !! fcn is called with iflag = 0 at the beginning of the first
                                 !! iteration and every nprint iterations thereafter and
                                 !! immediately prior to return, with x, fvec, and fjac
                                 !! available for printing. fvec and fjac should not be
                                 !! altered. if nprint is not positive, no special calls
                                 !! of fcn with iflag = 0 are made.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                !! terminated execution, info is set to the (negative)
                                !! value of iflag. see description of fcn. otherwise,
                                !! info is set as follows:
                                !!
                                !!  * ***info = 0***  improper input parameters.
                                !!  * ***info = 1***  both actual and predicted relative reductions
                                !!    in the sum of squares are at most ftol.
                                !!  * ***info = 2***  relative error between two consecutive iterates
                                !!    is at most xtol.
                                !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                !!  * ***info = 4***  the cosine of the angle between fvec and any
                                !!    column of the jacobian is at most gtol in
                                !!    absolute value.
                                !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                !!    reached maxfev.
                                !!  * ***info = 6***  ftol is too small. no further reduction in
                                !!    the sum of squares is possible.
                                !!  * ***info = 7***  xtol is too small. no further improvement in
                                !!    the approximate solution x is possible.
                                !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
                                !!    columns of the jacobian to machine precision.
        integer, intent(out) :: Nfev !! an integer output variable set to the number of
                                !! calls to fcn with iflag = 1.
        integer, intent(out) :: Njev !! an integer output variable set to the number of
                                !! calls to fcn with iflag = 2.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                   !! defines a permutation matrix p such that jac*p = q*r,
                                   !! where jac is the final calculated jacobian, q is
                                   !! orthogonal (not stored), and r is upper triangular
                                   !! with diagonal elements of nonincreasing magnitude.
                                   !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
                                !! occurs when both the actual and predicted relative
                                !! reductions in the sum of squares are at most ftol.
                                !! therefore, ftol measures the relative error desired
                                !! in the sum of squares.
        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
                                !! occurs when the relative error between two consecutive
                                !! iterates is at most xtol. therefore, xtol measures the
                                !! relative error desired in the approximate solution.
        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
                                !! occurs when the cosine of the angle between fvec and
                                !! any column of the jacobian is at most gtol in absolute
                                !! value. therefore, gtol measures the orthogonality
                                !! desired between the function vector and the columns
                                !! of the jacobian.
        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
                                  !! initial step bound. this bound is set to the product of
                                  !! factor and the euclidean norm of diag*x if nonzero, or else
                                  !! to factor itself. in most cases factor should lie in the
                                  !! interval (.1,100.).100. is a generally recommended value.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                   !! an initial estimate of the solution vector. on output x
                                   !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                    !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
                                            !! of fjac contains an upper triangular matrix r with
                                            !! diagonal elements of nonincreasing magnitude such that
                                            !!```
                                            !!        t     t           t
                                            !!       p *(jac *jac)*p = r *r,
                                            !!```
                                            !! where p is a permutation matrix and jac is the final
                                            !! calculated jacobian. column j of p is column ipvt(j)
                                            !! (see below) of the identity matrix. the lower trapezoidal
                                            !! part of fjac contains information generated during
                                            !! the computation of r.
        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
                                      !! below), diag is internally set. if mode = 2, diag
                                      !! must contain positive entries that serve as
                                      !! multiplicative scale factors for the variables.
        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
                                   !! the first n elements of the vector (q transpose)*fvec.
        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
        real(wp), intent(inout) :: Wa4(m) !! work array of length m.

        type(lmder_data) :: wrapper
        wrapper%fcn => fcn

        call lmder_(wrap_fcn_lmder, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
                    Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
                    Wa1, Wa2, Wa3, Wa4, wrapper)
    end subroutine lmder
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of lmder1 is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of the
!  levenberg-marquardt algorithm. this is done by using the more
!  general least-squares solver lmder. the user must provide a
!  subroutine which calculates the functions and the jacobian.

    subroutine lmder1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa)
        implicit none

        procedure(fcn_lmder) :: fcn !! user-supplied subroutine which
                                    !! calculates the functions and the jacobian.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                     !! which specifies the leading dimension of the array fjac.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows.
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  algorithm estimates that the relative error
                                    !!    in the sum of squares is at most tol.
                                    !!  * ***info = 2***  algorithm estimates that the relative error
                                    !!    between x and the solution is at most tol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
                                    !!    jacobian to machine precision.
                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                    !!    reached 100*(n+1).
                                    !!  * ***info = 6***  tol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  tol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than 5*n+m.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                       !! defines a permutation matrix p such that jac*p = q*r,
                                       !! where jac is the final calculated jacobian, q is
                                       !! orthogonal (not stored), and r is upper triangular
                                       !! with diagonal elements of nonincreasing magnitude.
                                       !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                   !! when the algorithm estimates either that the relative
                                   !! error in the sum of squares is at most tol or that
                                   !! the relative error between x and the solution is at
                                   !! most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
                                                !! of fjac contains an upper triangular matrix r with
                                                !! diagonal elements of nonincreasing magnitude such that
                                                !!```
                                                !!        t     t           t
                                                !!       p *(jac *jac)*p = r *r,
                                                !!```
                                                !! where p is a permutation matrix and jac is the final
                                                !! calculated jacobian. column j of p is column ipvt(j)
                                                !! (see below) of the identity matrix. the lower trapezoidal
                                                !! part of fjac contains information generated during
                                                !! the computation of r.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

        type(lmder_data) :: wrapper
        wrapper%fcn => fcn

        call lmder1_(wrap_fcn_lmder, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa, &
            & wrapper)
    end subroutine lmder1
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of lmdif is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm. the user must provide a
!  subroutine which calculates the functions. the jacobian is
!  then calculated by a forward-difference approximation.

    subroutine lmdif(fcn, m, n, x, Fvec, Ftol, Xtol, Gtol, Maxfev, Epsfcn, Diag, &
                     Mode, Factor, Nprint, Info, Nfev, Fjac, Ldfjac, Ipvt, &
                     Qtf, Wa1, Wa2, Wa3, Wa4)
        implicit none

        procedure(func2) :: fcn !! the user-supplied subroutine which
                                !! calculates the functions.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
                                     !! occurs when the number of calls to fcn is at least
                                     !! maxfev by the end of an iteration.
        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
                                   !! variables will be scaled internally. if mode = 2,
                                   !! the scaling is specified by the input diag. other
                                   !! values of mode are equivalent to mode = 1.
        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
                                     !! printing of iterates if it is positive. in this case,
                                     !! fcn is called with iflag = 0 at the beginning of the first
                                     !! iteration and every nprint iterations thereafter and
                                     !! immediately prior to return, with x and fvec available
                                     !! for printing. if nprint is not positive, no special calls
                                     !! of fcn with iflag = 0 are made.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows:
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  both actual and predicted relative reductions
                                    !!    in the sum of squares are at most ftol.
                                    !!  * ***info = 2***  relative error between two consecutive iterates
                                    !!    is at most xtol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  the cosine of the angle between fvec and any
                                    !!    column of the jacobian is at most gtol in
                                    !!    absolute value.
                                    !!  * ***info = 5***  number of calls to fcn has reached or
                                    !!    exceeded maxfev.
                                    !!  * ***info = 6***  ftol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  xtol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
                                    !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
                                    !!    columns of the jacobian to machine precision.
        integer, intent(out) :: Nfev !! an integer output variable set to the number of
                                    !! calls to fcn.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                     !! which specifies the leading dimension of the array fjac.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                       !! defines a permutation matrix p such that jac*p = q*r,
                                       !! where jac is the final calculated jacobian, q is
                                       !! orthogonal (not stored), and r is upper triangular
                                       !! with diagonal elements of nonincreasing magnitude.
                                       !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
                                    !! occurs when both the actual and predicted relative
                                    !! reductions in the sum of squares are at most ftol.
                                    !! therefore, ftol measures the relative error desired
                                    !! in the sum of squares.
        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
                                    !! occurs when the relative error between two consecutive
                                    !! iterates is at most xtol. therefore, xtol measures the
                                    !! relative error desired in the approximate solution.
        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
                                    !! occurs when the cosine of the angle between fvec and
                                    !! any column of the jacobian is at most gtol in absolute
                                    !! value. therefore, gtol measures the orthogonality
                                    !! desired between the function vector and the columns
                                    !! of the jacobian.
        real(wp), intent(in) :: Epsfcn !! an input variable used in determining a suitable
                                      !! step length for the forward-difference approximation. this
                                      !! approximation assumes that the relative errors in the
                                      !! functions are of the order of epsfcn. if epsfcn is less
                                      !! than the machine precision, it is assumed that the relative
                                      !! errors in the functions are of the order of the machine
                                      !! precision.
        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
                                      !! initial step bound. this bound is set to the product of
                                      !! factor and the euclidean norm of diag*x if nonzero, or else
                                      !! to factor itself. in most cases factor should lie in the
                                      !! interval (.1,100.). 100. is a generally recommended value.
        real(wp), intent(inout) :: x(n) !!  an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
                                          !! below), diag is internally set. if mode = 2, diag
                                          !! must contain positive entries that serve as
                                          !! multiplicative scale factors for the variables.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
                                                !! of fjac contains an upper triangular matrix r with
                                                !! diagonal elements of nonincreasing magnitude such that
                                                !!```
                                                !!        t     t           t
                                                !!       p *(jac *jac)*p = r *r,
                                                !!```
                                                !! where p is a permutation matrix and jac is the final
                                                !! calculated jacobian. column j of p is column ipvt(j)
                                                !! (see below) of the identity matrix. the lower trapezoidal
                                                !! part of fjac contains information generated during
                                                !! the computation of r.
        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
                                       !! the first n elements of the vector (q transpose)*fvec.
        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
        real(wp), intent(inout) :: Wa4(m) !! work array of length n.

        type(lmdif_data) :: wrapper
        wrapper%fcn => fcn

        call lmdif_(wrap_func2, m, n, x, Fvec, Ftol, Xtol, Gtol, Maxfev, Epsfcn, Diag, &
                    Mode, Factor, Nprint, Info, Nfev, Fjac, Ldfjac, Ipvt, &
                    Qtf, Wa1, Wa2, Wa3, Wa4, wrapper)
    end subroutine lmdif
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of lmdif1 is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of the
!  levenberg-marquardt algorithm. this is done by using the more
!  general least-squares solver lmdif. the user must provide a
!  subroutine which calculates the functions. the jacobian is
!  then calculated by a forward-difference approximation.

    subroutine lmdif1(fcn, m, n, x, Fvec, Tol, Info, Iwa, Wa, Lwa)
        implicit none

        procedure(func2) :: fcn !! the user-supplied subroutine which
                                !! calculates the functions.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows:
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  algorithm estimates that the relative error
                                    !!    in the sum of squares is at most tol.
                                    !!  * ***info = 2***  algorithm estimates that the relative error
                                    !!    between x and the solution is at most tol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
                                    !!    jacobian to machine precision.
                                    !!  * ***info = 5***  number of calls to fcn has reached or
                                    !!    exceeded 200*(n+1).
                                    !!  * ***info = 6***  tol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  tol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than
                                  !! m*n+5*n+m.
        integer, intent(inout) :: Iwa(n) !! an integer work array of length n.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                   !! when the algorithm estimates either that the relative
                                   !! error in the sum of squares is at most tol or that
                                   !! the relative error between x and the solution is at
                                   !! most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

        type(lmdif_data) :: wrapper
        wrapper%fcn => fcn

        call lmdif1_(wrap_func2, m, n, x, Fvec, Tol, Info, Iwa, Wa, Lwa, wrapper)
    end subroutine lmdif1
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of lmstr is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm which uses minimal storage.
!  the user must provide a subroutine which calculates the
!  functions and the rows of the jacobian.

    subroutine lmstr(fcn, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
                     Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
                     Wa1, Wa2, Wa3, Wa4)
        implicit none

        procedure(fcn_lmstr) :: fcn !! user-supplied subroutine which
                                    !! calculates the functions and the rows of the jacobian.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
                                     !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
                                     !! occurs when the number of calls to fcn with iflag = 1
                                     !! has reached maxfev.
        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
                                   !! variables will be scaled internally. if mode = 2,
                                   !! the scaling is specified by the input diag. other
                                   !! values of mode are equivalent to mode = 1.
        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
                                     !! printing of iterates if it is positive. in this case,
                                     !! fcn is called with iflag = 0 at the beginning of the first
                                     !! iteration and every nprint iterations thereafter and
                                     !! immediately prior to return, with x and fvec available
                                     !! for printing. if nprint is not positive, no special calls
                                     !! of fcn with iflag = 0 are made.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows:
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  both actual and predicted relative reductions
                                    !!    in the sum of squares are at most ftol.
                                    !!  * ***info = 2***  relative error between two consecutive iterates
                                    !!    is at most xtol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  the cosine of the angle between fvec and any
                                    !!    column of the jacobian is at most gtol in
                                    !!    absolute value.
                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                    !!    reached maxfev.
                                    !!  * ***info = 6***  ftol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  xtol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
                                    !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
                                    !!    columns of the jacobian to machine precision.
        integer, intent(out) :: Nfev !! an integer output variable set to the number of
                                    !! calls to fcn with iflag = 1.
        integer, intent(out) :: Njev !! an integer output variable set to the number of
                                    !! calls to fcn with iflag = 2.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                       !! defines a permutation matrix p such that jac*p = q*r,
                                       !! where jac is the final calculated jacobian, q is
                                       !! orthogonal (not stored), and r is upper triangular.
                                       !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
                                    !! occurs when both the actual and predicted relative
                                    !! reductions in the sum of squares are at most ftol.
                                    !! therefore, ftol measures the relative error desired
                                    !! in the sum of squares.
        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
                                    !! occurs when the relative error between two consecutive
                                    !! iterates is at most xtol. therefore, xtol measures the
                                    !! relative error desired in the approximate solution.
        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
                                    !! occurs when the cosine of the angle between fvec and
                                    !! any column of the jacobian is at most gtol in absolute
                                    !! value. therefore, gtol measures the orthogonality
                                    !! desired between the function vector and the columns
                                    !! of the jacobian.
        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
                                      !! initial step bound. this bound is set to the product of
                                      !! factor and the euclidean norm of diag*x if nonzero, or else
                                      !! to factor itself. in most cases factor should lie in the
                                      !! interval (.1,100.). 100. is a generally recommended value.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array. the upper triangle of fjac
                                                !! contains an upper triangular matrix r such that
                                                !!```
                                                !!        t     t           t
                                                !!       p *(jac *jac)*p = r *r,
                                                !!```
                                                !! where p is a permutation matrix and jac is the final
                                                !! calculated jacobian. column j of p is column ipvt(j)
                                                !! (see below) of the identity matrix. the lower triangular
                                                !! part of fjac contains information generated during
                                                !! the computation of r.
        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
                                          !! below), diag is internally set. if mode = 2, diag
                                          !! must contain positive entries that serve as
                                          !! multiplicative scale factors for the variables.
        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
                                       !! the first n elements of the vector (q transpose)*fvec.
        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
        real(wp), intent(inout) :: Wa4(m) !! work array of length m.

        type(lmstr_data) :: wrapper
        wrapper%fcn => fcn

        call lmstr_(wrap_fcn_lmstr, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
                    Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
                    Wa1, Wa2, Wa3, Wa4, wrapper)
    end subroutine lmstr
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of lmstr1 is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm which uses minimal storage.
!  this is done by using the more general least-squares solver
!  lmstr. the user must provide a subroutine which calculates
!  the functions and the rows of the jacobian.

    subroutine lmstr1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa)
        implicit none

        procedure(fcn_lmstr) :: fcn !! user-supplied subroutine which
                                    !! calculates the functions and the rows of the jacobian.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
                                     !! which specifies the leading dimension of the array fjac.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows:
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  algorithm estimates that the relative error
                                    !!           in the sum of squares is at most tol.
                                    !!  * ***info = 2***  algorithm estimates that the relative error
                                    !!           between x and the solution is at most tol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
                                    !!           jacobian to machine precision.
                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                    !!           reached 100*(n+1).
                                    !!  * ***info = 6***  tol is too small. no further reduction in
                                    !!           the sum of squares is possible.
                                    !!  * ***info = 7***  tol is too small. no further improvement in
                                    !!           the approximate solution x is possible.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than 5*n+m.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                       !! defines a permutation matrix p such that jac*p = q*r,
                                       !! where jac is the final calculated jacobian, q is
                                       !! orthogonal (not stored), and r is upper triangular.
                                       !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                   !! when the algorithm estimates either that the relative
                                   !! error in the sum of squares is at most tol or that
                                   !! the relative error between x and the solution is at
                                   !! most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array. the upper triangle of fjac
                                                !! contains an upper triangular matrix r such that
                                                !!```
                                                !!        t     t           t
                                                !!       p *(jac *jac)*p = r *r,
                                                !!```
                                                !! where p is a permutation matrix and jac is the final
                                                !! calculated jacobian. column j of p is column ipvt(j)
                                                !! (see below) of the identity matrix. the lower triangular
                                                !! part of fjac contains information generated during
                                                !! the computation of r.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

        type(lmstr_data) :: wrapper
        wrapper%fcn => fcn

        call lmstr1_(wrap_fcn_lmstr, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa, &
            & wrapper)
    end subroutine lmstr1
!*****************************************************************************************


    subroutine wrap_func(n, x, fvec, iflag, wrapper)
        !! user-supplied subroutine for [[hybrd]], [[hybrd1]], and [[fdjac1]]
        integer, intent(in) :: n !! the number of variables.
        real(wp), intent(in) :: x(n) !! independent variable vector
        real(wp), intent(out) :: fvec(n) !! value of function at `x`
        integer, intent(inout) :: iflag !! set to <0 to terminate execution
        class(*), intent(inout), optional :: wrapper

        if (.not.present(wrapper)) then
            iflag = -1
            return
        end if

        select type(wrapper)
        class is(hybrd_data)
            call wrapper%fcn(n, x, fvec, iflag)
        class default
            iflag = -1
        end select
    end subroutine wrap_func

    subroutine wrap_func2(m, n, x, fvec, iflag, wrapper)
        !! user-supplied subroutine for [[fdjac2]], [[lmdif]], and [[lmdif1]]
        integer, intent(in) :: m !! the number of functions.
        integer, intent(in) :: n !! the number of variables.
        real(wp), intent(in) :: x(n) !! independent variable vector
        real(wp), intent(out) :: fvec(m) !! value of function at `x`
        integer, intent(inout) :: iflag !! the value of iflag should not be changed unless
                                       !! the user wants to terminate execution of lmdif.
                                       !! in this case set iflag to a negative integer.
        class(*), intent(inout), optional :: wrapper

        if (.not.present(wrapper)) then
            iflag = -1
            return
        end if

        select type(wrapper)
        class is(lmdif_data)
            call wrapper%fcn(m, n, x, fvec, iflag)
        class default
            iflag = -1
        end select
    end subroutine wrap_func2

    subroutine wrap_fcn_hybrj(n, x, fvec, fjac, ldfjac, iflag, wrapper)
        !! user-supplied subroutine for [[hybrj]] and [[hybrj1]]
        integer, intent(in) :: n !! the number of variables.
        real(wp), dimension(n), intent(in) :: x !! independent variable vector
        integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
        real(wp), dimension(n), intent(inout) :: fvec !! value of function at `x`
        real(wp), dimension(ldfjac, n), intent(inout) :: fjac !! jacobian matrix at `x`
        integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                        !! return this vector in fvec. do not alter fjac.
                                        !! if iflag = 2 calculate the jacobian at x and
                                        !! return this matrix in fjac. do not alter fvec.
                                        !!
                                        !! the value of iflag should not be changed by fcn unless
                                        !! the user wants to terminate execution of hybrj.
                                        !! in this case set iflag to a negative integer.
        class(*), intent(inout), optional :: wrapper

        if (.not.present(wrapper)) then
            iflag = -1
            return
        end if

        select type(wrapper)
        class is(hybrj_data)
            call wrapper%fcn(n, x, fvec, fjac, ldfjac, iflag)
        class default
            iflag = -1
        end select
    end subroutine wrap_fcn_hybrj

    subroutine wrap_fcn_lmder(m, n, x, fvec, fjac, ldfjac, iflag, wrapper)
        !! user-supplied subroutine for [[lmder]] and [[lmder1]]
        integer, intent(in) :: m !! the number of functions.
        integer, intent(in) :: n !! the number of variables.
        integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
        integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                       !! return this vector in fvec. do not alter fjac.
                                       !! if iflag = 2 calculate the jacobian at x and
                                       !! return this matrix in fjac. do not alter fvec.
                                       !!
                                       !! the value of iflag should not be changed by fcn unless
                                       !! the user wants to terminate execution of lmder.
                                       !! in this case set iflag to a negative integer.
        real(wp), intent(in) :: x(n) !! independent variable vector
        real(wp), intent(inout) :: fvec(m) !! value of function at `x`
        real(wp), intent(inout) :: fjac(ldfjac, n) !! jacobian matrix at `x`
        class(*), intent(inout), optional :: wrapper

        if (.not.present(wrapper)) then
            iflag = -1
            return
        end if

        select type(wrapper)
        class is(lmder_data)
            call wrapper%fcn(m, n, x, fvec, fjac, ldfjac, iflag)
        class default
            iflag = -1
        end select
    end subroutine wrap_fcn_lmder

    subroutine wrap_fcn_lmstr(m, n, x, fvec, fjrow, iflag, wrapper)
        integer, intent(in) :: m !! the number of functions.
        integer, intent(in) :: n !! the number of variables.
        integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                   !! return this vector in fvec.
                                   !! if iflag = i calculate the (i-1)-st row of the
                                   !! jacobian at x and return this vector in fjrow.
                                   !!
                                   !! the value of iflag should not be changed by fcn unless
                                   !! the user wants to terminate execution of lmstr.
                                   !! in this case set iflag to a negative integer.
        real(wp), intent(in) :: x(n) !! independent variable vector
        real(wp), intent(inout) :: fvec(m) !! value of function at `x`
        real(wp), intent(inout) :: fjrow(n) !! jacobian row
        class(*), intent(inout), optional :: wrapper

        if (.not.present(wrapper)) then
            iflag = -1
            return
        end if

        select type(wrapper)
        class is(lmstr_data)
            call wrapper%fcn(m, n, x, fvec, fjrow, iflag)
        class default
            iflag = -1
        end select
    end subroutine wrap_fcn_lmstr

!*****************************************************************************************
end module minpack_legacy
!*****************************************************************************************
