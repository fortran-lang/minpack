
!*****************************************************************************************
!>
!  This program tests codes for the solution of n nonlinear
!  equations in n variables. it consists of a driver and an
!  interface subroutine fcn. the driver reads in data, calls the
!  nonlinear equation solver, and finally prints out information
!  on the performance of the solver. this is only a sample driver,
!  many other drivers are possible. the interface subroutine fcn
!  is necessary to take into account the forms of calling
!  sequences used by the function subroutines in the various
!  nonlinear equation solvers.
!
!### Author
!  * Argonne national laboratory. minpack project. march 1980.
!    burton s. garbow, kenneth e. hillstrom, jorge j. more

program test

    use minpack_module
    use iso_fortran_env, only: output_unit

    implicit none

    integer :: i, ic, info, k, lwa, n, NFEv, NPRob, &
               ntries, icase
    integer :: na(60), nf(60), np(60), nx(60)
    real(wp) :: factor, fnorm1, fnorm2, tol
    real(wp) :: fnm(60), fvec(40), wa(2660), x(40)

    integer, parameter :: nwrite = output_unit ! logical output unit
    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: ten = 10.0_wp

    tol = sqrt(dpmpar(1))
    lwa = 2660
    ic = 0
    n = 5
    ntries = 1
    do NPRob = 1, 16
        if (NPRob == 16) then
            write (nwrite, 99002) ic
99002       format('1SUMMARY OF ', i3, ' CALLS TO HYBRD1'/)
            write (nwrite, 99003)
99003       format(' NPROB   N    NFEV  INFO  FINAL L2 NORM'/)
            do i = 1, ic
                write (nwrite, 99004) np(i), na(i), nf(i), nx(i), fnm(i)
99004           format(i4, i6, i7, i6, 1x, d15.7)
            end do
            stop
        else
            factor = one
            do k = 1, ntries
                ic = ic + 1
                call initpt(n, x, NPRob, factor)
                call vecfcn(n, x, fvec, NPRob)
                fnorm1 = enorm(n, fvec)
                write (nwrite, 99005) NPRob, n
99005           format(////5x, ' PROBLEM', i5, 5x, ' DIMENSION', i5, 5x//)
                NFEv = 0
                call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
                fnorm2 = enorm(n, fvec)
                np(ic) = NPRob
                na(ic) = n
                nf(ic) = NFEv
                nx(ic) = info
                fnm(ic) = fnorm2
                write (nwrite, 99006) fnorm1, fnorm2, NFEv, info,        &
                                   & (x(i), i=1, n)
99006           format(5x, ' INITIAL L2 NORM OF THE RESIDUALS', d15.7//5x,   &
                                 &' FINAL L2 NORM OF THE RESIDUALS  ', d15.7//5x,      &
                                 &' NUMBER OF FUNCTION EVALUATIONS  ', i10//5x,        &
                                 &' EXIT PARAMETER', 18x, i10//5x,                      &
                                 &' FINAL APPROXIMATE SOLUTION'//(5x, 5d15.7))
                factor = ten*factor
            end do
        end if
    end do

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  The calling sequence of fcn should be identical to the
!  calling sequence of the function subroutine in the nonlinear
!  equation solver. fcn should only call the testing function
!  subroutine vecfcn with the appropriate value of problem
!  number (nprob).

    subroutine fcn(n, x, Fvec, Iflag)
        implicit none

        integer, intent(in) :: n !! the number of variables.
        real(wp), intent(in) :: x(n) !! independant variable vector
        real(wp), intent(out) :: fvec(n) !! value of function at `x`
        integer, intent(inout) :: iflag !! set to <0 to terminate execution

        call vecfcn(n, x, Fvec, NPRob)
        NFEv = NFEv + 1

    end subroutine fcn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replaced statement function in original code.

    pure elemental function dfloat(i) result(f)
        implicit none
        integer, intent(in) :: i
        real(wp) :: f
        f = real(i, wp)
    end function dfloat
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine defines fourteen test functions. the first
!  five test functions are of dimensions 2,4,2,4,3, respectively,
!  while the remaining test functions are of variable dimension
!  n for any n greater than or equal to 1 (problem 6 is an
!  exception to this, since it does not allow n = 1).

    subroutine vecfcn(n, x, Fvec, Nprob)
        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        integer, intent(in) :: nprob !! a positive integer input variable which defines the
                                     !! number of the problem. nprob must not exceed 14.
        real(wp), intent(in) :: x(n) !! an input array of length n.
        real(wp), intent(out) :: fvec(n) !! an output array of length n which contains the nprob
                                         !! function vector evaluated at x.

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: two = 2.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: five = 5.0_wp
        real(wp), parameter :: eight = 8.0_wp
        real(wp), parameter :: ten = 10.0_wp
        real(wp), parameter :: c1 = 1.0e4_wp
        real(wp), parameter :: c2 = 1.0001_wp
        real(wp), parameter :: c3 = 2.0e2_wp
        real(wp), parameter :: c4 = 2.02e1_wp
        real(wp), parameter :: c5 = 1.98e1_wp
        real(wp), parameter :: c6 = 1.8e2_wp
        real(wp), parameter :: c7 = 2.5e-1_wp
        real(wp), parameter :: c8 = 5.0e-1_wp
        real(wp), parameter :: c9 = 2.9e1_wp

        integer :: i, iev, ivar, j, k, k1, k2, kp1, ml, mu
        real(wp) :: h, prod, sum, sum1, sum2, temp, temp1, &
                    temp2, ti, tj, tk, tpi

        ! PROBLEM SELECTOR.

        select case (Nprob)
        case (2)
            ! POWELL SINGULAR FUNCTION.
            Fvec(1) = x(1) + ten*x(2)
            Fvec(2) = sqrt(five)*(x(3) - x(4))
            Fvec(3) = (x(2) - two*x(3))**2
            Fvec(4) = sqrt(ten)*(x(1) - x(4))**2
        case (3)
            ! POWELL BADLY SCALED FUNCTION.
            Fvec(1) = c1*x(1)*x(2) - one
            Fvec(2) = exp(-x(1)) + exp(-x(2)) - c2
        case (4)
            ! WOOD FUNCTION.
            temp1 = x(2) - x(1)**2
            temp2 = x(4) - x(3)**2
            Fvec(1) = -c3*x(1)*temp1 - (one - x(1))
            Fvec(2) = c3*temp1 + c4*(x(2) - one) + c5*(x(4) - one)
            Fvec(3) = -c6*x(3)*temp2 - (one - x(3))
            Fvec(4) = c6*temp2 + c4*(x(4) - one) + c5*(x(2) - one)
        case (5)
            ! HELICAL VALLEY FUNCTION.
            tpi = eight*atan(one)
            temp1 = sign(c7, x(2))
            if (x(1) > zero) temp1 = atan(x(2)/x(1))/tpi
            if (x(1) < zero) temp1 = atan(x(2)/x(1))/tpi + c8
            temp2 = sqrt(x(1)**2 + x(2)**2)
            Fvec(1) = ten*(x(3) - ten*temp1)
            Fvec(2) = ten*(temp2 - one)
            Fvec(3) = x(3)
        case (6)
            ! WATSON FUNCTION.
            do k = 1, n
                Fvec(k) = zero
            end do
            do i = 1, 29
                ti = dfloat(i)/c9
                sum1 = zero
                temp = one
                do j = 2, n
                    sum1 = sum1 + dfloat(j - 1)*temp*x(j)
                    temp = ti*temp
                end do
                sum2 = zero
                temp = one
                do j = 1, n
                    sum2 = sum2 + temp*x(j)
                    temp = ti*temp
                end do
                temp1 = sum1 - sum2**2 - one
                temp2 = two*ti*sum2
                temp = one/ti
                do k = 1, n
                    Fvec(k) = Fvec(k) + temp*(dfloat(k - 1) - temp2)*temp1
                    temp = ti*temp
                end do
            end do
            temp = x(2) - x(1)**2 - one
            Fvec(1) = Fvec(1) + x(1)*(one - two*temp)
            Fvec(2) = Fvec(2) + temp
        case (7)
            ! CHEBYQUAD FUNCTION.
            do k = 1, n
                Fvec(k) = zero
            end do
            do j = 1, n
                temp1 = one
                temp2 = two*x(j) - one
                temp = two*temp2
                do i = 1, n
                    Fvec(i) = Fvec(i) + temp2
                    ti = temp*temp2 - temp1
                    temp1 = temp2
                    temp2 = ti
                end do
            end do
            tk = one/dfloat(n)
            iev = -1
            do k = 1, n
                Fvec(k) = tk*Fvec(k)
                if (iev > 0) Fvec(k) = Fvec(k) + one/(dfloat(k)**2 - one)
                iev = -iev
            end do
        case (8)
            ! BROWN ALMOST-LINEAR FUNCTION.
            sum = -dfloat(n + 1)
            prod = one
            do j = 1, n
                sum = sum + x(j)
                prod = x(j)*prod
            end do
            do k = 1, n
                Fvec(k) = x(k) + sum
            end do
            Fvec(n) = prod - one
        case (9)
            ! DISCRETE BOUNDARY VALUE FUNCTION.
            h = one/dfloat(n + 1)
            do k = 1, n
                temp = (x(k) + dfloat(k)*h + one)**3
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                Fvec(k) = two*x(k) - temp1 - temp2 + temp*h**2/two
            end do
        case (10)
            ! DISCRETE INTEGRAL EQUATION FUNCTION.
            h = one/dfloat(n + 1)
            do k = 1, n
                tk = dfloat(k)*h
                sum1 = zero
                do j = 1, k
                    tj = dfloat(j)*h
                    temp = (x(j) + tj + one)**3
                    sum1 = sum1 + tj*temp
                end do
                sum2 = zero
                kp1 = k + 1
                if (n >= kp1) then
                    do j = kp1, n
                        tj = dfloat(j)*h
                        temp = (x(j) + tj + one)**3
                        sum2 = sum2 + (one - tj)*temp
                    end do
                end if
                Fvec(k) = x(k) + h*((one - tk)*sum1 + tk*sum2)/two
            end do
        case (11)
            ! TRIGONOMETRIC FUNCTION.
            sum = zero
            do j = 1, n
                Fvec(j) = cos(x(j))
                sum = sum + Fvec(j)
            end do
            do k = 1, n
                Fvec(k) = dfloat(n + k) - sin(x(k)) - sum - dfloat(k)*Fvec(k)
            end do
        case (12)
            ! VARIABLY DIMENSIONED FUNCTION.
            sum = zero
            do j = 1, n
                sum = sum + dfloat(j)*(x(j) - one)
            end do
            temp = sum*(one + two*sum**2)
            do k = 1, n
                Fvec(k) = x(k) - one + dfloat(k)*temp
            end do
        case (13)
            ! BROYDEN TRIDIAGONAL FUNCTION.
            do k = 1, n
                temp = (three - two*x(k))*x(k)
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                Fvec(k) = temp - temp1 - two*temp2 + one
            end do
        case (14)
            ! BROYDEN BANDED FUNCTION.
            ml = 5
            mu = 1
            do k = 1, n
                k1 = max(1, k - ml)
                k2 = min(k + mu, n)
                temp = zero
                do j = k1, k2
                    if (j /= k) temp = temp + x(j)*(one + x(j))
                end do
                Fvec(k) = x(k)*(two + five*x(k)**2) + one - temp
            end do
        case default
            ! ROSENBROCK FUNCTION.
            Fvec(1) = one - x(1)
            Fvec(2) = ten*(x(2) - x(1)**2)
        end select

    end subroutine vecfcn
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine specifies the standard starting points for
!  the functions defined by subroutine vecfcn. the subroutine
!  returns in x a multiple (factor) of the standard starting
!  point. for the sixth function the standard starting point is
!  zero, so in this case, if factor is not unity, then the
!  subroutine returns the vector  x(j) = factor, j=1,...,n.

    subroutine initpt(n, x, Nprob, Factor)

        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        real(wp), intent(out) :: x(n) !! an output array of length n which contains the standard
                                      !! starting point for problem nprob multiplied by factor.
        integer, intent(in) :: Nprob !! a positive integer input variable which defines the
                                     !! number of the problem. nprob must not exceed 14.
        real(wp), intent(in) :: Factor !! an input variable which specifies the multiple of
                                       !! the standard starting point. if factor is unity, no
                                       !! multiplication is performed.

        integer :: ivar, j
        real(wp) :: h, tj

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: half = 0.5_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: c1 = 1.2_wp

        ! selection of initial point.

        select case (nprob)
        case (2)
            ! powell singular function.
            x(1) = three
            x(2) = -one
            x(3) = zero
            x(4) = one
        case (3)
            ! powell badly scaled function.
            x(1) = zero
            x(2) = one
        case (4)
            ! wood function.
            x(1) = -three
            x(2) = -one
            x(3) = -three
            x(4) = -one
        case (5)
            ! helical valley function.
            x(1) = -one
            x(2) = zero
            x(3) = zero
        case (6)
            ! watson function.
            do j = 1, n
                x(j) = zero
            end do
        case (7)
            ! chebyquad function.
            h = one/dfloat(n + 1)
            do j = 1, n
                x(j) = dfloat(j)*h
            end do
        case (8)
            ! brown almost-linear function.
            do j = 1, n
                x(j) = half
            end do
        case (9, 10)
            ! discrete boundary value and integral equation functions.
            h = one/dfloat(n + 1)
            do j = 1, n
                tj = dfloat(j)*h
                x(j) = tj*(tj - one)
            end do
        case (11)
            ! trigonometric function.
            h = one/dfloat(n)
            do j = 1, n
                x(j) = h
            end do
        case (12)
            ! variably dimensioned function.
            h = one/dfloat(n)
            do j = 1, n
                x(j) = one - dfloat(j)*h
            end do
        case (13, 14)
            ! broyden tridiagonal and banded functions.
            do j = 1, n
                x(j) = -one
            end do
        case default
            ! rosenbrock function.
            x(1) = -c1
            x(2) = one
        end select

        ! compute multiple of initial point.

        if (factor /= one) then
            if (nprob == 6) then
                do j = 1, n
                    x(j) = factor
                end do
            else
                do j = 1, n
                    x(j) = factor*x(j)
                end do
            end if
        end if

    end subroutine initpt
!*****************************************************************************************

!*****************************************************************************************
end program test
!*****************************************************************************************