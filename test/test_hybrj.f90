!*****************************************************************************************
!>
!  This program tests codes for the solution of n nonlinear
!  equations in n variables. it consists of a driver and an
!  interface subroutine fcn. the driver reads in data, calls the
!  nonlinear equation solver, and finally prints out information
!  on the performance of the solver. this is only a sample driver,
!  many other drivers are possible. the interface subroutine fcn
!  is necessary to take into account the forms of calling
!  sequences used by the function and jacobian subroutines in
!  the various nonlinear equation solvers.

program test_hybrj

    use minpack_module
    use iso_fortran_env, only: nwrite => output_unit

    implicit none

    ! originally from file21
    integer,parameter :: ncases = 22
    integer,dimension(ncases),parameter :: nprobs  = [1,2,3,4,5,6,6,7,7,7,7,7,8,8,8,9,10,10,11,12,13,14]
    integer,dimension(ncases),parameter :: ns      = [2,4,2,4,3,6,9,5,6,7,8,9,10,30,40,10,1,10,10,10,10,10]
    integer,dimension(ncases),parameter :: ntriess = [3,3,2,3,3,2,2,3,3,3,1,1,3,1,1,3,3,3,3,3,3,3]

    integer :: i, ic, info, k, n, NFEv, NJEv, NPRob, ntries, icase, lwa, ldfjac
    integer :: na(60), nf(60), nj(60), np(60), nx(60)
    real(wp) :: fnm(60)
    real(wp) :: factor, fnorm1, fnorm2
    real(wp),allocatable :: fjac(:,:), fvec(:), wa(:), x(:)

    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: ten = 10.0_wp
    real(wp), parameter :: tol = sqrt(dpmpar(1))

    ic = 0

    do icase = 1, ncases+1
        if (icase == ncases+1) then
            write (nwrite, '(A,I3,A/)') '1SUMMARY OF ', ic, ' CALLS TO HYBRJ1'
            write (nwrite, '(A/)')      ' NPROB   N    NFEV   NJEV  INFO  FINAL L2 NORM'
            do i = 1, ic
                write (nwrite, '(I4,I6,2I7,I6,1X,D15.7)') np(i), na(i), nf(i), nj(i), nx(i), fnm(i)
            end do
            stop
        else
            nprob = nprobs(icase)
            n = ns(icase)
            ldfjac = n
            lwa = (n*(n+13))/2
            ntries = ntriess(icase)

            if (allocated(fjac)) deallocate(fjac)
            if (allocated(fvec)) deallocate(fvec)
            if (allocated(wa)) deallocate(wa)
            if (allocated(x)) deallocate(x)
            allocate(fjac(n,n))
            allocate(fvec(n))
            allocate(wa(lwa))
            allocate(x(n))

            factor = one
            do k = 1, ntries
                ic = ic + 1
                call initpt(n, x, NPRob, factor)
                call vecfcn(n, x, fvec, NPRob)
                fnorm1 = enorm(n, fvec)
                write (nwrite, '(////5X,A,I5,5X,A,I5,5X//)') ' PROBLEM', NPRob, ' DIMENSION', n
                NFEv = 0
                NJEv = 0
                call hybrj1(fcn, n, x, fvec, fjac, ldfjac, tol, info, wa, lwa)
                fnorm2 = enorm(n, fvec)
                np(ic) = NPRob
                na(ic) = n
                nf(ic) = NFEv
                nj(ic) = NJEv
                nx(ic) = info
                fnm(ic) = fnorm2
                write (nwrite, '(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,I10//5X,A,18X,I10//5X,A//*(5X,5D15.7/))') &
                        ' INITIAL L2 NORM OF THE RESIDUALS', fnorm1, &
                        ' FINAL L2 NORM OF THE RESIDUALS  ', fnorm2, &
                        ' NUMBER OF FUNCTION EVALUATIONS  ', NFEv,   &
                        ' NUMBER OF JACOBIAN EVALUATIONS  ', NJEv,   &
                        ' EXIT PARAMETER', info, &
                        ' FINAL APPROXIMATE SOLUTION', x(1:n)
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
!  and jacobian subroutines vecfcn and vecjac with the
!  appropriate value of problem number (nprob).

    subroutine fcn(n, x, Fvec, Fjac, Ldfjac, Iflag)
        implicit none

        integer,intent(in) :: n
        integer,intent(in) :: Ldfjac
        integer,intent(inout) :: Iflag
        real(wp),intent(in) :: x(n)
        real(wp),intent(inout) :: Fvec(n)
        real(wp),intent(inout) :: Fjac(Ldfjac, n)

        select case (iflag)
        case (1)
            call vecfcn(n, x, Fvec, NPRob)
            NFEv = NFEv + 1
        case(2)
            call vecjac(n, x, Fjac, Ldfjac, NPRob)
            NJEv = NJEv + 1
        case default
            error stop 'invalid iflag value'
        end select

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
!  This subroutine defines the jacobian matrices of fourteen
!  test functions. the problem dimensions are as described
!  in the prologue comments of vecfcn.

    subroutine vecjac(n, x, Fjac, Ldfjac, Nprob)
        implicit none

        integer,intent(in) :: n !! a positive integer variable.
        integer,intent(in) :: Ldfjac !! a positive integer variable not less than n
                                     !! which specifies the leading dimension of the array fjac.
        integer,intent(in) :: Nprob !! a positive integer variable which defines the
                                    !! number of the problem. nprob must not exceed 14.
        real(wp),intent(in) :: x(n) !! an array of length n.
        real(wp),intent(out) :: Fjac(Ldfjac, n) !! an n by n array. on output fjac contains the
                                                !! jacobian matrix of the nprob function evaluated at x.

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: two = 2.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: four = 4.0_wp
        real(wp), parameter :: five = 5.0_wp
        real(wp), parameter :: six = 6.0_wp
        real(wp), parameter :: eight = 8.0_wp
        real(wp), parameter :: ten = 10.0_wp
        real(wp), parameter :: fiftn = 15.0_wp
        real(wp), parameter :: twenty = 20.0_wp
        real(wp), parameter :: hundrd = 100.0_wp
        real(wp), parameter :: c1 = 1.0e4_wp
        real(wp), parameter :: c3 = 2.0e2_wp
        real(wp), parameter :: c4 = 2.02e1_wp
        real(wp), parameter :: c5 = 1.98e1_wp
        real(wp), parameter :: c6 = 1.8e2_wp
        real(wp), parameter :: c9 = 2.9e1_wp

        integer :: i, j, k, k1, k2, ml, mu
        real(wp) :: h, prod, sum, sum1, sum2, temp, temp1, temp2, &
                    temp3, temp4, ti, tj, tk, tpi

        Fjac(1:n,1:n) = zero

        ! jacobian routine selector.

        select case (nprob)
        case (2)
            ! powell singular function.
            do k = 1, 4
                do j = 1, 4
                    fjac(k, j) = zero
                end do
            end do
            fjac(1, 1) = one
            fjac(1, 2) = ten
            fjac(2, 3) = sqrt(five)
            fjac(2, 4) = -fjac(2, 3)
            fjac(3, 2) = two*(x(2) - two*x(3))
            fjac(3, 3) = -two*fjac(3, 2)
            fjac(4, 1) = two*sqrt(ten)*(x(1) - x(4))
            fjac(4, 4) = -fjac(4, 1)
        case (3)
            ! powell badly scaled function.
            fjac(1, 1) = c1*x(2)
            fjac(1, 2) = c1*x(1)
            fjac(2, 1) = -exp(-x(1))
            fjac(2, 2) = -exp(-x(2))
        case (4)
            ! wood function.
            do k = 1, 4
                do j = 1, 4
                    fjac(k, j) = zero
                end do
            end do
            temp1 = x(2) - three*x(1)**2
            temp2 = x(4) - three*x(3)**2
            fjac(1, 1) = -c3*temp1 + one
            fjac(1, 2) = -c3*x(1)
            fjac(2, 1) = -two*c3*x(1)
            fjac(2, 2) = c3 + c4
            fjac(2, 4) = c5
            fjac(3, 3) = -c6*temp2 + one
            fjac(3, 4) = -c6*x(3)
            fjac(4, 2) = c5
            fjac(4, 3) = -two*c6*x(3)
            fjac(4, 4) = c6 + c4
        case (5)
            ! helical valley function.
            tpi = eight*atan(one)
            temp = x(1)**2 + x(2)**2
            temp1 = tpi*temp
            temp2 = sqrt(temp)
            fjac(1, 1) = hundrd*x(2)/temp1
            fjac(1, 2) = -hundrd*x(1)/temp1
            fjac(1, 3) = ten
            fjac(2, 1) = ten*x(1)/temp2
            fjac(2, 2) = ten*x(2)/temp2
            fjac(2, 3) = zero
            fjac(3, 1) = zero
            fjac(3, 2) = zero
            fjac(3, 3) = one
        case (6)
            ! watson function.
            do k = 1, n
                do j = k, n
                    fjac(k, j) = zero
                end do
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
                temp1 = two*(sum1 - sum2**2 - one)
                temp2 = two*sum2
                temp = ti**2
                tk = one
                do k = 1, n
                    tj = tk
                    do j = k, n
                        fjac(k, j) = fjac(k, j) &
                                     + tj*((dfloat(k - 1)/ti - temp2)*(dfloat(j - 1) &
                                     /ti - temp2) - temp1)
                        tj = ti*tj
                    end do
                    tk = temp*tk
                end do
            end do
            fjac(1, 1) = fjac(1, 1) + six*x(1)**2 - two*x(2) + three
            fjac(1, 2) = fjac(1, 2) - two*x(1)
            fjac(2, 2) = fjac(2, 2) + one
            do k = 1, n
                do j = k, n
                    fjac(j, k) = fjac(k, j)
                end do
            end do
        case (7)
            ! chebyquad function.
            tk = one/dfloat(n)
            do j = 1, n
                temp1 = one
                temp2 = two*x(j) - one
                temp = two*temp2
                temp3 = zero
                temp4 = two
                do k = 1, n
                    fjac(k, j) = tk*temp4
                    ti = four*temp2 + temp*temp4 - temp3
                    temp3 = temp4
                    temp4 = ti
                    ti = temp*temp2 - temp1
                    temp1 = temp2
                    temp2 = ti
                end do
            end do
        case (8)
            ! brown almost-linear function.
            prod = one
            do j = 1, n
                prod = x(j)*prod
                do k = 1, n
                    fjac(k, j) = one
                end do
                fjac(j, j) = two
            end do
            do j = 1, n
                temp = x(j)
                if (temp == zero) then
                    temp = one
                    prod = one
                    do k = 1, n
                        if (k /= j) prod = x(k)*prod
                    end do
                end if
                fjac(n, j) = prod/temp
            end do
        case (9)
            ! discrete boundary value function.
            h = one/dfloat(n + 1)
            do k = 1, n
                temp = three*(x(k) + dfloat(k)*h + one)**2
                do j = 1, n
                    fjac(k, j) = zero
                end do
                fjac(k, k) = two + temp*h**2/two
                if (k /= 1) fjac(k, k - 1) = -one
                if (k /= n) fjac(k, k + 1) = -one
            end do
        case (10)
            ! discrete integral equation function.
            h = one/dfloat(n + 1)
            do k = 1, n
                tk = dfloat(k)*h
                do j = 1, n
                    tj = dfloat(j)*h
                    temp = three*(x(j) + tj + one)**2
                    fjac(k, j) = h*min(tj*(one - tk), tk*(one - tj))*temp/two
                end do
                fjac(k, k) = fjac(k, k) + one
            end do
        case (11)
            ! trigonometric function.
            do j = 1, n
                temp = sin(x(j))
                do k = 1, n
                    fjac(k, j) = temp
                end do
                fjac(j, j) = dfloat(j + 1)*temp - cos(x(j))
            end do
        case (12)
            ! variably dimensioned function.
            sum = zero
            do j = 1, n
                sum = sum + dfloat(j)*(x(j) - one)
            end do
            temp = one + six*sum**2
            do k = 1, n
                do j = k, n
                    fjac(k, j) = dfloat(k*j)*temp
                    fjac(j, k) = fjac(k, j)
                end do
                fjac(k, k) = fjac(k, k) + one
            end do
        case (13)
            ! broyden tridiagonal function.
            do k = 1, n
                do j = 1, n
                    fjac(k, j) = zero
                end do
                fjac(k, k) = three - four*x(k)
                if (k /= 1) fjac(k, k - 1) = -one
                if (k /= n) fjac(k, k + 1) = -two
            end do
        case (14)
            ! broyden banded function.
            ml = 5
            mu = 1
            do k = 1, n
                do j = 1, n
                    fjac(k, j) = zero
                end do
                k1 = max(1, k - ml)
                k2 = min(k + mu, n)
                do j = k1, k2
                    if (j /= k) fjac(k, j) = -(one + two*x(j))
                end do
                fjac(k, k) = two + fiftn*x(k)**2
            end do
        case default
            ! rosenbrock function.
            fjac(1, 1) = -one
            fjac(1, 2) = zero
            fjac(2, 1) = -twenty*x(1)
            fjac(2, 2) = ten
        end select

    end subroutine vecjac
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

        integer,intent(in) :: n !! a positive integer input variable.
        integer,intent(in) :: Nprob !! a positive integer input variable which defines the
                                    !! number of the problem. nprob must not exceed 14.
        real(wp),intent(in) :: Factor !! an input variable which specifies the multiple of
                                      !! the standard starting point. if factor is unity, no
                                      !! multiplication is performed.
        real(wp),intent(out) :: x(n) !! an output array of length n which contains the standard
                                     !! starting point for problem nprob multiplied by factor.

        integer :: j
        real(wp) :: h, tj

        real(wp),parameter :: zero = 0.0_wp
        real(wp),parameter :: half = 0.5_wp
        real(wp),parameter :: one = 1.0_wp
        real(wp),parameter :: three = 3.0_wp
        real(wp),parameter :: c1 = 1.2_wp

        x(1:n) = zero

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
!>
!  This subroutine defines fourteen test functions. the first
!  five test functions are of dimensions 2,4,2,4,3, respectively,
!  while the remaining test functions are of variable dimension
!  n for any n greater than or equal to 1 (problem 6 is an
!  exception to this, since it does not allow n = 1).

    subroutine vecfcn(n, x, Fvec, Nprob)
        implicit none

        integer,intent(in) :: n !! a positive integer input variable.
        real(wp),intent(in) :: x(n) !! an input array of length n.
        real(wp),intent(out) :: Fvec(n) !! an output array of length n which contains the nprob
                                        !! function vector evaluated at x.
        integer,intent(in) :: Nprob !! a positive integer input variable which defines the
                                    !! number of the problem. nprob must not exceed 14.

        real(wp),parameter :: zero = 0.0_wp
        real(wp),parameter :: one = 1.0_wp
        real(wp),parameter :: two = 2.0_wp
        real(wp),parameter :: three = 3.0_wp
        real(wp),parameter :: five = 5.0_wp
        real(wp),parameter :: eight = 8.0_wp
        real(wp),parameter :: ten = 10.0_wp

        real(wp),parameter :: c1 =  1.0e4_wp
        real(wp),parameter :: c2 =  1.0001_wp
        real(wp),parameter :: c3 =  2.0e2_wp
        real(wp),parameter :: c4 =  2.02e1_wp
        real(wp),parameter :: c5 =  1.98e1_wp
        real(wp),parameter :: c6 =  1.8e2_wp
        real(wp),parameter :: c7 =  2.5e-1_wp
        real(wp),parameter :: c8 =  5.0e-1_wp
        real(wp),parameter :: c9 =  2.9e1_wp

        real(wp),parameter :: tpi = eight*atan(one)

        integer :: i, iev, j, k, k1, k2, kp1, ml, mu
        real(wp) :: h, prod, sum, sum1, sum2, temp, temp1, &
                    temp2, ti, tj, tk

        Fvec(1:n) = zero

        ! problem selector.

        select case (nprob)
        case (2)
            ! powell singular function.
            fvec(1) = x(1) + ten*x(2)
            fvec(2) = sqrt(five)*(x(3) - x(4))
            fvec(3) = (x(2) - two*x(3))**2
            fvec(4) = sqrt(ten)*(x(1) - x(4))**2
        case (3)
            ! powell badly scaled function.
            fvec(1) = c1*x(1)*x(2) - one
            fvec(2) = exp(-x(1)) + exp(-x(2)) - c2
        case (4)
            ! wood function.
            temp1 = x(2) - x(1)**2
            temp2 = x(4) - x(3)**2
            fvec(1) = -c3*x(1)*temp1 - (one - x(1))
            fvec(2) = c3*temp1 + c4*(x(2) - one) + c5*(x(4) - one)
            fvec(3) = -c6*x(3)*temp2 - (one - x(3))
            fvec(4) = c6*temp2 + c4*(x(4) - one) + c5*(x(2) - one)
        case (5)
            ! helical valley function.
            temp1 = sign(c7, x(2))
            if (x(1) > zero) temp1 = atan(x(2)/x(1))/tpi
            if (x(1) < zero) temp1 = atan(x(2)/x(1))/tpi + c8
            temp2 = sqrt(x(1)**2 + x(2)**2)
            fvec(1) = ten*(x(3) - ten*temp1)
            fvec(2) = ten*(temp2 - one)
            fvec(3) = x(3)
        case (6)
            ! watson function.
            do k = 1, n
                fvec(k) = zero
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
                    fvec(k) = fvec(k) + temp*(dfloat(k - 1) - temp2)*temp1
                    temp = ti*temp
                end do
            end do
            temp = x(2) - x(1)**2 - one
            fvec(1) = fvec(1) + x(1)*(one - two*temp)
            fvec(2) = fvec(2) + temp
        case (7)
            ! chebyquad function.
            do k = 1, n
                fvec(k) = zero
            end do
            do j = 1, n
                temp1 = one
                temp2 = two*x(j) - one
                temp = two*temp2
                do i = 1, n
                    fvec(i) = fvec(i) + temp2
                    ti = temp*temp2 - temp1
                    temp1 = temp2
                    temp2 = ti
                end do
            end do
            tk = one/dfloat(n)
            iev = -1
            do k = 1, n
                fvec(k) = tk*fvec(k)
                if (iev > 0) fvec(k) = fvec(k) + one/(dfloat(k)**2 - one)
                iev = -iev
            end do
        case (8)
            ! brown almost-linear function.
            sum = -dfloat(n + 1)
            prod = one
            do j = 1, n
                sum = sum + x(j)
                prod = x(j)*prod
            end do
            do k = 1, n
                fvec(k) = x(k) + sum
            end do
            fvec(n) = prod - one
        case (9)
            ! discrete boundary value function.
            h = one/dfloat(n + 1)
            do k = 1, n
                temp = (x(k) + dfloat(k)*h + one)**3
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                fvec(k) = two*x(k) - temp1 - temp2 + temp*h**2/two
            end do
        case (10)
            ! discrete integral equation function.
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
                fvec(k) = x(k) + h*((one - tk)*sum1 + tk*sum2)/two
            end do
        case (11)
            ! trigonometric function.
            sum = zero
            do j = 1, n
                fvec(j) = cos(x(j))
                sum = sum + fvec(j)
            end do
            do k = 1, n
                fvec(k) = dfloat(n + k) - sin(x(k)) - sum - dfloat(k)*fvec(k)
            end do
        case (12)
            ! variably dimensioned function.
            sum = zero
            do j = 1, n
                sum = sum + dfloat(j)*(x(j) - one)
            end do
            temp = sum*(one + two*sum**2)
            do k = 1, n
                fvec(k) = x(k) - one + dfloat(k)*temp
            end do
        case (13)
            ! broyden tridiagonal function.
            do k = 1, n
                temp = (three - two*x(k))*x(k)
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                fvec(k) = temp - temp1 - two*temp2 + one
            end do
        case (14)
            ! broyden banded function.
            ml = 5
            mu = 1
            do k = 1, n
                k1 = max(1, k - ml)
                k2 = min(k + mu, n)
                temp = zero
                do j = k1, k2
                    if (j /= k) temp = temp + x(j)*(one + x(j))
                end do
                fvec(k) = x(k)*(two + five*x(k)**2) + one - temp
            end do
        case default
            ! rosenbrock function.
            fvec(1) = one - x(1)
            fvec(2) = ten*(x(2) - x(1)**2)
        end select

    end subroutine vecfcn

end program test_hybrj
