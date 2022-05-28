!> The problem is to determine the values of x(1), x(2), ..., x(9)
!>  which solve the system of tridiagonal equations
!>     (3-2*x(1))*x(1)           -2*x(2)                   = -1
!>             -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!>                                 -x(8) + (3-2*x(9))*x(9) = -1
program example_hybrd

    use minpack_legacy, only: hybrd
    use minpack_module, only: wp, enorm, dpmpar
    use iso_fortran_env, only: nwrite => output_unit

    implicit none

    integer,parameter :: n = 9
    integer,parameter :: ldfjac = n
    integer,parameter :: lr = (n*(n+1))/2

    integer :: maxfev, ml, mu, mode, nprint, info, nfev
    real(wp) :: epsfcn, factor, fnorm, xtol
    real(wp) :: x(n), fvec(n), diag(n), fjac(n, n), r(lr), qtf(n), &
                wa1(n), wa2(n), wa3(n), wa4(n)

    xtol = sqrt(dpmpar(1))  ! square root of the machine precision.
    maxfev = 2000
    ml = 1
    mu = 1
    epsfcn = 0.0_wp
    mode = 2
    factor = 100.0_wp
    nprint = 0
    diag = 1.0_wp
    x = -1.0_wp  !  starting values to provide a rough solution.

    call hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, &
               mode, factor, nprint, info, nfev, fjac, ldfjac, &
               r, lr, qtf, wa1, wa2, wa3, wa4)
    fnorm = enorm(n, fvec)

    write (nwrite, '(5x,a,d15.7//5x,a,i10//5x,a,16x,i10//5x,a//(5x,3d15.7))') &
           "FINAL L2 NORM OF THE RESIDUALS", fnorm, &
           "NUMBER OF FUNCTION EVALUATIONS", nfev, &
           "EXIT PARAMETER", info, &
           "FINAL APPROXIMATE SOLUTION", x

    !> Results obtained with different compilers or machines
    !>  may be slightly different.
    !>
    !>> FINAL L2 NORM OF THE RESIDUALS  0.1192636D-07
    !>>
    !>> NUMBER OF FUNCTION EVALUATIONS        14
    !>>
    !>> EXIT PARAMETER                         1
    !>>
    !>> FINAL APPROXIMATE SOLUTION
    !>>
    !>>  -0.5706545D+00 -0.6816283D+00 -0.7017325D+00
    !>>  -0.7042129D+00 -0.7013690D+00 -0.6918656D+00
    !>>  -0.6657920D+00 -0.5960342D+00 -0.4164121D+00

contains

    !> Subroutine fcn for hybrd example.
    subroutine fcn(n, x, fvec, iflag)

        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: iflag
        real(wp), intent(in) :: x(n)
        real(wp), intent(out) :: fvec(n)

        integer :: k !! counter
        real(wp) :: temp, temp1, temp2

        real(wp),parameter :: zero = 0.0_wp
        real(wp),parameter :: one = 1.0_wp
        real(wp),parameter :: two = 2.0_wp
        real(wp),parameter :: three = 3.0_wp

        if (iflag == 0) then
            !! Insert print statements here when nprint is positive.
        else
            do k = 1, n
                temp = (three - two*x(k))*x(k)
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                fvec(k) = temp - temp1 - two*temp2 + one
            end do
        end if

    end subroutine fcn

end program example_hybrd
