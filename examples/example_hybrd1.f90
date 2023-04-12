!> the problem is to determine the values of x(1), x(2), ..., x(9),
!>  which solve the system of tridiagonal equations.
!>
!>  (3-2*x(1))*x(1)           -2*x(2)                   = -1
!>          -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!>                              -x(8) + (3-2*x(9))*x(9) = -1
program example_hybrd1

    use minpack_module, only: wp, hybrd1, dpmpar, enorm
    use iso_fortran_env, only: nwrite => output_unit

    implicit none

    integer,parameter :: n = 9
    integer,parameter :: lwa = (n*(3*n+13))/2

    integer :: info
    real(wp) :: tol, fnorm
    real(wp) :: x(n), fvec(n), wa(lwa)

    !> The following starting values provide a rough solution.
    x = -1.0_wp

    !> Set tol to the square root of the machine precision.
    !>  unless high precision solutions are required,
    !>  this is the recommended setting.
    tol = sqrt(dpmpar(1))

    call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
    fnorm = enorm(n, fvec)

    write (nwrite, '(5x,a,d15.7//5x,a,16x,i10//5x,a//(5x,3d15.7))') &
            "FINAL L2 NORM OF THE RESIDUALS", fnorm, &
            "EXIT PARAMETER", info, &
            "FINAL APPROXIMATE SOLUTION", x

    !> Results obtained with different compilers or machines
    !>  may be slightly different.
    !>
    !>> FINAL L2 NORM OF THE RESIDUALS  0.1192636D-07
    !>>
    !>> EXIT PARAMETER                         1
    !>>
    !>> FINAL APPROXIMATE SOLUTION
    !>>
    !>>  -0.5706545D+00 -0.6816283D+00 -0.7017325D+00
    !>>  -0.7042129D+00 -0.7013690D+00 -0.6918656D+00
    !>>  -0.6657920D+00 -0.5960342D+00 -0.4164121D+00

contains

    !> Subroutine fcn for hybrd1 example.
    subroutine fcn(n, x, fvec, iflag)

        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: iflag
        real(wp), intent(in) :: x(n)
        real(wp), intent(out) :: fvec(n)

        integer :: k
        real(wp) :: temp, temp1, temp2

        real(wp),parameter :: zero = 0.0_wp
        real(wp),parameter :: one = 1.0_wp
        real(wp),parameter :: two = 2.0_wp
        real(wp),parameter :: three = 3.0_wp

        do k = 1, n
            temp = (three - two*x(k))*x(k)
            temp1 = zero
            if (k /= 1) temp1 = x(k - 1)
            temp2 = zero
            if (k /= n) temp2 = x(k + 1)
            fvec(k) = temp - temp1 - two*temp2 + one
        end do

    end subroutine fcn

end program example_hybrd1
