program example_lmdif1

use minpack_legacy, only: lmdif1
use minpack_module, only: wp, enorm
use iso_fortran_env, only: nwrite => output_unit

implicit none

integer, parameter :: n = 3
integer, parameter :: m = 15
integer, parameter :: lwa = m*n+5*n+m

integer :: info
real(wp) :: tol, x(n), fvec(m), wa(lwa)
integer :: iwa(n)

! The following starting values provide a rough fit.
x = [1.0_wp, 1.0_wp, 1.0_wp]

! Set tol to the square root of the machine precision. Unless high precision
! solutions are required, this is the recommended setting.
tol = sqrt(epsilon(1.0_wp))

call lmdif1(fcn, m, n, x, fvec, tol, info, iwa, wa, lwa)

write(nwrite, '(5x,a,d15.7//,5x,a,16x,i10//,5x,a//(5x,3d15.7))') &
        'FINAL L2 NORM OF THE RESIDUALS', enorm(m, fvec), &
        'EXIT PARAMETER', info, &
        'FINAL APPROXIMATE SOLUTION', x

contains

subroutine fcn(m, n, x, fvec, iflag)

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n)
    real(wp), intent(out) :: fvec(m)
    integer, intent(inout) :: iflag

    integer :: i
    real(wp) :: tmp1, tmp2, tmp3

    real(wp),parameter :: y(15) = [1.4e-1_wp, 1.8e-1_wp, 2.2e-1_wp, 2.5e-1_wp, 2.9e-1_wp, &
                                   3.2e-1_wp, 3.5e-1_wp, 3.9e-1_wp, 3.7e-1_wp, 5.8e-1_wp, &
                                   7.3e-1_wp, 9.6e-1_wp, 1.34e0_wp, 2.1e0_wp, 4.39e0_wp]

    if (iflag > 0) then ! just to avoid the compiler warning
        do i = 1, 15
            tmp1 = i
            tmp2 = 16 - i
            tmp3 = tmp1
            if (i > 8) tmp3 = tmp2
            fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
        end do
    end if

end subroutine fcn

end program example_lmdif1
