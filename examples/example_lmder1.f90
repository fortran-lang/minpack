program example_lmder1

use minpack_module, only: wp, enorm, lmder1, chkder
use iso_fortran_env, only: nwrite => output_unit

implicit none

integer, parameter :: n = 3
integer, parameter :: m = 15
integer, parameter :: lwa = 5*n+m

integer :: info
real(wp) :: tol, x(n), fvec(m), fjac(m,n)
integer :: ipvt(n)
real(wp) :: wa(lwa)

! The following starting values provide a rough fit.
x = [1.0_wp, 1.0_wp, 1.0_wp]

call check_deriv()

! Set tol to the square root of the machine precision. Unless high precision
! solutions are required, this is the recommended setting.
tol = sqrt(epsilon(1._wp))

call lmder1(fcn, m, n, x, fvec, fjac, m, tol, info, ipvt, wa, lwa)

write(nwrite, '(5x,a,d15.7//,5x,a,16x,i10//,5x,a//(5x,3d15.7))') &
        'FINAL L2 NORM OF THE RESIDUALS', enorm(m, fvec), &
        'EXIT PARAMETER', info, &
        'FINAL APPROXIMATE SOLUTION', x

contains

subroutine check_deriv()

    integer :: iflag
    real(wp) :: xp(n), fvecp(m), err(m)

    call chkder(m, n, x, fvec, fjac, m, xp, fvecp, 1, err)
    iflag = 1
    call fcn(m, n, x, fvec, fjac, m, iflag)
    iflag = 2
    call fcn(m, n, x, fvec, fjac, m, iflag)
    iflag = 1
    call fcn(m, n, xp, fvecp, fjac, m, iflag)
    call chkder(m, n, x, fvec, fjac, m, xp, fvecp, 2, err)

    write(nwrite, '(a)') 'Derivatives check (1.0 is correct, 0.0 is incorrect):'
    write(nwrite,'(1p,(5x,3d15.7))') err
    if (any(abs(err-1.0_wp)>epsilon(1.0_wp))) error stop 'Derivative check failed'

end subroutine check_deriv

subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag)

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n)
    real(wp), intent(inout) :: fvec(m)
    real(wp), intent(inout) :: fjac(ldfjac, n)
    integer, intent(in) :: ldfjac
    integer, intent(inout) :: iflag

    integer :: i
    real(wp) :: tmp1, tmp2, tmp3, tmp4

    real(wp),parameter :: y(15) = [1.4e-1_wp, 1.8e-1_wp, 2.2e-1_wp, 2.5e-1_wp, 2.9e-1_wp, &
                                   3.2e-1_wp, 3.5e-1_wp, 3.9e-1_wp, 3.7e-1_wp, 5.8e-1_wp, &
                                   7.3e-1_wp, 9.6e-1_wp, 1.34e0_wp, 2.1e0_wp,  4.39e0_wp]

    if (iflag == 1) then
        do i = 1, 15
            tmp1 = i
            tmp2 = 16 - i
            tmp3 = tmp1
            if (i > 8) tmp3 = tmp2
            fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
        end do
    else
        do i = 1, 15
            tmp1 = i
            tmp2 = 16 - i
            tmp3 = tmp1
            if (i > 8) tmp3 = tmp2
            tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
            fjac(i,1) = -1.0_wp
            fjac(i,2) = tmp1*tmp2/tmp4
            fjac(i,3) = tmp1*tmp3/tmp4
        end do
    end if

    end subroutine fcn

end program
