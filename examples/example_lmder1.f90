module testmod_der1
implicit none
private
public fcn, dp

integer, parameter :: dp=kind(0d0)

contains

subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag)
integer, intent(in) :: m, n, ldfjac
integer, intent(inout) :: iflag
real(dp), intent(in) :: x(n)
real(dp), intent(inout) :: fvec(m), fjac(ldfjac, n)

integer :: i
real(dp) :: tmp1, tmp2, tmp3, tmp4, y(15)
y = [1.4D-1, 1.8D-1, 2.2D-1, 2.5D-1, 2.9D-1, 3.2D-1, 3.5D-1, 3.9D-1, &
        3.7D-1, 5.8D-1, 7.3D-1, 9.6D-1, 1.34D0, 2.1D0, 4.39D0]

if (iflag == 1) then
    do i = 1, 15
        tmp1 = i
        tmp2 = 16 - i
        tmp3 = tmp1
        if (i .gt. 8) tmp3 = tmp2
        fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
    end do
else
    do i = 1, 15
        tmp1 = i
        tmp2 = 16 - i
        tmp3 = tmp1
        if (i .gt. 8) tmp3 = tmp2
        tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
        fjac(i,1) = -1.D0
        fjac(i,2) = tmp1*tmp2/tmp4
        fjac(i,3) = tmp1*tmp3/tmp4
    end do
end if
end subroutine

end module


program example_lmder1
use minpack_module, only: lmder1, chkder
use testmod_der1, only: dp, fcn
implicit none

integer :: info
real(dp) :: tol, x(3), fvec(15), fjac(size(fvec), size(x))
integer :: ipvt(size(x))
real(dp), allocatable :: wa(:)

! The following starting values provide a rough fit.
x = [1._dp, 1._dp, 1._dp]

call check_deriv()

! Set tol to the square root of the machine precision. Unless high precision
! solutions are required, this is the recommended setting.
tol = sqrt(epsilon(1._dp))

allocate(wa(5*size(x) + size(fvec)))
call lmder1(fcn, size(fvec), size(x), x, fvec, fjac, size(fjac, 1), tol, &
    info, ipvt, wa, size(wa))
print 1000, norm2(fvec), info, x
1000 format(5x, 'FINAL L2 NORM OF THE RESIDUALS', d15.7 // &
            5x, 'EXIT PARAMETER', 16x, i10              // &
            5x, 'FINAL APPROXIMATE SOLUTION'            // &
            5x, 3d15.7)

contains

subroutine check_deriv()
integer :: iflag
real(dp) :: xp(size(x)), fvecp(size(fvec)), err(size(fvec))
call chkder(size(fvec), size(x), x, fvec, fjac, size(fjac, 1), xp, fvecp, &
    1, err)
iflag = 1
call fcn(size(fvec), size(x), x, fvec, fjac, size(fjac, 1), iflag)
iflag = 2
call fcn(size(fvec), size(x), x, fvec, fjac, size(fjac, 1), iflag)
iflag = 1
call fcn(size(fvec), size(x), xp, fvecp, fjac, size(fjac, 1), iflag)
call chkder(size(fvec), size(x), x, fvec, fjac, size(fjac, 1), xp, fvecp, &
    2, err)
print *, "Derivatives check (1.0 is correct, 0.0 is incorrect):"
print *, err
end subroutine

end program
