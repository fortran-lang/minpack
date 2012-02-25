module testmod_dif
implicit none
private
public fcn, dp

integer, parameter :: dp=kind(0d0)

contains

subroutine fcn(m, n, x, fvec, iflag)
integer, intent(in) :: m, n, iflag
real(dp), intent(in) :: x(n)
real(dp), intent(out) :: fvec(m)

integer :: i
real(dp) :: tmp1, tmp2, tmp3, y(15)
! Suppress compiler warning:
y(1) = iflag
y = [1.4D-1, 1.8D-1, 2.2D-1, 2.5D-1, 2.9D-1, 3.2D-1, 3.5D-1, 3.9D-1, &
        3.7D-1, 5.8D-1, 7.3D-1, 9.6D-1, 1.34D0, 2.1D0, 4.39D0]

do i = 1, 15
    tmp1 = i
    tmp2 = 16 - i
    tmp3 = tmp1
    if (i .gt. 8) tmp3 = tmp2
    fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
end do
end subroutine

end module


program example_lmdif1
use minpack, only: enorm, dpmpar, lmdif1
use testmod_dif, only: dp, fcn
implicit none

integer :: info, m, n
real(dp) :: tol, x(3), fvec(15)
integer :: iwa(size(x))
real(dp), allocatable :: wa(:)

! The following starting values provide a rough fit.
x = [1._dp, 1._dp, 1._dp]

! Set tol to the square root of the machine precision. Unless high precision
! solutions are required, this is the recommended setting.
tol = sqrt(dpmpar(1))

m = size(fvec)
n = size(x)
allocate(wa(m*n + 5*n + m))
call lmdif1(fcn, size(fvec), size(x), x, fvec, tol, info, iwa, wa, size(wa))
print 1000, enorm(size(fvec), fvec), info, x
1000 format(5x, 'FINAL L2 NORM OF THE RESIDUALS', d15.7 // &
            5x, 'EXIT PARAMETER', 16x, i10              // &
            5x, 'FINAL APPROXIMATE SOLUTION'            // &
            5x, 3d15.7)
end program
