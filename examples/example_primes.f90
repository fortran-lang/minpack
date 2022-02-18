module types
implicit none
private
public dp

integer, parameter :: dp=kind(0d0)

end module

module find_fit_module

! This module contains a general function find_fit() for a nonlinear least
! squares fitting. The function can fit any nonlinear expression to any data.

use minpack_module, only: lmdif1
use types, only: dp
implicit none
private
public find_fit

contains

subroutine find_fit(data_x, data_y, expr, pars)
! Fits the (data_x, data_y) arrays with the function expr(x, pars).
! The user can provide any nonlinear function 'expr' depending on any number of
! parameters 'pars' and it must return the evaluated expression on the
! array 'x'. The arrays 'data_x' and 'data_y' must have the same
! length.
real(dp), intent(in) :: data_x(:), data_y(:)
interface
    function expr(x, pars) result(y)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(:), pars(:)
    real(dp) :: y(size(x))
    end function
end interface
real(dp), intent(inout) :: pars(:)

real(dp) :: tol, fvec(size(data_x))
integer :: iwa(size(pars)), info, m, n
real(dp), allocatable :: wa(:)

tol = sqrt(epsilon(1._dp))
m = size(fvec)
n = size(pars)
allocate(wa(m*n + 5*n + m))
call lmdif1(fcn, m, n, pars, fvec, tol, info, iwa, wa, size(wa))
if (info /= 1) stop "failed to converge"

contains

subroutine fcn(m, n, x, fvec, iflag)
integer, intent(in) :: m, n
integer, intent(inout) :: iflag
real(dp), intent(in) :: x(n)
real(dp), intent(out) :: fvec(m)
! Suppress compiler warning:
fvec(1) = iflag
fvec = data_y - expr(data_x, x)
end subroutine

end subroutine

end module


program example_primes

! Find a nonlinear fit of the form a*x*log(b + c*x) to a list of primes.

use find_fit_module, only: find_fit
use types, only: dp
implicit none

real(dp) :: pars(3)
real(dp), parameter :: y(*) = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, &
    37, 41, 43, 47, 53, 59, 61, 67, 71]
integer :: i
pars = [1._dp, 1._dp, 1._dp]
call find_fit([(real(i, dp), i=1,size(y))], y, expression, pars)
print *, pars

contains

function expression(x, pars) result(y)
real(dp), intent(in) :: x(:), pars(:)
real(dp) :: y(size(x))
real(dp) :: a, b, c
a = pars(1)
b = pars(2)
c = pars(3)
y = a*x*log(b + c*x)
end function

end program
