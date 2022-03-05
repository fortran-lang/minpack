!*****************************************************************************************
!>
!  CHKDER data for test cases.
!
!  These were from `file23` in the original code:
!  https://netlib.org/minpack/ex/file23

module file23_module

    implicit none

    public

    integer,parameter :: ncases = 14
    integer,dimension(ncases),parameter :: nprobs  = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    integer,dimension(ncases),parameter :: ns      = [2,4,2,4,3,9,7,10,10,10,10,10,10,10]

end module file23_module