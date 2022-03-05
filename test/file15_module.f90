!*****************************************************************************************
!>
!  hybrd and hybrj data for test cases.
!
!  These were from `file15` in the original code:
!  https://netlib.org/minpack/ex/file15

module file15_module

    implicit none

    public

    integer,parameter :: ncases = 22
    integer,dimension(ncases),parameter :: nprobs  = [1,2,3,4,5,6,6,7,7,7,7,7,8,8,8,9,10,10,11,12,13,14]
    integer,dimension(ncases),parameter :: ns      = [2,4,2,4,3,6,9,5,6,7,8,9,10,30,40,10,1,10,10,10,10,10]
    integer,dimension(ncases),parameter :: ntriess = [3,3,2,3,3,2,2,3,3,3,1,1,3,1,1,3,3,3,3,3,3,3]

end module file15_module