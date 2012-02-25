module minpack
implicit none

interface

    double precision function dpmpar(i)
    integer i
    end function

    double precision function enorm(n,x)
    integer n
    double precision x(n)
    end function

    subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
    integer m,n,ldfjac,info,lwa
    integer ipvt(n)
    double precision tol
    double precision x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
    interface
        subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
        implicit none
        integer, intent(in) :: m,n,ldfjac,iflag
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: fvec(m),fjac(ldfjac,n)
        end subroutine
    end interface
    end subroutine

    subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
    integer m,n,info,lwa
    integer iwa(n)
    double precision tol
    double precision x(n),fvec(m),wa(lwa)
    interface
        subroutine fcn(m,n,x,fvec,iflag)
        implicit none
        integer, intent(in) :: m,n,iflag
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: fvec(m)
        end subroutine
    end interface
    end subroutine

    subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
    integer m,n,ldfjac,mode
    double precision x(n),fvec(m),fjac(ldfjac,n),xp(n),fvecp(m),err(m)
    end subroutine

end interface

contains

end module
