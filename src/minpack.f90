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

    !> The purpose of `hybrd` is to find a zero of a system of N non-
    !>  linear functions in N variables by a modification of the Powell
    !>  hybrid method.  The user must provide a subroutine which calcu-
    !>  lates the functions.  The Jacobian is then calculated by a for-
    !>  ward-difference approximation.
    subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
                         mode,factor,nprint,info,nfev,fjac,ldfjac, &
                         r,lr,qtf,wa1,wa2,wa3,wa4)
        integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
        double precision xtol,epsfcn,factor
        double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),qtf(n), &
                         wa1(n),wa2(n),wa3(n),wa4(n)
        interface
            subroutine fcn(n,x,fvec,iflag)
                integer n,iflag
                double precision x(n),fvec(n)
            end subroutine fcn
        end interface
    end subroutine hybrd

    !> The purpose of `hybrd1` is to find a zero of a system of
    !>  n nonlinear functions in n variables by a modification
    !>  of the powell hybrid method.
    subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
        integer n,info,lwa
        double precision tol
        double precision x(n),fvec(n),wa(lwa)
        interface
            subroutine fcn(n,x,fvec,iflag)
                integer n,iflag
                double precision x(n),fvec(n)
            end subroutine fcn
        end interface
    end subroutine hybrd1

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
