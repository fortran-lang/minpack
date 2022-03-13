module minpack_capi
    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_funptr, c_ptr, c_f_procpointer
    use minpack_module
    implicit none
    private

    public :: minpack_hybrd,minpack_hybrd1, minpack_hybrj, minpack_hybrj1, &
        & minpack_lmdif, minpack_lmdif1, minpack_lmder, minpack_lmder1, &
        & minpack_chkder

    public :: minpack_dpmpar

    public :: minpack_func, minpack_func2, minpack_fcn_hybrj, minpack_fcn_lmder, &
        & minpack_fcn_lmstr

    abstract interface
        subroutine minpack_func(n, x, fvec, iflag, udata) bind(c)
            import :: c_int, c_double, c_ptr
            implicit none
            integer(c_int), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out) :: fvec(n)
            integer(c_int), intent(inout) :: iflag
            type(c_ptr), value :: udata
        end subroutine minpack_func

        subroutine minpack_func2(m, n, x, fvec, iflag, udata) bind(c)
            import :: c_int, c_double, c_ptr
            implicit none
            integer(c_int), value :: m
            integer(c_int), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out) :: fvec(m)
            integer(c_int), intent(inout) :: iflag
            type(c_ptr), value :: udata
        end subroutine minpack_func2

        subroutine minpack_fcn_hybrj(n, x, fvec, fjac, ldfjac, iflag, udata) bind(c)
            import :: c_int, c_double, c_ptr
            implicit none
            integer(c_int), value :: n
            real(c_double), intent(in) :: x(n)
            integer(c_int), value :: ldfjac
            real(c_double), intent(out) :: fvec(n)
            real(c_double), intent(out) :: fjac(ldfjac, n)
            integer(c_int), intent(inout) :: iflag
            type(c_ptr), value :: udata
        end subroutine minpack_fcn_hybrj

        subroutine minpack_fcn_lmder(m, n, x, fvec, fjac, ldfjac, iflag, udata) bind(c)
            import :: c_int, c_double, c_ptr
            implicit none
            integer(c_int), value :: m
            integer(c_int), value :: n
            integer(c_int), value :: ldfjac
            integer(c_int), intent(inout) :: iflag
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(inout) :: fvec(m)
            real(c_double), intent(inout) :: fjac(ldfjac, n)
            type(c_ptr), value :: udata
        end subroutine minpack_fcn_lmder

        subroutine minpack_fcn_lmstr(m, n, x, fvec, fjrow, iflag, udata) bind(c)
            import :: c_int, c_double, c_ptr
            implicit none
            integer(c_int), value :: m
            integer(c_int), value :: n
            integer(c_int), intent(inout) :: iflag
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(inout) :: fvec(m)
            real(c_double), intent(inout) :: fjrow(n)
            type(c_ptr), value :: udata
        end subroutine minpack_fcn_lmstr
    end interface

contains

    function minpack_dpmpar(i) result(par) bind(c)
        integer(c_int), value :: i
        real(c_double) :: par
        if (i > 0_c_int .and. i <= 3_c_int) then
            par = dpmpar(i)
        end if
    end function minpack_dpmpar

    subroutine minpack_hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
            & factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4, &
            & udata) &
            & bind(c)
        procedure(minpack_func) :: fcn
        integer(c_int), value :: n
        integer(c_int), value :: maxfev
        integer(c_int), value :: ml
        integer(c_int), value :: mu
        integer(c_int), value :: mode
        integer(c_int), value :: nprint
        integer(c_int), intent(out) :: info
        integer(c_int), intent(out) :: nfev
        integer(c_int), value :: ldfjac
        integer(c_int), value :: lr
        real(c_double), value :: xtol
        real(c_double), value :: epsfcn
        real(c_double), value :: factor
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(n)
        real(c_double), intent(inout) :: diag(n)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(out) :: r(lr)
        real(c_double), intent(out) :: qtf(n)
        real(c_double), intent(inout) :: wa1(n)
        real(c_double), intent(inout) :: wa2(n)
        real(c_double), intent(inout) :: wa3(n)
        real(c_double), intent(inout) :: wa4(n)
        type(c_ptr), value :: udata

        call hybrd(wrap_fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
            & factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4)

    contains
        subroutine wrap_fcn(n, x, fvec, iflag)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: fvec(n)
            integer, intent(inout) :: iflag

            call fcn(n, x, fvec, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_hybrd

    subroutine minpack_hybrd1(fcn, n, x, Fvec, Tol, Info, Wa, Lwa, udata) &
          & bind(c)
        procedure(minpack_func) :: fcn
        integer(c_int), value :: n
        integer(c_int), intent(out) :: info
        real(c_double), value :: tol
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(n)
        integer(c_int), value :: Lwa
        real(c_double), intent(inout) :: Wa(Lwa)
        type(c_ptr), value :: udata

        call hybrd1(wrap_fcn, n, x, fvec, tol, info, Wa, Lwa)

    contains
        subroutine wrap_fcn(n, x, fvec, iflag)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: fvec(n)
            integer, intent(inout) :: iflag

            call fcn(n, x, fvec, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_hybrd1

    subroutine minpack_hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
            & factor, nprint, info, nfev, njev, r, lr, qtf, wa1, wa2, wa3, wa4, udata) &
            & bind(c)
        procedure(minpack_fcn_hybrj) :: fcn
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), value :: maxfev
        integer(c_int), value :: mode
        integer(c_int), value :: nprint
        integer(c_int), intent(out) :: info
        integer(c_int), intent(out) :: nfev
        integer(c_int), intent(out) :: njev
        integer(c_int), value :: lr
        real(c_double), value :: xtol
        real(c_double), value :: factor
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(n)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(inout) :: diag(n)
        real(c_double), intent(out) :: r(lr)
        real(c_double), intent(out) :: qtf(n)
        real(c_double), intent(inout) :: wa1(n)
        real(c_double), intent(inout) :: wa2(n)
        real(c_double), intent(inout) :: wa3(n)
        real(c_double), intent(inout) :: wa4(n)
        type(c_ptr), value :: udata

        call hybrj(wrap_fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
            & factor, nprint, info, nfev, njev, r, lr, qtf, wa1, wa2, wa3, wa4)

    contains
        subroutine wrap_fcn(n, x, fvec, fjac, ldfjac, iflag)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            integer, intent(in) :: ldfjac
            real(wp), intent(inout) :: fvec(n)
            real(wp), intent(inout) :: fjac(ldfjac, n)
            integer, intent(inout) :: iflag

            call fcn(n, x, fvec, fjac, ldfjac, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_hybrj

    subroutine minpack_hybrj1(fcn, n, x, fvec, fjac, ldfjac, tol, info, wa, lwa, udata) &
            & bind(c)
        procedure(minpack_fcn_hybrj) :: fcn
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), intent(out) :: info
        integer(c_int), value :: lwa
        real(c_double), value :: tol
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(n)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(inout) :: wa(lwa)
        type(c_ptr), value :: udata

        call hybrj1(wrap_fcn, n, x, fvec, fjac, ldfjac, tol, info, wa, lwa)

    contains
        subroutine wrap_fcn(n, x, fvec, fjac, ldfjac, iflag)
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            integer, intent(in) :: ldfjac
            real(wp), intent(inout) :: fvec(n)
            real(wp), intent(inout) :: fjac(ldfjac, n)
            integer, intent(inout) :: iflag

            call fcn(n, x, fvec, fjac, ldfjac, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_hybrj1

    subroutine minpack_lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, &
            & mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, &
            & udata) &
            & bind(c)
        procedure(minpack_func2) :: fcn
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), value :: maxfev
        integer(c_int), value :: mode
        integer(c_int), value :: nprint
        integer(c_int), intent(out) :: info
        integer(c_int), intent(out) :: nfev
        integer(c_int), value :: ldfjac
        integer(c_int), intent(out) :: ipvt(n)
        real(c_double), value :: ftol
        real(c_double), value :: xtol
        real(c_double), value :: gtol
        real(c_double), value :: epsfcn
        real(c_double), value :: factor
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(m)
        real(c_double), intent(inout) :: diag(n)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(out) :: qtf(n)
        real(c_double), intent(inout) :: wa1(n)
        real(c_double), intent(inout) :: wa2(n)
        real(c_double), intent(inout) :: wa3(n)
        real(c_double), intent(inout) :: wa4(m)
        type(c_ptr), value :: udata

        call lmdif(wrap_fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, &
            & mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)

    contains
        subroutine wrap_fcn(m, n, x, fvec, iflag)
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: fvec(m)
            integer, intent(inout) :: iflag

            call fcn(m, n, x, fvec, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_lmdif

    subroutine minpack_lmdif1(fcn, m, n, x, fvec, tol, info, iwa, wa, lwa, udata) &
            & bind(c)
        procedure(minpack_func2) :: fcn
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), intent(out) :: info
        integer(c_int), value :: lwa
        integer(c_int), intent(inout) :: iwa(n)
        real(c_double), value :: tol
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(inout) :: fvec(m)
        real(c_double), intent(inout) :: wa(lwa)
        type(c_ptr), value :: udata

        call lmdif1(wrap_fcn, m, n, x, fvec, tol, info, iwa, wa, lwa)

    contains
        subroutine wrap_fcn(m, n, x, fvec, iflag)
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(wp), intent(in) :: x(n)
            real(wp), intent(out) :: fvec(m)
            integer, intent(inout) :: iflag

            call fcn(m, n, x, fvec, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_lmdif1

    subroutine minpack_lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
            & diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4, &
            & udata) &
            & bind(c)
        procedure(minpack_fcn_lmder) :: fcn
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), value :: maxfev
        integer(c_int), value :: mode
        integer(c_int), value :: nprint
        integer(c_int), intent(out) :: info
        integer(c_int), intent(out) :: nfev
        integer(c_int), intent(out) :: njev
        integer(c_int), intent(out) :: ipvt(n)
        real(c_double), value :: ftol
        real(c_double), value :: xtol
        real(c_double), value :: gtol
        real(c_double), value :: factor
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(m)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(inout) :: diag(n)
        real(c_double), intent(out) :: qtf(n)
        real(c_double), intent(inout) :: wa1(n)
        real(c_double), intent(inout) :: wa2(n)
        real(c_double), intent(inout) :: wa3(n)
        real(c_double), intent(inout) :: wa4(m)
        type(c_ptr), value :: udata

        call lmder(wrap_fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
            & diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4)

    contains
        subroutine wrap_fcn(m, n, x, fvec, fjac, ldfjac, iflag)
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: ldfjac
            integer, intent(inout) :: iflag
            real(wp), intent(in) :: x(n)
            real(wp), intent(inout) :: fvec(m)
            real(wp), intent(inout) :: fjac(ldfjac, n)

            call fcn(m, n, x, fvec, fjac, ldfjac, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_lmder

    subroutine minpack_lmder1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa, &
            & udata) &
            & bind(c)
        procedure(minpack_fcn_lmder) :: fcn
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), intent(out) :: info
        integer(c_int), value :: lwa
        integer(c_int), intent(out) :: ipvt(n)
        real(c_double), value :: tol
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(m)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(inout) :: wa(lwa)
        type(c_ptr), value :: udata

        call lmder1(wrap_fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa)

    contains
        subroutine wrap_fcn(m, n, x, fvec, fjac, ldfjac, iflag)
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: ldfjac
            integer, intent(inout) :: iflag
            real(wp), intent(in) :: x(n)
            real(wp), intent(inout) :: fvec(m)
            real(wp), intent(inout) :: fjac(ldfjac, n)

            call fcn(m, n, x, fvec, fjac, ldfjac, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_lmder1

    subroutine minpack_lmstr(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
            & diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4, &
            & udata) &
            & bind(c)
        procedure(minpack_fcn_lmstr) :: fcn
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), value :: maxfev
        integer(c_int), value :: mode
        integer(c_int), value :: nprint
        integer(c_int), intent(out) :: info
        integer(c_int), intent(out) :: nfev
        integer(c_int), intent(out) :: njev
        integer(c_int), intent(out) :: ipvt(n)
        real(c_double), value :: ftol
        real(c_double), value :: xtol
        real(c_double), value :: gtol
        real(c_double), value :: factor
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(m)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(inout) :: diag(n)
        real(c_double), intent(out) :: qtf(n)
        real(c_double), intent(inout) :: wa1(n)
        real(c_double), intent(inout) :: wa2(n)
        real(c_double), intent(inout) :: wa3(n)
        real(c_double), intent(inout) :: wa4(m)
        type(c_ptr), value :: udata

        call lmstr(wrap_fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
            & diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4)
    contains
        subroutine wrap_fcn(m, n, x, fvec, fjrow, iflag)
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(inout) :: iflag
            real(wp), intent(in) :: x(n)
            real(wp), intent(inout) :: fvec(m)
            real(wp), intent(inout) :: fjrow(n)

            call fcn(m, n, x, fvec, fjrow, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_lmstr

    subroutine minpack_lmstr1(fcn, m, n, x, fvec, fjac, ldfjac, tol, info, ipvt, wa, lwa, &
            & udata) &
            & bind(c)
        procedure(minpack_fcn_lmstr) :: fcn
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), intent(out) :: info
        integer(c_int), value :: lwa
        integer(c_int), intent(out) :: ipvt(n)
        real(c_double), value :: tol
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: fvec(m)
        real(c_double), intent(out) :: fjac(ldfjac, n)
        real(c_double), intent(inout) :: wa(lwa)
        type(c_ptr), value :: udata

        call lmstr1(wrap_fcn, m, n, x, fvec, fjac, ldfjac, tol, info, ipvt, wa, lwa)
    contains
        subroutine wrap_fcn(m, n, x, fvec, fjrow, iflag)
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(inout) :: iflag
            real(wp), intent(in) :: x(n)
            real(wp), intent(inout) :: fvec(m)
            real(wp), intent(inout) :: fjrow(n)

            call fcn(m, n, x, fvec, fjrow, iflag, udata)
        end subroutine wrap_fcn
    end subroutine minpack_lmstr1

    subroutine minpack_chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err) &
         & bind(c)
        integer(c_int), value :: m
        integer(c_int), value :: n
        integer(c_int), value :: ldfjac
        integer(c_int), value :: mode
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: fvec(m)
        real(c_double), intent(in) :: fjac(ldfjac, n)
        real(c_double), intent(out) :: xp(n)
        real(c_double), intent(in) :: fvecp(m)
        real(c_double), intent(out) :: err(m)

        call chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err)
    end subroutine minpack_chkder

end module minpack_capi
