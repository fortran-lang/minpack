!*****************************************************************************************
!>
!  given an n-vector x, this function calculates the
!  euclidean norm of x.
!
!  the euclidean norm is computed by accumulating the sum of
!  squares in three different sums. the sums of squares for the
!  small and large components are scaled so that no overflows
!  occur. non-destructive underflows are permitted. underflows
!  and overflows do not occur in the computation of the unscaled
!  sum of squares for the intermediate components.
!  the definitions of small, intermediate and large components
!  depend on two constants, rdwarf and rgiant. the main
!  restrictions on these constants are that rdwarf**2 not
!  underflow and rgiant**2 not overflow. the constants
!  given here are suitable for every known computer.

    pure real(wp) function enorm(n, x)
        use iso_fortran_env, only: wp => real64
        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        real(wp), intent(in) :: x(n) !! an input array of length n.

        integer :: i
        real(wp) :: agiant, s1, s2, s3, xabs, x1max, x3max

        real(wp), parameter :: rdwarf = 3.834e-20_wp
        real(wp), parameter :: rgiant = 1.304e19_wp

        s1 = zero
        s2 = zero
        s3 = zero
        x1max = zero
        x3max = zero
        agiant = rgiant/real(n, wp)
        do i = 1, n
            xabs = abs(x(i))
            if (xabs > rdwarf .and. xabs < agiant) then
                ! sum for intermediate components.
                s2 = s2 + xabs**2
            elseif (xabs <= rdwarf) then
                ! sum for small components.
                if (xabs <= x3max) then
                    if (xabs /= zero) s3 = s3 + (xabs/x3max)**2
                else
                    s3 = one + s3*(x3max/xabs)**2
                    x3max = xabs
                end if
                ! sum for large components.
            elseif (xabs <= x1max) then
                s1 = s1 + (xabs/x1max)**2
            else
                s1 = one + s1*(x1max/xabs)**2
                x1max = xabs
            end if
        end do

        ! calculation of norm.

        if (s1 /= zero) then
            enorm = x1max*sqrt(s1 + (s2/x1max)/x1max)
        elseif (s2 == zero) then
            enorm = x3max*sqrt(s3)
        else
            if (s2 >= x3max) enorm = sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
            if (s2 < x3max) enorm = sqrt(x3max*((s2/x3max) + (x3max*s3)))
        end if

    end function enorm
!*****************************************************************************************