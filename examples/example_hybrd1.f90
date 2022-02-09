!> the problem is to determine the values of x(1), x(2), ..., x(9),
!>  which solve the system of tridiagonal equations.
!>
!>  (3-2*x(1))*x(1)           -2*x(2)                   = -1
!>          -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!>                              -x(8) + (3-2*x(9))*x(9) = -1
program example_hybrd1

    use minpack_module, only: hybrd1, dpmpar, enorm
    implicit none
    integer j, n, info, lwa, nwrite
    double precision tol, fnorm
    double precision x(9), fvec(9), wa(180)

    !> Logical output unit is assumed to be number 6.
    data nwrite/6/

    n = 9

    !> The following starting values provide a rough solution.
    do j = 1, 9
        x(j) = -1.d0
    end do

    lwa = 180

    !> Set tol to the square root of the machine precision.
    !>  unless high precision solutions are required,
    !>  this is the recommended setting.
    tol = dsqrt(dpmpar(1))

    call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
    fnorm = enorm(n, fvec)
    write (nwrite, 1000) fnorm, info, (x(j), j=1, n)

1000 format(5x, "FINAL L2 NORM OF THE RESIDUALS", d15.7// &
           5x, "EXIT PARAMETER", 16x, i10// &
           5x, "FINAL APPROXIMATE SOLUTION"// &
           (5x, 3d15.7))

    !> Results obtained with different compilers or machines
    !>  may be slightly different.
    !>
    !>> FINAL L2 NORM OF THE RESIDUALS  0.1192636D-07
    !>>
    !>> EXIT PARAMETER                         1
    !>>
    !>> FINAL APPROXIMATE SOLUTION
    !>>
    !>>  -0.5706545D+00 -0.6816283D+00 -0.7017325D+00
    !>>  -0.7042129D+00 -0.7013690D+00 -0.6918656D+00
    !>>  -0.6657920D+00 -0.5960342D+00 -0.4164121D+00

contains

    !> Subroutine fcn for hybrd1 example.
    subroutine fcn(n, x, fvec, iflag)

        implicit none
        integer n, iflag
        double precision x(n), fvec(n)

        integer k
        double precision one, temp, temp1, temp2, three, two, zero
        data zero, one, two, three/0.0d0, 1.0d0, 2.0d0, 3.0d0/

        do k = 1, n
            temp = (three - two*x(k))*x(k)
            temp1 = zero
            if (k /= 1) temp1 = x(k - 1)
            temp2 = zero
            if (k /= n) temp2 = x(k + 1)
            fvec(k) = temp - temp1 - two*temp2 + one
        end do

    end subroutine fcn

end program example_hybrd1
