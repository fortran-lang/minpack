
      program test
      implicit none

!     **********
!
!     THIS PROGRAM TESTS CODES FOR THE SOLUTION OF N NONLINEAR
!     EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER AND AN
!     INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA, CALLS THE
!     NONLINEAR EQUATION SOLVER, AND FINALLY PRINTS OUT INFORMATION
!     ON THE PERFORMANCE OF THE SOLVER. THIS IS ONLY A SAMPLE DRIVER,
!     MANY OTHER DRIVERS ARE POSSIBLE. THE INTERFACE SUBROUTINE FCN
!     IS NECESSARY TO TAKE INTO ACCOUNT THE FORMS OF CALLING
!     SEQUENCES USED BY THE FUNCTION SUBROUTINES IN THE VARIOUS
!     NONLINEAR EQUATION SOLVERS.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... FCN
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,HYBRD1,INITPT,VECFCN
!
!       FORTRAN-SUPPLIED ... DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer i , ic , info , k , lwa , n , NFEv , NPRob , nread ,      &
              ntries , nwrite, icase
      integer na(60) , nf(60) , np(60) , nx(60)
      double precision factor , fnorm1 , fnorm2 , one , ten , tol
      double precision fnm(60) , fvec(40) , wa(2660) , x(40)
      double precision dpmpar , enorm
      external fcn
      common /refnum/ NPRob , NFEv
!
!     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.
!     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.
!
      data nread , nwrite/5 , 6/
!
      data one , ten/1.0d0 , 1.0d1/
      tol = dsqrt(dpmpar(1))
      lwa = 2660
      ic = 0
      n = 5
      ntries = 1
      do NPRob = 1, 16
      if ( NPRob==16 ) then
         write (nwrite,99002) ic
99002    format ('1SUMMARY OF ',i3,' CALLS TO HYBRD1'/)
         write (nwrite,99003)
99003    format (' NPROB   N    NFEV  INFO  FINAL L2 NORM'/)
         do i = 1 , ic
            write (nwrite,99004) np(i) , na(i) , nf(i) , nx(i) , fnm(i)
99004       format (i4,i6,i7,i6,1x,d15.7)
         enddo
         stop
      else
         factor = one
         do k = 1 , ntries
            ic = ic + 1
            call initpt(n,x,NPRob,factor)
            call vecfcn(n,x,fvec,NPRob)
            fnorm1 = enorm(n,fvec)
            write (nwrite,99005) NPRob , n
99005       format (////5x,' PROBLEM',i5,5x,' DIMENSION',i5,5x//)
            NFEv = 0
            call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
            fnorm2 = enorm(n,fvec)
            np(ic) = NPRob
            na(ic) = n
            nf(ic) = NFEv
            nx(ic) = info
            fnm(ic) = fnorm2
            write (nwrite,99006) fnorm1 , fnorm2 , NFEv , info ,        &
                               & (x(i),i=1,n)
99006       format (5x,' INITIAL L2 NORM OF THE RESIDUALS',d15.7//5x,   &
                   &' FINAL L2 NORM OF THE RESIDUALS  ',d15.7//5x,      &
                   &' NUMBER OF FUNCTION EVALUATIONS  ',i10//5x,        &
                   &' EXIT PARAMETER',18x,i10//5x,                      &
                   &' FINAL APPROXIMATE SOLUTION'//(5x,5d15.7))
            factor = ten*factor
         enddo
      endif
      end do
      end program test

      subroutine fcn(n,x,Fvec,Iflag)
      implicit none

      integer n , Iflag
      double precision x(n) , Fvec(n)
!     **********
!
!     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
!     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
!     EQUATION SOLVER. FCN SHOULD ONLY CALL THE TESTING FUNCTION
!     SUBROUTINE VECFCN WITH THE APPROPRIATE VALUE OF PROBLEM
!     NUMBER (NPROB).
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... VECFCN
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer NPRob , NFEv
      common /refnum/ NPRob , NFEv
      call vecfcn(n,x,Fvec,NPRob)
      NFEv = NFEv + 1
!
!     LAST CARD OF INTERFACE SUBROUTINE FCN.
!
      end

      subroutine vecfcn(n,x,Fvec,Nprob)
      implicit none

      integer n , Nprob
      double precision x(n) , Fvec(n)
!     **********
!
!     SUBROUTINE VECFCN
!
!     THIS SUBROUTINE DEFINES FOURTEEN TEST FUNCTIONS. THE FIRST
!     FIVE TEST FUNCTIONS ARE OF DIMENSIONS 2,4,2,4,3, RESPECTIVELY,
!     WHILE THE REMAINING TEST FUNCTIONS ARE OF VARIABLE DIMENSION
!     N FOR ANY N GREATER THAN OR EQUAL TO 1 (PROBLEM 6 IS AN
!     EXCEPTION TO THIS, SINCE IT DOES NOT ALLOW N = 1).
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE VECFCN(N,X,FVEC,NPROB)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!       FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE NPROB
!         FUNCTION VECTOR EVALUATED AT X.
!
!       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIGN,DSIN,DSQRT,
!                            MAX0,MIN0
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer i , iev , ivar , j , k , k1 , k2 , kp1 , ml , mu
      double precision c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 ,     &
                     & eight , five , h , one , prod , sum , sum1 ,     &
                     & sum2 , temp , temp1 , temp2 , ten , three , ti , &
                     & tj , tk , tpi , two , zero
      double precision dfloat
      data zero , one , two , three , five , eight , ten/0.0d0 , 1.0d0 ,&
         & 2.0d0 , 3.0d0 , 5.0d0 , 8.0d0 , 1.0d1/
      data c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9/1.0d4 , 1.0001d0 ,&
         & 2.0d2 , 2.02d1 , 1.98d1 , 1.8d2 , 2.5d-1 , 5.0d-1 , 2.9d1/
      dfloat(ivar) = ivar
!
!     PROBLEM SELECTOR.
!
      select case (Nprob)
      case (2)
!
!     POWELL SINGULAR FUNCTION.
!
         Fvec(1) = x(1) + ten*x(2)
         Fvec(2) = dsqrt(five)*(x(3)-x(4))
         Fvec(3) = (x(2)-two*x(3))**2
         Fvec(4) = dsqrt(ten)*(x(1)-x(4))**2
      case (3)
!
!     POWELL BADLY SCALED FUNCTION.
!
         Fvec(1) = c1*x(1)*x(2) - one
         Fvec(2) = dexp(-x(1)) + dexp(-x(2)) - c2
      case (4)
!
!     WOOD FUNCTION.
!
         temp1 = x(2) - x(1)**2
         temp2 = x(4) - x(3)**2
         Fvec(1) = -c3*x(1)*temp1 - (one-x(1))
         Fvec(2) = c3*temp1 + c4*(x(2)-one) + c5*(x(4)-one)
         Fvec(3) = -c6*x(3)*temp2 - (one-x(3))
         Fvec(4) = c6*temp2 + c4*(x(4)-one) + c5*(x(2)-one)
      case (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*datan(one)
         temp1 = dsign(c7,x(2))
         if ( x(1)>zero ) temp1 = datan(x(2)/x(1))/tpi
         if ( x(1)<zero ) temp1 = datan(x(2)/x(1))/tpi + c8
         temp2 = dsqrt(x(1)**2+x(2)**2)
         Fvec(1) = ten*(x(3)-ten*temp1)
         Fvec(2) = ten*(temp2-one)
         Fvec(3) = x(3)
      case (6)
!
!     WATSON FUNCTION.
!
         do k = 1 , n
            Fvec(k) = zero
         enddo
         do i = 1 , 29
            ti = dfloat(i)/c9
            sum1 = zero
            temp = one
            do j = 2 , n
               sum1 = sum1 + dfloat(j-1)*temp*x(j)
               temp = ti*temp
            enddo
            sum2 = zero
            temp = one
            do j = 1 , n
               sum2 = sum2 + temp*x(j)
               temp = ti*temp
            enddo
            temp1 = sum1 - sum2**2 - one
            temp2 = two*ti*sum2
            temp = one/ti
            do k = 1 , n
               Fvec(k) = Fvec(k) + temp*(dfloat(k-1)-temp2)*temp1
               temp = ti*temp
            enddo
         enddo
         temp = x(2) - x(1)**2 - one
         Fvec(1) = Fvec(1) + x(1)*(one-two*temp)
         Fvec(2) = Fvec(2) + temp
      case (7)
!
!     CHEBYQUAD FUNCTION.
!
         do k = 1 , n
            Fvec(k) = zero
         enddo
         do j = 1 , n
            temp1 = one
            temp2 = two*x(j) - one
            temp = two*temp2
            do i = 1 , n
               Fvec(i) = Fvec(i) + temp2
               ti = temp*temp2 - temp1
               temp1 = temp2
               temp2 = ti
            enddo
         enddo
         tk = one/dfloat(n)
         iev = -1
         do k = 1 , n
            Fvec(k) = tk*Fvec(k)
            if ( iev>0 ) Fvec(k) = Fvec(k) + one/(dfloat(k)**2-one)
            iev = -iev
         enddo
      case (8)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         sum = -dfloat(n+1)
         prod = one
         do j = 1 , n
            sum = sum + x(j)
            prod = x(j)*prod
         enddo
         do k = 1 , n
            Fvec(k) = x(k) + sum
         enddo
         Fvec(n) = prod - one
      case (9)
!
!     DISCRETE BOUNDARY VALUE FUNCTION.
!
         h = one/dfloat(n+1)
         do k = 1 , n
            temp = (x(k)+dfloat(k)*h+one)**3
            temp1 = zero
            if ( k/=1 ) temp1 = x(k-1)
            temp2 = zero
            if ( k/=n ) temp2 = x(k+1)
            Fvec(k) = two*x(k) - temp1 - temp2 + temp*h**2/two
         enddo
      case (10)
!
!     DISCRETE INTEGRAL EQUATION FUNCTION.
!
         h = one/dfloat(n+1)
         do k = 1 , n
            tk = dfloat(k)*h
            sum1 = zero
            do j = 1 , k
               tj = dfloat(j)*h
               temp = (x(j)+tj+one)**3
               sum1 = sum1 + tj*temp
            enddo
            sum2 = zero
            kp1 = k + 1
            if ( n>=kp1 ) then
               do j = kp1 , n
                  tj = dfloat(j)*h
                  temp = (x(j)+tj+one)**3
                  sum2 = sum2 + (one-tj)*temp
               enddo
            endif
            Fvec(k) = x(k) + h*((one-tk)*sum1+tk*sum2)/two
         enddo
      case (11)
!
!     TRIGONOMETRIC FUNCTION.
!
         sum = zero
         do j = 1 , n
            Fvec(j) = dcos(x(j))
            sum = sum + Fvec(j)
         enddo
         do k = 1 , n
            Fvec(k) = dfloat(n+k) - dsin(x(k)) - sum - dfloat(k)*Fvec(k)
         enddo
      case (12)
!
!     VARIABLY DIMENSIONED FUNCTION.
!
         sum = zero
         do j = 1 , n
            sum = sum + dfloat(j)*(x(j)-one)
         enddo
         temp = sum*(one+two*sum**2)
         do k = 1 , n
            Fvec(k) = x(k) - one + dfloat(k)*temp
         enddo
      case (13)
!
!     BROYDEN TRIDIAGONAL FUNCTION.
!
         do k = 1 , n
            temp = (three-two*x(k))*x(k)
            temp1 = zero
            if ( k/=1 ) temp1 = x(k-1)
            temp2 = zero
            if ( k/=n ) temp2 = x(k+1)
            Fvec(k) = temp - temp1 - two*temp2 + one
         enddo
      case (14)
!
!     BROYDEN BANDED FUNCTION.
!
         ml = 5
         mu = 1
         do k = 1 , n
            k1 = max0(1,k-ml)
            k2 = min0(k+mu,n)
            temp = zero
            do j = k1 , k2
               if ( j/=k ) temp = temp + x(j)*(one+x(j))
            enddo
            Fvec(k) = x(k)*(two+five*x(k)**2) + one - temp
         enddo
      case default
!
!     ROSENBROCK FUNCTION.
!
         Fvec(1) = one - x(1)
         Fvec(2) = ten*(x(2)-x(1)**2)
      end select
!
!     LAST CARD OF SUBROUTINE VECFCN.
!
      end

      subroutine initpt(n,x,Nprob,Factor)
      implicit none

      integer n , Nprob
      double precision Factor
      double precision x(n)
!     **********
!
!     SUBROUTINE INITPT
!
!     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR
!     THE FUNCTIONS DEFINED BY SUBROUTINE VECFCN. THE SUBROUTINE
!     RETURNS IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING
!     POINT. FOR THE SIXTH FUNCTION THE STANDARD STARTING POINT IS
!     ZERO, SO IN THIS CASE, IF FACTOR IS NOT UNITY, THEN THE
!     SUBROUTINE RETURNS THE VECTOR  X(J) = FACTOR, J=1,...,N.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE.
!
!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
!         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.
!
!       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
!
!       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
!         THE STANDARD STARTING POINT. IF FACTOR IS UNITY, NO
!         MULTIPLICATION IS PERFORMED.
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer ivar , j
      double precision c1 , h , half , one , three , tj , zero
      double precision dfloat
      data zero , half , one , three , c1/0.0d0 , 5.0d-1 , 1.0d0 ,      &
         & 3.0d0 , 1.2d0/
      dfloat(ivar) = ivar
!
!     SELECTION OF INITIAL POINT.
!
      select case (Nprob)
      case (2)
!
!     POWELL SINGULAR FUNCTION.
!
         x(1) = three
         x(2) = -one
         x(3) = zero
         x(4) = one
      case (3)
!
!     POWELL BADLY SCALED FUNCTION.
!
         x(1) = zero
         x(2) = one
      case (4)
!
!     WOOD FUNCTION.
!
         x(1) = -three
         x(2) = -one
         x(3) = -three
         x(4) = -one
      case (5)
!
!     HELICAL VALLEY FUNCTION.
!
         x(1) = -one
         x(2) = zero
         x(3) = zero
      case (6)
!
!     WATSON FUNCTION.
!
         do j = 1 , n
            x(j) = zero
         enddo
      case (7)
!
!     CHEBYQUAD FUNCTION.
!
         h = one/dfloat(n+1)
         do j = 1 , n
            x(j) = dfloat(j)*h
         enddo
      case (8)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         do j = 1 , n
            x(j) = half
         enddo
      case (9,10)
!
!     DISCRETE BOUNDARY VALUE AND INTEGRAL EQUATION FUNCTIONS.
!
         h = one/dfloat(n+1)
         do j = 1 , n
            tj = dfloat(j)*h
            x(j) = tj*(tj-one)
         enddo
      case (11)
!
!     TRIGONOMETRIC FUNCTION.
!
         h = one/dfloat(n)
         do j = 1 , n
            x(j) = h
         enddo
      case (12)
!
!     VARIABLY DIMENSIONED FUNCTION.
!
         h = one/dfloat(n)
         do j = 1 , n
            x(j) = one - dfloat(j)*h
         enddo
      case (13,14)
!
!     BROYDEN TRIDIAGONAL AND BANDED FUNCTIONS.
!
         do j = 1 , n
            x(j) = -one
         enddo
      case default
!
!     ROSENBROCK FUNCTION.
!
         x(1) = -c1
         x(2) = one
      end select
!
!     COMPUTE MULTIPLE OF INITIAL POINT.
!
      if ( Factor/=one ) then
         if ( Nprob==6 ) then
            do j = 1 , n
               x(j) = Factor
            enddo
         else
            do j = 1 , n
               x(j) = Factor*x(j)
            enddo
         endif
      endif
!
!     LAST CARD OF SUBROUTINE INITPT.
!
      end
