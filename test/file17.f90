
      program test
      implicit none

!     **********
!
!     THIS PROGRAM TESTS CODES FOR THE LEAST-SQUARES SOLUTION OF
!     M NONLINEAR EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER
!     AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,
!     CALLS THE NONLINEAR LEAST-SQUARES SOLVER, AND FINALLY PRINTS
!     OUT INFORMATION ON THE PERFORMANCE OF THE SOLVER. THIS IS
!     ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE
!     INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE
!     FORMS OF CALLING SEQUENCES USED BY THE FUNCTION AND JACOBIAN
!     SUBROUTINES IN THE VARIOUS NONLINEAR LEAST-SQUARES SOLVERS.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... FCN
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,INITPT,LMDER1,SSQFCN
!
!       FORTRAN-SUPPLIED ... DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer i , ic , info , k , ldfjac , lwa , m , n , NFEv , NJEv ,  &
            & NPRob , nread , ntries , nwrite
      integer iwa(40) , ma(60) , na(60) , nf(60) , nj(60) , np(60) ,    &
            & nx(60)
      double precision factor , fnorm1 , fnorm2 , one , ten , tol
      double precision fjac(65,40) , fnm(60) , fvec(65) , wa(265) ,     &
                     & x(40)
      double precision dpmpar , enorm
      external fcn
      common /refnum/ NPRob , NFEv , NJEv
!
!     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.
!     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.
!
      data nread , nwrite/5 , 6/
!
      data one , ten/1.0d0 , 1.0d1/
      tol = dsqrt(dpmpar(1))
      ldfjac = 65
      lwa = 265
      ic = 0
      n = 40
      m = 65
      ntries = 1
      do NPRob = 1, 20
      if ( NPRob==20) then
         write (nwrite,99002) ic
99002    format ('1SUMMARY OF ',i3,' CALLS TO LMDER1'/)
         write (nwrite,99003)
99003    format (' NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'/)
         do i = 1 , ic
            write (nwrite,99004) np(i) , na(i) , ma(i) , nf(i) , nj(i) ,&
                               & nx(i) , fnm(i)
99004       format (3i5,3i6,1x,d15.7)
         enddo
         stop
      else
         factor = one
         do k = 1 , ntries
            ic = ic + 1
            call initpt(n,x,NPRob,factor)
            call ssqfcn(m,n,x,fvec,NPRob)
            fnorm1 = enorm(m,fvec)
            write (nwrite,99005) NPRob , n , m
99005       format (////5x,' PROBLEM',i5,5x,' DIMENSIONS',2i5,5x//)
            NFEv = 0
            NJEv = 0
            call lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,iwa,wa,lwa)
            call ssqfcn(m,n,x,fvec,NPRob)
            fnorm2 = enorm(m,fvec)
            np(ic) = NPRob
            na(ic) = n
            ma(ic) = m
            nf(ic) = NFEv
            nj(ic) = NJEv
            nx(ic) = info
            fnm(ic) = fnorm2
            write (nwrite,99006) fnorm1 , fnorm2 , NFEv , NJEv , info , &
                               & (x(i),i=1,n)
99006       format (5x,' INITIAL L2 NORM OF THE RESIDUALS',d15.7//5x,   &
                   &' FINAL L2 NORM OF THE RESIDUALS  ',d15.7//5x,      &
                   &' NUMBER OF FUNCTION EVALUATIONS  ',i10//5x,        &
                   &' NUMBER OF JACOBIAN EVALUATIONS  ',i10//5x,        &
                   &' EXIT PARAMETER',18x,i10//5x,                      &
                   &' FINAL APPROXIMATE SOLUTION'//(5x,5d15.7))
            factor = ten*factor
         enddo
      endif
      end do
      end program test

      subroutine fcn(m,n,x,Fvec,Fjac,Ldfjac,Iflag)
      implicit none

      integer m , n , Ldfjac , Iflag
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n)
!     **********
!
!     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
!     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
!     LEAST-SQUARES SOLVER. FCN SHOULD ONLY CALL THE TESTING
!     FUNCTION AND JACOBIAN SUBROUTINES SSQFCN AND SSQJAC WITH
!     THE APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... SSQFCN,SSQJAC
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer NPRob , NFEv , NJEv
      common /refnum/ NPRob , NFEv , NJEv
      if ( Iflag==1 ) call ssqfcn(m,n,x,Fvec,NPRob)
      if ( Iflag==2 ) call ssqjac(m,n,x,Fjac,Ldfjac,NPRob)
      if ( Iflag==1 ) NFEv = NFEv + 1
      if ( Iflag==2 ) NJEv = NJEv + 1
!
!     LAST CARD OF INTERFACE SUBROUTINE FCN.
!
      end

      subroutine ssqjac(m,n,x,Fjac,Ldfjac,Nprob)
      implicit none

      integer m , n , Ldfjac , Nprob
      double precision x(n) , Fjac(Ldfjac,n)
!     **********
!
!     SUBROUTINE SSQJAC
!
!     THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF EIGHTEEN
!     NONLINEAR LEAST SQUARES PROBLEMS. THE PROBLEM DIMENSIONS ARE
!     AS DESCRIBED IN THE PROLOGUE COMMENTS OF SSQFCN.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)
!
!     WHERE
!
!       M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT
!         EXCEED M.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!       FJAC IS AN M BY N OUTPUT ARRAY WHICH CONTAINS THE JACOBIAN
!         MATRIX OF THE NPROB FUNCTION EVALUATED AT X.
!
!       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
!
!       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIN,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer i , ivar , j , k , mm1 , nm1
      double precision c14 , c20 , c29 , c45 , c100 , div , dx , eight ,&
                     & five , four , one , prod , s2 , temp , ten ,     &
                     & three , ti , tmp1 , tmp2 , tmp3 , tmp4 , tpi ,   &
                     & two , zero
      double precision v(11)
      double precision dfloat
      data zero , one , two , three , four , five , eight , ten , c14 , &
         & c20 , c29 , c45 , c100/0.0d0 , 1.0d0 , 2.0d0 , 3.0d0 ,       &
         & 4.0d0 , 5.0d0 , 8.0d0 , 1.0d1 , 1.4d1 , 2.0d1 , 2.9d1 ,      &
         & 4.5d1 , 1.0d2/
      data v(1) , v(2) , v(3) , v(4) , v(5) , v(6) , v(7) , v(8) ,      &
         & v(9) , v(10) , v(11)/4.0d0 , 2.0d0 , 1.0d0 , 5.0d-1 ,        &
         & 2.5d-1 , 1.67d-1 , 1.25d-1 , 1.0d-1 , 8.33d-2 , 7.14d-2 ,    &
         & 6.25d-2/
      dfloat(ivar) = ivar
!
!     JACOBIAN ROUTINE SELECTOR.
!
      select case (Nprob)
      case (2)
!
!     LINEAR FUNCTION - RANK 1.
!
         do j = 1 , n
            do i = 1 , m
               Fjac(i,j) = dfloat(i)*dfloat(j)
            enddo
         enddo
      case (3)
!
!     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
!
         do j = 1 , n
            do i = 1 , m
               Fjac(i,j) = zero
            enddo
         enddo
         nm1 = n - 1
         mm1 = m - 1
         if ( nm1>=2 ) then
            do j = 2 , nm1
               do i = 2 , mm1
                  Fjac(i,j) = dfloat(i-1)*dfloat(j)
               enddo
            enddo
         endif
      case (4)
!
!     ROSENBROCK FUNCTION.
!
         Fjac(1,1) = -c20*x(1)
         Fjac(1,2) = ten
         Fjac(2,1) = -one
         Fjac(2,2) = zero
      case (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*datan(one)
         temp = x(1)**2 + x(2)**2
         tmp1 = tpi*temp
         tmp2 = dsqrt(temp)
         Fjac(1,1) = c100*x(2)/tmp1
         Fjac(1,2) = -c100*x(1)/tmp1
         Fjac(1,3) = ten
         Fjac(2,1) = ten*x(1)/tmp2
         Fjac(2,2) = ten*x(2)/tmp2
         Fjac(2,3) = zero
         Fjac(3,1) = zero
         Fjac(3,2) = zero
         Fjac(3,3) = one
      case (6)
!
!     POWELL SINGULAR FUNCTION.
!
         do j = 1 , 4
            do i = 1 , 4
               Fjac(i,j) = zero
            enddo
         enddo
         Fjac(1,1) = one
         Fjac(1,2) = ten
         Fjac(2,3) = dsqrt(five)
         Fjac(2,4) = -Fjac(2,3)
         Fjac(3,2) = two*(x(2)-two*x(3))
         Fjac(3,3) = -two*Fjac(3,2)
         Fjac(4,1) = two*dsqrt(ten)*(x(1)-x(4))
         Fjac(4,4) = -Fjac(4,1)
      case (7)
!
!     FREUDENSTEIN AND ROTH FUNCTION.
!
         Fjac(1,1) = one
         Fjac(1,2) = x(2)*(ten-three*x(2)) - two
         Fjac(2,1) = one
         Fjac(2,2) = x(2)*(two+three*x(2)) - c14
      case (8)
!
!     BARD FUNCTION.
!
         do i = 1 , 15
            tmp1 = dfloat(i)
            tmp2 = dfloat(16-i)
            tmp3 = tmp1
            if ( i>8 ) tmp3 = tmp2
            tmp4 = (x(2)*tmp2+x(3)*tmp3)**2
            Fjac(i,1) = -one
            Fjac(i,2) = tmp1*tmp2/tmp4
            Fjac(i,3) = tmp1*tmp3/tmp4
         enddo
      case (9)
!
!     KOWALIK AND OSBORNE FUNCTION.
!
         do i = 1 , 11
            tmp1 = v(i)*(v(i)+x(2))
            tmp2 = v(i)*(v(i)+x(3)) + x(4)
            Fjac(i,1) = -tmp1/tmp2
            Fjac(i,2) = -v(i)*x(1)/tmp2
            Fjac(i,3) = Fjac(i,1)*Fjac(i,2)
            Fjac(i,4) = Fjac(i,3)/v(i)
         enddo
      case (10)
!
!     MEYER FUNCTION.
!
         do i = 1 , 16
            temp = five*dfloat(i) + c45 + x(3)
            tmp1 = x(2)/temp
            tmp2 = dexp(tmp1)
            Fjac(i,1) = tmp2
            Fjac(i,2) = x(1)*tmp2/temp
            Fjac(i,3) = -tmp1*Fjac(i,2)
         enddo
      case (11)
!
!     WATSON FUNCTION.
!
         do i = 1 , 29
            div = dfloat(i)/c29
            s2 = zero
            dx = one
            do j = 1 , n
               s2 = s2 + dx*x(j)
               dx = div*dx
            enddo
            temp = two*div*s2
            dx = one/div
            do j = 1 , n
               Fjac(i,j) = dx*(dfloat(j-1)-temp)
               dx = div*dx
            enddo
         enddo
         do j = 1 , n
            do i = 30 , 31
               Fjac(i,j) = zero
            enddo
         enddo
         Fjac(30,1) = one
         Fjac(31,1) = -two*x(1)
         Fjac(31,2) = one
      case (12)
!
!     BOX 3-DIMENSIONAL FUNCTION.
!
         do i = 1 , m
            temp = dfloat(i)
            tmp1 = temp/ten
            Fjac(i,1) = -tmp1*dexp(-tmp1*x(1))
            Fjac(i,2) = tmp1*dexp(-tmp1*x(2))
            Fjac(i,3) = dexp(-temp) - dexp(-tmp1)
         enddo
      case (13)
!
!     JENNRICH AND SAMPSON FUNCTION.
!
         do i = 1 , m
            temp = dfloat(i)
            Fjac(i,1) = -temp*dexp(temp*x(1))
            Fjac(i,2) = -temp*dexp(temp*x(2))
         enddo
      case (14)
!
!     BROWN AND DENNIS FUNCTION.
!
         do i = 1 , m
            temp = dfloat(i)/five
            ti = dsin(temp)
            tmp1 = x(1) + temp*x(2) - dexp(temp)
            tmp2 = x(3) + ti*x(4) - dcos(temp)
            Fjac(i,1) = two*tmp1
            Fjac(i,2) = temp*Fjac(i,1)
            Fjac(i,3) = two*tmp2
            Fjac(i,4) = ti*Fjac(i,3)
         enddo
      case (15)
!
!     CHEBYQUAD FUNCTION.
!
         dx = one/dfloat(n)
         do j = 1 , n
            tmp1 = one
            tmp2 = two*x(j) - one
            temp = two*tmp2
            tmp3 = zero
            tmp4 = two
            do i = 1 , m
               Fjac(i,j) = dx*tmp4
               ti = four*tmp2 + temp*tmp4 - tmp3
               tmp3 = tmp4
               tmp4 = ti
               ti = temp*tmp2 - tmp1
               tmp1 = tmp2
               tmp2 = ti
            enddo
         enddo
      case (16)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         prod = one
         do j = 1 , n
            prod = x(j)*prod
            do i = 1 , n
               Fjac(i,j) = one
            enddo
            Fjac(j,j) = two
         enddo
         do j = 1 , n
            temp = x(j)
            if ( temp==zero ) then
               temp = one
               prod = one
               do k = 1 , n
                  if ( k/=j ) prod = x(k)*prod
               enddo
            endif
            Fjac(n,j) = prod/temp
         enddo
      case (17)
!
!     OSBORNE 1 FUNCTION.
!
         do i = 1 , 33
            temp = ten*dfloat(i-1)
            tmp1 = dexp(-x(4)*temp)
            tmp2 = dexp(-x(5)*temp)
            Fjac(i,1) = -one
            Fjac(i,2) = -tmp1
            Fjac(i,3) = -tmp2
            Fjac(i,4) = temp*x(2)*tmp1
            Fjac(i,5) = temp*x(3)*tmp2
         enddo
      case (18)
!
!     OSBORNE 2 FUNCTION.
!
         do i = 1 , 65
            temp = dfloat(i-1)/ten
            tmp1 = dexp(-x(5)*temp)
            tmp2 = dexp(-x(6)*(temp-x(9))**2)
            tmp3 = dexp(-x(7)*(temp-x(10))**2)
            tmp4 = dexp(-x(8)*(temp-x(11))**2)
            Fjac(i,1) = -tmp1
            Fjac(i,2) = -tmp2
            Fjac(i,3) = -tmp3
            Fjac(i,4) = -tmp4
            Fjac(i,5) = temp*x(1)*tmp1
            Fjac(i,6) = x(2)*(temp-x(9))**2*tmp2
            Fjac(i,7) = x(3)*(temp-x(10))**2*tmp3
            Fjac(i,8) = x(4)*(temp-x(11))**2*tmp4
            Fjac(i,9) = -two*x(2)*x(6)*(temp-x(9))*tmp2
            Fjac(i,10) = -two*x(3)*x(7)*(temp-x(10))*tmp3
            Fjac(i,11) = -two*x(4)*x(8)*(temp-x(11))*tmp4
         enddo
      case default
!
!     LINEAR FUNCTION - FULL RANK.
!
         temp = two/dfloat(m)
         do j = 1 , n
            do i = 1 , m
               Fjac(i,j) = -temp
            enddo
            Fjac(j,j) = Fjac(j,j) + one
         enddo
      end select
!
!     LAST CARD OF SUBROUTINE SSQJAC.
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
!     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
!     FUNCTIONS DEFINED BY SUBROUTINE SSQFCN. THE SUBROUTINE RETURNS
!     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
!     THE 11TH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
!     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
!     THE VECTOR  X(J) = FACTOR, J=1,...,N.
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
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
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
      double precision c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 ,     &
                     & c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 ,  &
                     & five , h , half , one , seven , ten , three ,    &
                     & twenty , twntf , two , zero
      double precision dfloat
      data zero , half , one , two , three , five , seven , ten ,       &
         & twenty , twntf/0.0d0 , 5.0d-1 , 1.0d0 , 2.0d0 , 3.0d0 ,      &
         & 5.0d0 , 7.0d0 , 1.0d1 , 2.0d1 , 2.5d1/
      data c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 ,     &
         & c12 , c13 , c14 , c15 , c16 , c17/1.2d0 , 2.5d-1 , 3.9d-1 ,  &
         & 4.15d-1 , 2.0d-2 , 4.0d3 , 2.5d2 , 3.0d-1 , 4.0d-1 , 1.5d0 , &
         & 1.0d-2 , 1.3d0 , 6.5d-1 , 7.0d-1 , 6.0d-1 , 4.5d0 , 5.5d0/
      dfloat(ivar) = ivar
!
!     SELECTION OF INITIAL POINT.
!
      select case (Nprob)
      case (4)
!
!     ROSENBROCK FUNCTION.
!
         x(1) = -c1
         x(2) = one
      case (5)
!
!     HELICAL VALLEY FUNCTION.
!
         x(1) = -one
         x(2) = zero
         x(3) = zero
      case (6)
!
!     POWELL SINGULAR FUNCTION.
!
         x(1) = three
         x(2) = -one
         x(3) = zero
         x(4) = one
      case (7)
!
!     FREUDENSTEIN AND ROTH FUNCTION.
!
         x(1) = half
         x(2) = -two
      case (8)
!
!     BARD FUNCTION.
!
         x(1) = one
         x(2) = one
         x(3) = one
      case (9)
!
!     KOWALIK AND OSBORNE FUNCTION.
!
         x(1) = c2
         x(2) = c3
         x(3) = c4
         x(4) = c3
      case (10)
!
!     MEYER FUNCTION.
!
         x(1) = c5
         x(2) = c6
         x(3) = c7
      case (11)
!
!     WATSON FUNCTION.
!
         do j = 1 , n
            x(j) = zero
         enddo
      case (12)
!
!     BOX 3-DIMENSIONAL FUNCTION.
!
         x(1) = zero
         x(2) = ten
         x(3) = twenty
      case (13)
!
!     JENNRICH AND SAMPSON FUNCTION.
!
         x(1) = c8
         x(2) = c9
      case (14)
!
!     BROWN AND DENNIS FUNCTION.
!
         x(1) = twntf
         x(2) = five
         x(3) = -five
         x(4) = -one
      case (15)
!
!     CHEBYQUAD FUNCTION.
!
         h = one/dfloat(n+1)
         do j = 1 , n
            x(j) = dfloat(j)*h
         enddo
      case (16)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         do j = 1 , n
            x(j) = half
         enddo
      case (17)
!
!     OSBORNE 1 FUNCTION.
!
         x(1) = half
         x(2) = c10
         x(3) = -one
         x(4) = c11
         x(5) = c5
      case (18)
!
!     OSBORNE 2 FUNCTION.
!
         x(1) = c12
         x(2) = c13
         x(3) = c13
         x(4) = c14
         x(5) = c15
         x(6) = three
         x(7) = five
         x(8) = seven
         x(9) = two
         x(10) = c16
         x(11) = c17
      case default
!
!     LINEAR FUNCTION - FULL RANK OR RANK 1.
!
         do j = 1 , n
            x(j) = one
         enddo
      end select
!
!     COMPUTE MULTIPLE OF INITIAL POINT.
!
      if ( Factor/=one ) then
         if ( Nprob==11 ) then
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

      subroutine ssqfcn(m,n,x,Fvec,Nprob)
      implicit none

      integer m , n , Nprob
      double precision x(n) , Fvec(m)
!     **********
!
!     SUBROUTINE SSQFCN
!
!     THIS SUBROUTINE DEFINES THE FUNCTIONS OF EIGHTEEN NONLINEAR
!     LEAST SQUARES PROBLEMS. THE ALLOWABLE VALUES OF (M,N) FOR
!     FUNCTIONS 1,2 AND 3 ARE VARIABLE BUT WITH M .GE. N.
!     FOR FUNCTIONS 4,5,6,7,8,9 AND 10 THE VALUES OF (M,N) ARE
!     (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) AND (16,3), RESPECTIVELY.
!     FUNCTION 11 (WATSON) HAS M = 31 WITH N USUALLY 6 OR 9.
!     HOWEVER, ANY N, N = 2,...,31, IS PERMITTED.
!     FUNCTIONS 12,13 AND 14 HAVE N = 3,2 AND 4, RESPECTIVELY, BUT
!     ALLOW ANY M .GE. N, WITH THE USUAL CHOICES BEING 10,10 AND 20.
!     FUNCTION 15 (CHEBYQUAD) ALLOWS M AND N VARIABLE WITH M .GE. N.
!     FUNCTION 16 (BROWN) ALLOWS N VARIABLE WITH M = N.
!     FOR FUNCTIONS 17 AND 18, THE VALUES OF (M,N) ARE
!     (33,5) AND (65,11), RESPECTIVELY.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE SSQFCN(M,N,X,FVEC,NPROB)
!
!     WHERE
!
!       M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT
!         EXCEED M.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS THE NPROB
!         FUNCTION EVALUATED AT X.
!
!       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIN,DSQRT,DSIGN
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      integer i , iev , ivar , j , nm1
      double precision c13 , c14 , c29 , c45 , div , dx , eight , five ,&
                     & one , prod , sum , s1 , s2 , temp , ten , ti ,   &
                     & tmp1 , tmp2 , tmp3 , tmp4 , tpi , two , zero ,   &
                     & zp25 , zp5
      double precision v(11) , y1(15) , y2(11) , y3(16) , y4(33) ,      &
                     & y5(65)
      double precision dfloat
      data zero , zp25 , zp5 , one , two , five , eight , ten , c13 ,   &
         & c14 , c29 , c45/0.0d0 , 2.5d-1 , 5.0d-1 , 1.0d0 , 2.0d0 ,    &
         & 5.0d0 , 8.0d0 , 1.0d1 , 1.3d1 , 1.4d1 , 2.9d1 , 4.5d1/
      data v(1) , v(2) , v(3) , v(4) , v(5) , v(6) , v(7) , v(8) ,      &
         & v(9) , v(10) , v(11)/4.0d0 , 2.0d0 , 1.0d0 , 5.0d-1 ,        &
         & 2.5d-1 , 1.67d-1 , 1.25d-1 , 1.0d-1 , 8.33d-2 , 7.14d-2 ,    &
         & 6.25d-2/
      data y1(1) , y1(2) , y1(3) , y1(4) , y1(5) , y1(6) , y1(7) ,      &
         & y1(8) , y1(9) , y1(10) , y1(11) , y1(12) , y1(13) , y1(14) , &
         & y1(15)/1.4d-1 , 1.8d-1 , 2.2d-1 , 2.5d-1 , 2.9d-1 , 3.2d-1 , &
         & 3.5d-1 , 3.9d-1 , 3.7d-1 , 5.8d-1 , 7.3d-1 , 9.6d-1 ,        &
         & 1.34d0 , 2.1d0 , 4.39d0/
      data y2(1) , y2(2) , y2(3) , y2(4) , y2(5) , y2(6) , y2(7) ,      &
         & y2(8) , y2(9) , y2(10) , y2(11)/1.957d-1 , 1.947d-1 ,        &
         & 1.735d-1 , 1.6d-1 , 8.44d-2 , 6.27d-2 , 4.56d-2 , 3.42d-2 ,  &
         & 3.23d-2 , 2.35d-2 , 2.46d-2/
      data y3(1) , y3(2) , y3(3) , y3(4) , y3(5) , y3(6) , y3(7) ,      &
         & y3(8) , y3(9) , y3(10) , y3(11) , y3(12) , y3(13) , y3(14) , &
         & y3(15) , y3(16)/3.478d4 , 2.861d4 , 2.365d4 , 1.963d4 ,      &
         & 1.637d4 , 1.372d4 , 1.154d4 , 9.744d3 , 8.261d3 , 7.03d3 ,   &
         & 6.005d3 , 5.147d3 , 4.427d3 , 3.82d3 , 3.307d3 , 2.872d3/
      data y4(1) , y4(2) , y4(3) , y4(4) , y4(5) , y4(6) , y4(7) ,      &
         & y4(8) , y4(9) , y4(10) , y4(11) , y4(12) , y4(13) , y4(14) , &
         & y4(15) , y4(16) , y4(17) , y4(18) , y4(19) , y4(20) , y4(21) &
         & , y4(22) , y4(23) , y4(24) , y4(25) , y4(26) , y4(27) ,      &
         & y4(28) , y4(29) , y4(30) , y4(31) , y4(32) , y4(33)/8.44d-1 ,&
         & 9.08d-1 , 9.32d-1 , 9.36d-1 , 9.25d-1 , 9.08d-1 , 8.81d-1 ,  &
         & 8.5d-1 , 8.18d-1 , 7.84d-1 , 7.51d-1 , 7.18d-1 , 6.85d-1 ,   &
         & 6.58d-1 , 6.28d-1 , 6.03d-1 , 5.8d-1 , 5.58d-1 , 5.38d-1 ,   &
         & 5.22d-1 , 5.06d-1 , 4.9d-1 , 4.78d-1 , 4.67d-1 , 4.57d-1 ,   &
         & 4.48d-1 , 4.38d-1 , 4.31d-1 , 4.24d-1 , 4.2d-1 , 4.14d-1 ,   &
         & 4.11d-1 , 4.06d-1/
      data y5(1) , y5(2) , y5(3) , y5(4) , y5(5) , y5(6) , y5(7) ,      &
         & y5(8) , y5(9) , y5(10) , y5(11) , y5(12) , y5(13) , y5(14) , &
         & y5(15) , y5(16) , y5(17) , y5(18) , y5(19) , y5(20) , y5(21) &
         & , y5(22) , y5(23) , y5(24) , y5(25) , y5(26) , y5(27) ,      &
         & y5(28) , y5(29) , y5(30) , y5(31) , y5(32) , y5(33) , y5(34) &
         & , y5(35) , y5(36) , y5(37) , y5(38) , y5(39) , y5(40) ,      &
         & y5(41) , y5(42) , y5(43) , y5(44) , y5(45) , y5(46) , y5(47) &
         & , y5(48) , y5(49) , y5(50) , y5(51) , y5(52) , y5(53) ,      &
         & y5(54) , y5(55) , y5(56) , y5(57) , y5(58) , y5(59) , y5(60) &
         & , y5(61) , y5(62) , y5(63) , y5(64) , y5(65)/1.366d0 ,       &
         & 1.191d0 , 1.112d0 , 1.013d0 , 9.91d-1 , 8.85d-1 , 8.31d-1 ,  &
         & 8.47d-1 , 7.86d-1 , 7.25d-1 , 7.46d-1 , 6.79d-1 , 6.08d-1 ,  &
         & 6.55d-1 , 6.16d-1 , 6.06d-1 , 6.02d-1 , 6.26d-1 , 6.51d-1 ,  &
         & 7.24d-1 , 6.49d-1 , 6.49d-1 , 6.94d-1 , 6.44d-1 , 6.24d-1 ,  &
         & 6.61d-1 , 6.12d-1 , 5.58d-1 , 5.33d-1 , 4.95d-1 , 5.0d-1 ,   &
         & 4.23d-1 , 3.95d-1 , 3.75d-1 , 3.72d-1 , 3.91d-1 , 3.96d-1 ,  &
         & 4.05d-1 , 4.28d-1 , 4.29d-1 , 5.23d-1 , 5.62d-1 , 6.07d-1 ,  &
         & 6.53d-1 , 6.72d-1 , 7.08d-1 , 6.33d-1 , 6.68d-1 , 6.45d-1 ,  &
         & 6.32d-1 , 5.91d-1 , 5.59d-1 , 5.97d-1 , 6.25d-1 , 7.39d-1 ,  &
         & 7.1d-1 , 7.29d-1 , 7.2d-1 , 6.36d-1 , 5.81d-1 , 4.28d-1 ,    &
         & 2.92d-1 , 1.62d-1 , 9.8d-2 , 5.4d-2/
      dfloat(ivar) = ivar
!
!     FUNCTION ROUTINE SELECTOR.
!
      select case (Nprob)
      case (2)
!
!     LINEAR FUNCTION - RANK 1.
!
         sum = zero
         do j = 1 , n
            sum = sum + dfloat(j)*x(j)
         enddo
         do i = 1 , m
            Fvec(i) = dfloat(i)*sum - one
         enddo
      case (3)
!
!     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
!
         sum = zero
         nm1 = n - 1
         if ( nm1>=2 ) then
            do j = 2 , nm1
               sum = sum + dfloat(j)*x(j)
            enddo
         endif
         do i = 1 , m
            Fvec(i) = dfloat(i-1)*sum - one
         enddo
         Fvec(m) = -one
      case (4)
!
!     ROSENBROCK FUNCTION.
!
         Fvec(1) = ten*(x(2)-x(1)**2)
         Fvec(2) = one - x(1)
      case (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*datan(one)
         tmp1 = dsign(zp25,x(2))
         if ( x(1)>zero ) tmp1 = datan(x(2)/x(1))/tpi
         if ( x(1)<zero ) tmp1 = datan(x(2)/x(1))/tpi + zp5
         tmp2 = dsqrt(x(1)**2+x(2)**2)
         Fvec(1) = ten*(x(3)-ten*tmp1)
         Fvec(2) = ten*(tmp2-one)
         Fvec(3) = x(3)
      case (6)
!
!     POWELL SINGULAR FUNCTION.
!
         Fvec(1) = x(1) + ten*x(2)
         Fvec(2) = dsqrt(five)*(x(3)-x(4))
         Fvec(3) = (x(2)-two*x(3))**2
         Fvec(4) = dsqrt(ten)*(x(1)-x(4))**2
      case (7)
!
!     FREUDENSTEIN AND ROTH FUNCTION.
!
         Fvec(1) = -c13 + x(1) + ((five-x(2))*x(2)-two)*x(2)
         Fvec(2) = -c29 + x(1) + ((one+x(2))*x(2)-c14)*x(2)
      case (8)
!
!     BARD FUNCTION.
!
         do i = 1 , 15
            tmp1 = dfloat(i)
            tmp2 = dfloat(16-i)
            tmp3 = tmp1
            if ( i>8 ) tmp3 = tmp2
            Fvec(i) = y1(i) - (x(1)+tmp1/(x(2)*tmp2+x(3)*tmp3))
         enddo
      case (9)
!
!     KOWALIK AND OSBORNE FUNCTION.
!
         do i = 1 , 11
            tmp1 = v(i)*(v(i)+x(2))
            tmp2 = v(i)*(v(i)+x(3)) + x(4)
            Fvec(i) = y2(i) - x(1)*tmp1/tmp2
         enddo
      case (10)
!
!     MEYER FUNCTION.
!
         do i = 1 , 16
            temp = five*dfloat(i) + c45 + x(3)
            tmp1 = x(2)/temp
            tmp2 = dexp(tmp1)
            Fvec(i) = x(1)*tmp2 - y3(i)
         enddo
      case (11)
!
!     WATSON FUNCTION.
!
         do i = 1 , 29
            div = dfloat(i)/c29
            s1 = zero
            dx = one
            do j = 2 , n
               s1 = s1 + dfloat(j-1)*dx*x(j)
               dx = div*dx
            enddo
            s2 = zero
            dx = one
            do j = 1 , n
               s2 = s2 + dx*x(j)
               dx = div*dx
            enddo
            Fvec(i) = s1 - s2**2 - one
         enddo
         Fvec(30) = x(1)
         Fvec(31) = x(2) - x(1)**2 - one
      case (12)
!
!     BOX 3-DIMENSIONAL FUNCTION.
!
         do i = 1 , m
            temp = dfloat(i)
            tmp1 = temp/ten
            Fvec(i) = dexp(-tmp1*x(1)) - dexp(-tmp1*x(2))               &
                    & + (dexp(-temp)-dexp(-tmp1))*x(3)
         enddo
      case (13)
!
!     JENNRICH AND SAMPSON FUNCTION.
!
         do i = 1 , m
            temp = dfloat(i)
            Fvec(i) = two + two*temp - dexp(temp*x(1)) - dexp(temp*x(2))
         enddo
      case (14)
!
!     BROWN AND DENNIS FUNCTION.
!
         do i = 1 , m
            temp = dfloat(i)/five
            tmp1 = x(1) + temp*x(2) - dexp(temp)
            tmp2 = x(3) + dsin(temp)*x(4) - dcos(temp)
            Fvec(i) = tmp1**2 + tmp2**2
         enddo
      case (15)
!
!     CHEBYQUAD FUNCTION.
!
         do i = 1 , m
            Fvec(i) = zero
         enddo
         do j = 1 , n
            tmp1 = one
            tmp2 = two*x(j) - one
            temp = two*tmp2
            do i = 1 , m
               Fvec(i) = Fvec(i) + tmp2
               ti = temp*tmp2 - tmp1
               tmp1 = tmp2
               tmp2 = ti
            enddo
         enddo
         dx = one/dfloat(n)
         iev = -1
         do i = 1 , m
            Fvec(i) = dx*Fvec(i)
            if ( iev>0 ) Fvec(i) = Fvec(i) + one/(dfloat(i)**2-one)
            iev = -iev
         enddo
      case (16)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         sum = -dfloat(n+1)
         prod = one
         do j = 1 , n
            sum = sum + x(j)
            prod = x(j)*prod
         enddo
         do i = 1 , n
            Fvec(i) = x(i) + sum
         enddo
         Fvec(n) = prod - one
      case (17)
!
!     OSBORNE 1 FUNCTION.
!
         do i = 1 , 33
            temp = ten*dfloat(i-1)
            tmp1 = dexp(-x(4)*temp)
            tmp2 = dexp(-x(5)*temp)
            Fvec(i) = y4(i) - (x(1)+x(2)*tmp1+x(3)*tmp2)
         enddo
      case (18)
!
!     OSBORNE 2 FUNCTION.
!
         do i = 1 , 65
            temp = dfloat(i-1)/ten
            tmp1 = dexp(-x(5)*temp)
            tmp2 = dexp(-x(6)*(temp-x(9))**2)
            tmp3 = dexp(-x(7)*(temp-x(10))**2)
            tmp4 = dexp(-x(8)*(temp-x(11))**2)
            Fvec(i) = y5(i) - (x(1)*tmp1+x(2)*tmp2+x(3)*tmp3+x(4)*tmp4)
         enddo
      case default
!
!     LINEAR FUNCTION - FULL RANK.
!
         sum = zero
         do j = 1 , n
            sum = sum + x(j)
         enddo
         temp = two*sum/dfloat(m) + one
         do i = 1 , m
            Fvec(i) = -temp
            if ( i<=n ) Fvec(i) = Fvec(i) + x(i)
         enddo
      end select
!
!     LAST CARD OF SUBROUTINE SSQFCN.
!
      end
