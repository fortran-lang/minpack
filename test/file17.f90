!*==AA0001.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      PROGRAM TEST
      IMPLICIT NONE
!*--AA00013
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
      INTEGER i , ic , info , k , ldfjac , lwa , m , n , NFEv , NJEv ,  &
            & NPRob , nread , ntries , nwrite
      INTEGER iwa(40) , ma(60) , na(60) , nf(60) , nj(60) , np(60) ,    &
            & nx(60)
      DOUBLE PRECISION factor , fnorm1 , fnorm2 , one , ten , tol
      DOUBLE PRECISION fjac(65,40) , fnm(60) , fvec(65) , wa(265) ,     &
                     & x(40)
      DOUBLE PRECISION DPMPAR , ENORM
      EXTERNAL FCN
      COMMON /REFNUM/ NPRob , NFEv , NJEv
!
!     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.
!     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.
!
      DATA nread , nwrite/5 , 6/
!
      DATA one , ten/1.0D0 , 1.0D1/
      tol = DSQRT(DPMPAR(1))
      ldfjac = 65
      lwa = 265
      ic = 0
 100  READ (nread,99001) NPRob , n , m , ntries
99001 FORMAT (4I5)
      IF ( NPRob<=0 ) THEN
         WRITE (nwrite,99002) ic
99002    FORMAT ('1SUMMARY OF ',I3,' CALLS TO LMDER1'/)
         WRITE (nwrite,99003)
99003    FORMAT (' NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'/)
         DO i = 1 , ic
            WRITE (nwrite,99004) np(i) , na(i) , ma(i) , nf(i) , nj(i) ,&
                               & nx(i) , fnm(i)
99004       FORMAT (3I5,3I6,1X,D15.7)
         ENDDO
         STOP
      ELSE
         factor = one
         DO k = 1 , ntries
            ic = ic + 1
            CALL INITPT(n,x,NPRob,factor)
            CALL SSQFCN(m,n,x,fvec,NPRob)
            fnorm1 = ENORM(m,fvec)
            WRITE (nwrite,99005) NPRob , n , m
99005       FORMAT (////5X,' PROBLEM',I5,5X,' DIMENSIONS',2I5,5X//)
            NFEv = 0
            NJEv = 0
            CALL LMDER1(FCN,m,n,x,fvec,fjac,ldfjac,tol,info,iwa,wa,lwa)
            CALL SSQFCN(m,n,x,fvec,NPRob)
            fnorm2 = ENORM(m,fvec)
            np(ic) = NPRob
            na(ic) = n
            ma(ic) = m
            nf(ic) = NFEv
            nj(ic) = NJEv
            nx(ic) = info
            fnm(ic) = fnorm2
            WRITE (nwrite,99006) fnorm1 , fnorm2 , NFEv , NJEv , info , &
                               & (x(i),i=1,n)
99006       FORMAT (5X,' INITIAL L2 NORM OF THE RESIDUALS',D15.7//5X,   &
                   &' FINAL L2 NORM OF THE RESIDUALS  ',D15.7//5X,      &
                   &' NUMBER OF FUNCTION EVALUATIONS  ',I10//5X,        &
                   &' NUMBER OF JACOBIAN EVALUATIONS  ',I10//5X,        &
                   &' EXIT PARAMETER',18X,I10//5X,                      &
                   &' FINAL APPROXIMATE SOLUTION'//(5X,5D15.7))
            factor = ten*factor
         ENDDO
         GOTO 100
      ENDIF
!
!     LAST CARD OF DRIVER.
!
      END
!*==FCN.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE FCN(M,N,X,Fvec,Fjac,Ldfjac,Iflag)
      IMPLICIT NONE
!*--FCN105
      INTEGER M , N , Ldfjac , Iflag
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N)
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
      INTEGER NPRob , NFEv , NJEv
      COMMON /REFNUM/ NPRob , NFEv , NJEv
      IF ( Iflag==1 ) CALL SSQFCN(M,N,X,Fvec,NPRob)
      IF ( Iflag==2 ) CALL SSQJAC(M,N,X,Fjac,Ldfjac,NPRob)
      IF ( Iflag==1 ) NFEv = NFEv + 1
      IF ( Iflag==2 ) NJEv = NJEv + 1
!
!     LAST CARD OF INTERFACE SUBROUTINE FCN.
!
      END
!*==SSQJAC.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE SSQJAC(M,N,X,Fjac,Ldfjac,Nprob)
      IMPLICIT NONE
!*--SSQJAC137
      INTEGER M , N , Ldfjac , Nprob
      DOUBLE PRECISION X(N) , Fjac(Ldfjac,N)
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
      INTEGER i , ivar , j , k , mm1 , nm1
      DOUBLE PRECISION c14 , c20 , c29 , c45 , c100 , div , dx , eight ,&
                     & five , four , one , prod , s2 , temp , ten ,     &
                     & three , ti , tmp1 , tmp2 , tmp3 , tmp4 , tpi ,   &
                     & two , zero
      DOUBLE PRECISION v(11)
      DOUBLE PRECISION DFLOAT
      DATA zero , one , two , three , four , five , eight , ten , c14 , &
         & c20 , c29 , c45 , c100/0.0D0 , 1.0D0 , 2.0D0 , 3.0D0 ,       &
         & 4.0D0 , 5.0D0 , 8.0D0 , 1.0D1 , 1.4D1 , 2.0D1 , 2.9D1 ,      &
         & 4.5D1 , 1.0D2/
      DATA v(1) , v(2) , v(3) , v(4) , v(5) , v(6) , v(7) , v(8) ,      &
         & v(9) , v(10) , v(11)/4.0D0 , 2.0D0 , 1.0D0 , 5.0D-1 ,        &
         & 2.5D-1 , 1.67D-1 , 1.25D-1 , 1.0D-1 , 8.33D-2 , 7.14D-2 ,    &
         & 6.25D-2/
      DFLOAT(ivar) = ivar
!
!     JACOBIAN ROUTINE SELECTOR.
!
      SELECT CASE (Nprob)
      CASE (2)
!
!     LINEAR FUNCTION - RANK 1.
!
         DO j = 1 , N
            DO i = 1 , M
               Fjac(i,j) = DFLOAT(i)*DFLOAT(j)
            ENDDO
         ENDDO
      CASE (3)
!
!     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
!
         DO j = 1 , N
            DO i = 1 , M
               Fjac(i,j) = zero
            ENDDO
         ENDDO
         nm1 = N - 1
         mm1 = M - 1
         IF ( nm1>=2 ) THEN
            DO j = 2 , nm1
               DO i = 2 , mm1
                  Fjac(i,j) = DFLOAT(i-1)*DFLOAT(j)
               ENDDO
            ENDDO
         ENDIF
      CASE (4)
!
!     ROSENBROCK FUNCTION.
!
         Fjac(1,1) = -c20*X(1)
         Fjac(1,2) = ten
         Fjac(2,1) = -one
         Fjac(2,2) = zero
      CASE (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*DATAN(one)
         temp = X(1)**2 + X(2)**2
         tmp1 = tpi*temp
         tmp2 = DSQRT(temp)
         Fjac(1,1) = c100*X(2)/tmp1
         Fjac(1,2) = -c100*X(1)/tmp1
         Fjac(1,3) = ten
         Fjac(2,1) = ten*X(1)/tmp2
         Fjac(2,2) = ten*X(2)/tmp2
         Fjac(2,3) = zero
         Fjac(3,1) = zero
         Fjac(3,2) = zero
         Fjac(3,3) = one
      CASE (6)
!
!     POWELL SINGULAR FUNCTION.
!
         DO j = 1 , 4
            DO i = 1 , 4
               Fjac(i,j) = zero
            ENDDO
         ENDDO
         Fjac(1,1) = one
         Fjac(1,2) = ten
         Fjac(2,3) = DSQRT(five)
         Fjac(2,4) = -Fjac(2,3)
         Fjac(3,2) = two*(X(2)-two*X(3))
         Fjac(3,3) = -two*Fjac(3,2)
         Fjac(4,1) = two*DSQRT(ten)*(X(1)-X(4))
         Fjac(4,4) = -Fjac(4,1)
      CASE (7)
!
!     FREUDENSTEIN AND ROTH FUNCTION.
!
         Fjac(1,1) = one
         Fjac(1,2) = X(2)*(ten-three*X(2)) - two
         Fjac(2,1) = one
         Fjac(2,2) = X(2)*(two+three*X(2)) - c14
      CASE (8)
!
!     BARD FUNCTION.
!
         DO i = 1 , 15
            tmp1 = DFLOAT(i)
            tmp2 = DFLOAT(16-i)
            tmp3 = tmp1
            IF ( i>8 ) tmp3 = tmp2
            tmp4 = (X(2)*tmp2+X(3)*tmp3)**2
            Fjac(i,1) = -one
            Fjac(i,2) = tmp1*tmp2/tmp4
            Fjac(i,3) = tmp1*tmp3/tmp4
         ENDDO
      CASE (9)
!
!     KOWALIK AND OSBORNE FUNCTION.
!
         DO i = 1 , 11
            tmp1 = v(i)*(v(i)+X(2))
            tmp2 = v(i)*(v(i)+X(3)) + X(4)
            Fjac(i,1) = -tmp1/tmp2
            Fjac(i,2) = -v(i)*X(1)/tmp2
            Fjac(i,3) = Fjac(i,1)*Fjac(i,2)
            Fjac(i,4) = Fjac(i,3)/v(i)
         ENDDO
      CASE (10)
!
!     MEYER FUNCTION.
!
         DO i = 1 , 16
            temp = five*DFLOAT(i) + c45 + X(3)
            tmp1 = X(2)/temp
            tmp2 = DEXP(tmp1)
            Fjac(i,1) = tmp2
            Fjac(i,2) = X(1)*tmp2/temp
            Fjac(i,3) = -tmp1*Fjac(i,2)
         ENDDO
      CASE (11)
!
!     WATSON FUNCTION.
!
         DO i = 1 , 29
            div = DFLOAT(i)/c29
            s2 = zero
            dx = one
            DO j = 1 , N
               s2 = s2 + dx*X(j)
               dx = div*dx
            ENDDO
            temp = two*div*s2
            dx = one/div
            DO j = 1 , N
               Fjac(i,j) = dx*(DFLOAT(j-1)-temp)
               dx = div*dx
            ENDDO
         ENDDO
         DO j = 1 , N
            DO i = 30 , 31
               Fjac(i,j) = zero
            ENDDO
         ENDDO
         Fjac(30,1) = one
         Fjac(31,1) = -two*X(1)
         Fjac(31,2) = one
      CASE (12)
!
!     BOX 3-DIMENSIONAL FUNCTION.
!
         DO i = 1 , M
            temp = DFLOAT(i)
            tmp1 = temp/ten
            Fjac(i,1) = -tmp1*DEXP(-tmp1*X(1))
            Fjac(i,2) = tmp1*DEXP(-tmp1*X(2))
            Fjac(i,3) = DEXP(-temp) - DEXP(-tmp1)
         ENDDO
      CASE (13)
!
!     JENNRICH AND SAMPSON FUNCTION.
!
         DO i = 1 , M
            temp = DFLOAT(i)
            Fjac(i,1) = -temp*DEXP(temp*X(1))
            Fjac(i,2) = -temp*DEXP(temp*X(2))
         ENDDO
      CASE (14)
!
!     BROWN AND DENNIS FUNCTION.
!
         DO i = 1 , M
            temp = DFLOAT(i)/five
            ti = DSIN(temp)
            tmp1 = X(1) + temp*X(2) - DEXP(temp)
            tmp2 = X(3) + ti*X(4) - DCOS(temp)
            Fjac(i,1) = two*tmp1
            Fjac(i,2) = temp*Fjac(i,1)
            Fjac(i,3) = two*tmp2
            Fjac(i,4) = ti*Fjac(i,3)
         ENDDO
      CASE (15)
!
!     CHEBYQUAD FUNCTION.
!
         dx = one/DFLOAT(N)
         DO j = 1 , N
            tmp1 = one
            tmp2 = two*X(j) - one
            temp = two*tmp2
            tmp3 = zero
            tmp4 = two
            DO i = 1 , M
               Fjac(i,j) = dx*tmp4
               ti = four*tmp2 + temp*tmp4 - tmp3
               tmp3 = tmp4
               tmp4 = ti
               ti = temp*tmp2 - tmp1
               tmp1 = tmp2
               tmp2 = ti
            ENDDO
         ENDDO
      CASE (16)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         prod = one
         DO j = 1 , N
            prod = X(j)*prod
            DO i = 1 , N
               Fjac(i,j) = one
            ENDDO
            Fjac(j,j) = two
         ENDDO
         DO j = 1 , N
            temp = X(j)
            IF ( temp==zero ) THEN
               temp = one
               prod = one
               DO k = 1 , N
                  IF ( k/=j ) prod = X(k)*prod
               ENDDO
            ENDIF
            Fjac(N,j) = prod/temp
         ENDDO
      CASE (17)
!
!     OSBORNE 1 FUNCTION.
!
         DO i = 1 , 33
            temp = ten*DFLOAT(i-1)
            tmp1 = DEXP(-X(4)*temp)
            tmp2 = DEXP(-X(5)*temp)
            Fjac(i,1) = -one
            Fjac(i,2) = -tmp1
            Fjac(i,3) = -tmp2
            Fjac(i,4) = temp*X(2)*tmp1
            Fjac(i,5) = temp*X(3)*tmp2
         ENDDO
      CASE (18)
!
!     OSBORNE 2 FUNCTION.
!
         DO i = 1 , 65
            temp = DFLOAT(i-1)/ten
            tmp1 = DEXP(-X(5)*temp)
            tmp2 = DEXP(-X(6)*(temp-X(9))**2)
            tmp3 = DEXP(-X(7)*(temp-X(10))**2)
            tmp4 = DEXP(-X(8)*(temp-X(11))**2)
            Fjac(i,1) = -tmp1
            Fjac(i,2) = -tmp2
            Fjac(i,3) = -tmp3
            Fjac(i,4) = -tmp4
            Fjac(i,5) = temp*X(1)*tmp1
            Fjac(i,6) = X(2)*(temp-X(9))**2*tmp2
            Fjac(i,7) = X(3)*(temp-X(10))**2*tmp3
            Fjac(i,8) = X(4)*(temp-X(11))**2*tmp4
            Fjac(i,9) = -two*X(2)*X(6)*(temp-X(9))*tmp2
            Fjac(i,10) = -two*X(3)*X(7)*(temp-X(10))*tmp3
            Fjac(i,11) = -two*X(4)*X(8)*(temp-X(11))*tmp4
         ENDDO
      CASE DEFAULT
!
!     LINEAR FUNCTION - FULL RANK.
!
         temp = two/DFLOAT(M)
         DO j = 1 , N
            DO i = 1 , M
               Fjac(i,j) = -temp
            ENDDO
            Fjac(j,j) = Fjac(j,j) + one
         ENDDO
      END SELECT
!
!     LAST CARD OF SUBROUTINE SSQJAC.
!
      END
!*==INITPT.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE INITPT(N,X,Nprob,Factor)
      IMPLICIT NONE
!*--INITPT471
      INTEGER N , Nprob
      DOUBLE PRECISION Factor
      DOUBLE PRECISION X(N)
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
      INTEGER ivar , j
      DOUBLE PRECISION c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 ,     &
                     & c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 ,  &
                     & five , h , half , one , seven , ten , three ,    &
                     & twenty , twntf , two , zero
      DOUBLE PRECISION DFLOAT
      DATA zero , half , one , two , three , five , seven , ten ,       &
         & twenty , twntf/0.0D0 , 5.0D-1 , 1.0D0 , 2.0D0 , 3.0D0 ,      &
         & 5.0D0 , 7.0D0 , 1.0D1 , 2.0D1 , 2.5D1/
      DATA c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 ,     &
         & c12 , c13 , c14 , c15 , c16 , c17/1.2D0 , 2.5D-1 , 3.9D-1 ,  &
         & 4.15D-1 , 2.0D-2 , 4.0D3 , 2.5D2 , 3.0D-1 , 4.0D-1 , 1.5D0 , &
         & 1.0D-2 , 1.3D0 , 6.5D-1 , 7.0D-1 , 6.0D-1 , 4.5D0 , 5.5D0/
      DFLOAT(ivar) = ivar
!
!     SELECTION OF INITIAL POINT.
!
      SELECT CASE (Nprob)
      CASE (4)
!
!     ROSENBROCK FUNCTION.
!
         X(1) = -c1
         X(2) = one
      CASE (5)
!
!     HELICAL VALLEY FUNCTION.
!
         X(1) = -one
         X(2) = zero
         X(3) = zero
      CASE (6)
!
!     POWELL SINGULAR FUNCTION.
!
         X(1) = three
         X(2) = -one
         X(3) = zero
         X(4) = one
      CASE (7)
!
!     FREUDENSTEIN AND ROTH FUNCTION.
!
         X(1) = half
         X(2) = -two
      CASE (8)
!
!     BARD FUNCTION.
!
         X(1) = one
         X(2) = one
         X(3) = one
      CASE (9)
!
!     KOWALIK AND OSBORNE FUNCTION.
!
         X(1) = c2
         X(2) = c3
         X(3) = c4
         X(4) = c3
      CASE (10)
!
!     MEYER FUNCTION.
!
         X(1) = c5
         X(2) = c6
         X(3) = c7
      CASE (11)
!
!     WATSON FUNCTION.
!
         DO j = 1 , N
            X(j) = zero
         ENDDO
      CASE (12)
!
!     BOX 3-DIMENSIONAL FUNCTION.
!
         X(1) = zero
         X(2) = ten
         X(3) = twenty
      CASE (13)
!
!     JENNRICH AND SAMPSON FUNCTION.
!
         X(1) = c8
         X(2) = c9
      CASE (14)
!
!     BROWN AND DENNIS FUNCTION.
!
         X(1) = twntf
         X(2) = five
         X(3) = -five
         X(4) = -one
      CASE (15)
!
!     CHEBYQUAD FUNCTION.
!
         h = one/DFLOAT(N+1)
         DO j = 1 , N
            X(j) = DFLOAT(j)*h
         ENDDO
      CASE (16)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         DO j = 1 , N
            X(j) = half
         ENDDO
      CASE (17)
!
!     OSBORNE 1 FUNCTION.
!
         X(1) = half
         X(2) = c10
         X(3) = -one
         X(4) = c11
         X(5) = c5
      CASE (18)
!
!     OSBORNE 2 FUNCTION.
!
         X(1) = c12
         X(2) = c13
         X(3) = c13
         X(4) = c14
         X(5) = c15
         X(6) = three
         X(7) = five
         X(8) = seven
         X(9) = two
         X(10) = c16
         X(11) = c17
      CASE DEFAULT
!
!     LINEAR FUNCTION - FULL RANK OR RANK 1.
!
         DO j = 1 , N
            X(j) = one
         ENDDO
      END SELECT
!
!     COMPUTE MULTIPLE OF INITIAL POINT.
!
      IF ( Factor/=one ) THEN
         IF ( Nprob==11 ) THEN
            DO j = 1 , N
               X(j) = Factor
            ENDDO
         ELSE
            DO j = 1 , N
               X(j) = Factor*X(j)
            ENDDO
         ENDIF
      ENDIF
!
!     LAST CARD OF SUBROUTINE INITPT.
!
      END
!*==SSQFCN.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE SSQFCN(M,N,X,Fvec,Nprob)
      IMPLICIT NONE
!*--SSQFCN671
      INTEGER M , N , Nprob
      DOUBLE PRECISION X(N) , Fvec(M)
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
      INTEGER i , iev , ivar , j , nm1
      DOUBLE PRECISION c13 , c14 , c29 , c45 , div , dx , eight , five ,&
                     & one , prod , sum , s1 , s2 , temp , ten , ti ,   &
                     & tmp1 , tmp2 , tmp3 , tmp4 , tpi , two , zero ,   &
                     & zp25 , zp5
      DOUBLE PRECISION v(11) , y1(15) , y2(11) , y3(16) , y4(33) ,      &
                     & y5(65)
      DOUBLE PRECISION DFLOAT
      DATA zero , zp25 , zp5 , one , two , five , eight , ten , c13 ,   &
         & c14 , c29 , c45/0.0D0 , 2.5D-1 , 5.0D-1 , 1.0D0 , 2.0D0 ,    &
         & 5.0D0 , 8.0D0 , 1.0D1 , 1.3D1 , 1.4D1 , 2.9D1 , 4.5D1/
      DATA v(1) , v(2) , v(3) , v(4) , v(5) , v(6) , v(7) , v(8) ,      &
         & v(9) , v(10) , v(11)/4.0D0 , 2.0D0 , 1.0D0 , 5.0D-1 ,        &
         & 2.5D-1 , 1.67D-1 , 1.25D-1 , 1.0D-1 , 8.33D-2 , 7.14D-2 ,    &
         & 6.25D-2/
      DATA y1(1) , y1(2) , y1(3) , y1(4) , y1(5) , y1(6) , y1(7) ,      &
         & y1(8) , y1(9) , y1(10) , y1(11) , y1(12) , y1(13) , y1(14) , &
         & y1(15)/1.4D-1 , 1.8D-1 , 2.2D-1 , 2.5D-1 , 2.9D-1 , 3.2D-1 , &
         & 3.5D-1 , 3.9D-1 , 3.7D-1 , 5.8D-1 , 7.3D-1 , 9.6D-1 ,        &
         & 1.34D0 , 2.1D0 , 4.39D0/
      DATA y2(1) , y2(2) , y2(3) , y2(4) , y2(5) , y2(6) , y2(7) ,      &
         & y2(8) , y2(9) , y2(10) , y2(11)/1.957D-1 , 1.947D-1 ,        &
         & 1.735D-1 , 1.6D-1 , 8.44D-2 , 6.27D-2 , 4.56D-2 , 3.42D-2 ,  &
         & 3.23D-2 , 2.35D-2 , 2.46D-2/
      DATA y3(1) , y3(2) , y3(3) , y3(4) , y3(5) , y3(6) , y3(7) ,      &
         & y3(8) , y3(9) , y3(10) , y3(11) , y3(12) , y3(13) , y3(14) , &
         & y3(15) , y3(16)/3.478D4 , 2.861D4 , 2.365D4 , 1.963D4 ,      &
         & 1.637D4 , 1.372D4 , 1.154D4 , 9.744D3 , 8.261D3 , 7.03D3 ,   &
         & 6.005D3 , 5.147D3 , 4.427D3 , 3.82D3 , 3.307D3 , 2.872D3/
      DATA y4(1) , y4(2) , y4(3) , y4(4) , y4(5) , y4(6) , y4(7) ,      &
         & y4(8) , y4(9) , y4(10) , y4(11) , y4(12) , y4(13) , y4(14) , &
         & y4(15) , y4(16) , y4(17) , y4(18) , y4(19) , y4(20) , y4(21) &
         & , y4(22) , y4(23) , y4(24) , y4(25) , y4(26) , y4(27) ,      &
         & y4(28) , y4(29) , y4(30) , y4(31) , y4(32) , y4(33)/8.44D-1 ,&
         & 9.08D-1 , 9.32D-1 , 9.36D-1 , 9.25D-1 , 9.08D-1 , 8.81D-1 ,  &
         & 8.5D-1 , 8.18D-1 , 7.84D-1 , 7.51D-1 , 7.18D-1 , 6.85D-1 ,   &
         & 6.58D-1 , 6.28D-1 , 6.03D-1 , 5.8D-1 , 5.58D-1 , 5.38D-1 ,   &
         & 5.22D-1 , 5.06D-1 , 4.9D-1 , 4.78D-1 , 4.67D-1 , 4.57D-1 ,   &
         & 4.48D-1 , 4.38D-1 , 4.31D-1 , 4.24D-1 , 4.2D-1 , 4.14D-1 ,   &
         & 4.11D-1 , 4.06D-1/
      DATA y5(1) , y5(2) , y5(3) , y5(4) , y5(5) , y5(6) , y5(7) ,      &
         & y5(8) , y5(9) , y5(10) , y5(11) , y5(12) , y5(13) , y5(14) , &
         & y5(15) , y5(16) , y5(17) , y5(18) , y5(19) , y5(20) , y5(21) &
         & , y5(22) , y5(23) , y5(24) , y5(25) , y5(26) , y5(27) ,      &
         & y5(28) , y5(29) , y5(30) , y5(31) , y5(32) , y5(33) , y5(34) &
         & , y5(35) , y5(36) , y5(37) , y5(38) , y5(39) , y5(40) ,      &
         & y5(41) , y5(42) , y5(43) , y5(44) , y5(45) , y5(46) , y5(47) &
         & , y5(48) , y5(49) , y5(50) , y5(51) , y5(52) , y5(53) ,      &
         & y5(54) , y5(55) , y5(56) , y5(57) , y5(58) , y5(59) , y5(60) &
         & , y5(61) , y5(62) , y5(63) , y5(64) , y5(65)/1.366D0 ,       &
         & 1.191D0 , 1.112D0 , 1.013D0 , 9.91D-1 , 8.85D-1 , 8.31D-1 ,  &
         & 8.47D-1 , 7.86D-1 , 7.25D-1 , 7.46D-1 , 6.79D-1 , 6.08D-1 ,  &
         & 6.55D-1 , 6.16D-1 , 6.06D-1 , 6.02D-1 , 6.26D-1 , 6.51D-1 ,  &
         & 7.24D-1 , 6.49D-1 , 6.49D-1 , 6.94D-1 , 6.44D-1 , 6.24D-1 ,  &
         & 6.61D-1 , 6.12D-1 , 5.58D-1 , 5.33D-1 , 4.95D-1 , 5.0D-1 ,   &
         & 4.23D-1 , 3.95D-1 , 3.75D-1 , 3.72D-1 , 3.91D-1 , 3.96D-1 ,  &
         & 4.05D-1 , 4.28D-1 , 4.29D-1 , 5.23D-1 , 5.62D-1 , 6.07D-1 ,  &
         & 6.53D-1 , 6.72D-1 , 7.08D-1 , 6.33D-1 , 6.68D-1 , 6.45D-1 ,  &
         & 6.32D-1 , 5.91D-1 , 5.59D-1 , 5.97D-1 , 6.25D-1 , 7.39D-1 ,  &
         & 7.1D-1 , 7.29D-1 , 7.2D-1 , 6.36D-1 , 5.81D-1 , 4.28D-1 ,    &
         & 2.92D-1 , 1.62D-1 , 9.8D-2 , 5.4D-2/
      DFLOAT(ivar) = ivar
!
!     FUNCTION ROUTINE SELECTOR.
!
      SELECT CASE (Nprob)
      CASE (2)
!
!     LINEAR FUNCTION - RANK 1.
!
         sum = zero
         DO j = 1 , N
            sum = sum + DFLOAT(j)*X(j)
         ENDDO
         DO i = 1 , M
            Fvec(i) = DFLOAT(i)*sum - one
         ENDDO
      CASE (3)
!
!     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
!
         sum = zero
         nm1 = N - 1
         IF ( nm1>=2 ) THEN
            DO j = 2 , nm1
               sum = sum + DFLOAT(j)*X(j)
            ENDDO
         ENDIF
         DO i = 1 , M
            Fvec(i) = DFLOAT(i-1)*sum - one
         ENDDO
         Fvec(M) = -one
      CASE (4)
!
!     ROSENBROCK FUNCTION.
!
         Fvec(1) = ten*(X(2)-X(1)**2)
         Fvec(2) = one - X(1)
      CASE (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*DATAN(one)
         tmp1 = DSIGN(zp25,X(2))
         IF ( X(1)>zero ) tmp1 = DATAN(X(2)/X(1))/tpi
         IF ( X(1)<zero ) tmp1 = DATAN(X(2)/X(1))/tpi + zp5
         tmp2 = DSQRT(X(1)**2+X(2)**2)
         Fvec(1) = ten*(X(3)-ten*tmp1)
         Fvec(2) = ten*(tmp2-one)
         Fvec(3) = X(3)
      CASE (6)
!
!     POWELL SINGULAR FUNCTION.
!
         Fvec(1) = X(1) + ten*X(2)
         Fvec(2) = DSQRT(five)*(X(3)-X(4))
         Fvec(3) = (X(2)-two*X(3))**2
         Fvec(4) = DSQRT(ten)*(X(1)-X(4))**2
      CASE (7)
!
!     FREUDENSTEIN AND ROTH FUNCTION.
!
         Fvec(1) = -c13 + X(1) + ((five-X(2))*X(2)-two)*X(2)
         Fvec(2) = -c29 + X(1) + ((one+X(2))*X(2)-c14)*X(2)
      CASE (8)
!
!     BARD FUNCTION.
!
         DO i = 1 , 15
            tmp1 = DFLOAT(i)
            tmp2 = DFLOAT(16-i)
            tmp3 = tmp1
            IF ( i>8 ) tmp3 = tmp2
            Fvec(i) = y1(i) - (X(1)+tmp1/(X(2)*tmp2+X(3)*tmp3))
         ENDDO
      CASE (9)
!
!     KOWALIK AND OSBORNE FUNCTION.
!
         DO i = 1 , 11
            tmp1 = v(i)*(v(i)+X(2))
            tmp2 = v(i)*(v(i)+X(3)) + X(4)
            Fvec(i) = y2(i) - X(1)*tmp1/tmp2
         ENDDO
      CASE (10)
!
!     MEYER FUNCTION.
!
         DO i = 1 , 16
            temp = five*DFLOAT(i) + c45 + X(3)
            tmp1 = X(2)/temp
            tmp2 = DEXP(tmp1)
            Fvec(i) = X(1)*tmp2 - y3(i)
         ENDDO
      CASE (11)
!
!     WATSON FUNCTION.
!
         DO i = 1 , 29
            div = DFLOAT(i)/c29
            s1 = zero
            dx = one
            DO j = 2 , N
               s1 = s1 + DFLOAT(j-1)*dx*X(j)
               dx = div*dx
            ENDDO
            s2 = zero
            dx = one
            DO j = 1 , N
               s2 = s2 + dx*X(j)
               dx = div*dx
            ENDDO
            Fvec(i) = s1 - s2**2 - one
         ENDDO
         Fvec(30) = X(1)
         Fvec(31) = X(2) - X(1)**2 - one
      CASE (12)
!
!     BOX 3-DIMENSIONAL FUNCTION.
!
         DO i = 1 , M
            temp = DFLOAT(i)
            tmp1 = temp/ten
            Fvec(i) = DEXP(-tmp1*X(1)) - DEXP(-tmp1*X(2))               &
                    & + (DEXP(-temp)-DEXP(-tmp1))*X(3)
         ENDDO
      CASE (13)
!
!     JENNRICH AND SAMPSON FUNCTION.
!
         DO i = 1 , M
            temp = DFLOAT(i)
            Fvec(i) = two + two*temp - DEXP(temp*X(1)) - DEXP(temp*X(2))
         ENDDO
      CASE (14)
!
!     BROWN AND DENNIS FUNCTION.
!
         DO i = 1 , M
            temp = DFLOAT(i)/five
            tmp1 = X(1) + temp*X(2) - DEXP(temp)
            tmp2 = X(3) + DSIN(temp)*X(4) - DCOS(temp)
            Fvec(i) = tmp1**2 + tmp2**2
         ENDDO
      CASE (15)
!
!     CHEBYQUAD FUNCTION.
!
         DO i = 1 , M
            Fvec(i) = zero
         ENDDO
         DO j = 1 , N
            tmp1 = one
            tmp2 = two*X(j) - one
            temp = two*tmp2
            DO i = 1 , M
               Fvec(i) = Fvec(i) + tmp2
               ti = temp*tmp2 - tmp1
               tmp1 = tmp2
               tmp2 = ti
            ENDDO
         ENDDO
         dx = one/DFLOAT(N)
         iev = -1
         DO i = 1 , M
            Fvec(i) = dx*Fvec(i)
            IF ( iev>0 ) Fvec(i) = Fvec(i) + one/(DFLOAT(i)**2-one)
            iev = -iev
         ENDDO
      CASE (16)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         sum = -DFLOAT(N+1)
         prod = one
         DO j = 1 , N
            sum = sum + X(j)
            prod = X(j)*prod
         ENDDO
         DO i = 1 , N
            Fvec(i) = X(i) + sum
         ENDDO
         Fvec(N) = prod - one
      CASE (17)
!
!     OSBORNE 1 FUNCTION.
!
         DO i = 1 , 33
            temp = ten*DFLOAT(i-1)
            tmp1 = DEXP(-X(4)*temp)
            tmp2 = DEXP(-X(5)*temp)
            Fvec(i) = y4(i) - (X(1)+X(2)*tmp1+X(3)*tmp2)
         ENDDO
      CASE (18)
!
!     OSBORNE 2 FUNCTION.
!
         DO i = 1 , 65
            temp = DFLOAT(i-1)/ten
            tmp1 = DEXP(-X(5)*temp)
            tmp2 = DEXP(-X(6)*(temp-X(9))**2)
            tmp3 = DEXP(-X(7)*(temp-X(10))**2)
            tmp4 = DEXP(-X(8)*(temp-X(11))**2)
            Fvec(i) = y5(i) - (X(1)*tmp1+X(2)*tmp2+X(3)*tmp3+X(4)*tmp4)
         ENDDO
      CASE DEFAULT
!
!     LINEAR FUNCTION - FULL RANK.
!
         sum = zero
         DO j = 1 , N
            sum = sum + X(j)
         ENDDO
         temp = two*sum/DFLOAT(M) + one
         DO i = 1 , M
            Fvec(i) = -temp
            IF ( i<=N ) Fvec(i) = Fvec(i) + X(i)
         ENDDO
      END SELECT
!
!     LAST CARD OF SUBROUTINE SSQFCN.
!
      END
