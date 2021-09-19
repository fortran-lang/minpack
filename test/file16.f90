!*==AA0001.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      PROGRAM TEST
      IMPLICIT NONE
!*--AA00013
!     **********
!
!     THIS PROGRAM TESTS CODES FOR THE SOLUTION OF N NONLINEAR
!     EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER AND AN
!     INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA, CALLS THE
!     NONLINEAR EQUATION SOLVER, AND FINALLY PRINTS OUT INFORMATION
!     ON THE PERFORMANCE OF THE SOLVER. THIS IS ONLY A SAMPLE DRIVER,
!     MANY OTHER DRIVERS ARE POSSIBLE. THE INTERFACE SUBROUTINE FCN
!     IS NECESSARY TO TAKE INTO ACCOUNT THE FORMS OF CALLING
!     SEQUENCES USED BY THE FUNCTION AND JACOBIAN SUBROUTINES IN
!     THE VARIOUS NONLINEAR EQUATION SOLVERS.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... FCN
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,HYBRJ1,INITPT,VECFCN
!
!       FORTRAN-SUPPLIED ... DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER i , ic , info , k , ldfjac , lwa , n , NFEv , NJEv ,      &
            & NPRob , nread , ntries , nwrite
      INTEGER na(60) , nf(60) , nj(60) , np(60) , nx(60)
      DOUBLE PRECISION factor , fnorm1 , fnorm2 , one , ten , tol
      DOUBLE PRECISION fnm(60) , fjac(40,40) , fvec(40) , wa(1060) ,    &
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
      ldfjac = 40
      lwa = 1060
      ic = 0
 100  READ (nread,99001) NPRob , n , ntries
99001 FORMAT (3I5)
      IF ( NPRob<=0 ) THEN
         WRITE (nwrite,99002) ic
99002    FORMAT ('1SUMMARY OF ',I3,' CALLS TO HYBRJ1'/)
         WRITE (nwrite,99003)
99003    FORMAT (' NPROB   N    NFEV   NJEV  INFO  FINAL L2 NORM'/)
         DO i = 1 , ic
            WRITE (nwrite,99004) np(i) , na(i) , nf(i) , nj(i) , nx(i) ,&
                               & fnm(i)
99004       FORMAT (I4,I6,2I7,I6,1X,D15.7)
         ENDDO
         STOP
      ELSE
         factor = one
         DO k = 1 , ntries
            ic = ic + 1
            CALL INITPT(n,x,NPRob,factor)
            CALL VECFCN(n,x,fvec,NPRob)
            fnorm1 = ENORM(n,fvec)
            WRITE (nwrite,99005) NPRob , n
99005       FORMAT (////5X,' PROBLEM',I5,5X,' DIMENSION',I5,5X//)
            NFEv = 0
            NJEv = 0
            CALL HYBRJ1(FCN,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
            fnorm2 = ENORM(n,fvec)
            np(ic) = NPRob
            na(ic) = n
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
      SUBROUTINE FCN(N,X,Fvec,Fjac,Ldfjac,Iflag)
      IMPLICIT NONE
!*--FCN102
      INTEGER N , Ldfjac , Iflag
      DOUBLE PRECISION X(N) , Fvec(N) , Fjac(Ldfjac,N)
!     **********
!
!     THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
!     CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
!     EQUATION SOLVER. FCN SHOULD ONLY CALL THE TESTING FUNCTION
!     AND JACOBIAN SUBROUTINES VECFCN AND VECJAC WITH THE
!     APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... VECFCN,VECJAC
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER NPRob , NFEv , NJEv
      COMMON /REFNUM/ NPRob , NFEv , NJEv
      IF ( Iflag==1 ) CALL VECFCN(N,X,Fvec,NPRob)
      IF ( Iflag==2 ) CALL VECJAC(N,X,Fjac,Ldfjac,NPRob)
      IF ( Iflag==1 ) NFEv = NFEv + 1
      IF ( Iflag==2 ) NJEv = NJEv + 1
!
!     LAST CARD OF INTERFACE SUBROUTINE FCN.
!
      END
!*==VECJAC.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE VECJAC(N,X,Fjac,Ldfjac,Nprob)
      IMPLICIT NONE
!*--VECJAC134
      INTEGER N , Ldfjac , Nprob
      DOUBLE PRECISION X(N) , Fjac(Ldfjac,N)
!     **********
!
!     SUBROUTINE VECJAC
!
!     THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF FOURTEEN
!     TEST FUNCTIONS. THE PROBLEM DIMENSIONS ARE AS DESCRIBED
!     IN THE PROLOGUE COMMENTS OF VECFCN.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE VECJAC(N,X,FJAC,LDFJAC,NPROB)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER VARIABLE.
!
!       X IS AN ARRAY OF LENGTH N.
!
!       FJAC IS AN N BY N ARRAY. ON OUTPUT FJAC CONTAINS THE
!         JACOBIAN MATRIX OF THE NPROB FUNCTION EVALUATED AT X.
!
!       LDFJAC IS A POSITIVE INTEGER VARIABLE NOT LESS THAN N
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
!
!       NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 14.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DMIN1,DSIN,DSQRT,
!                            MAX0,MIN0
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER i , ivar , j , k , k1 , k2 , ml , mu
      DOUBLE PRECISION c1 , c3 , c4 , c5 , c6 , c9 , eight , fiftn ,    &
                     & five , four , h , hundrd , one , prod , six ,    &
                     & sum , sum1 , sum2 , temp , temp1 , temp2 ,       &
                     & temp3 , temp4 , ten , three , ti , tj , tk ,     &
                     & tpi , twenty , two , zero
      DOUBLE PRECISION DFLOAT
      DATA zero , one , two , three , four , five , six , eight , ten , &
         & fiftn , twenty , hundrd/0.0D0 , 1.0D0 , 2.0D0 , 3.0D0 ,      &
         & 4.0D0 , 5.0D0 , 6.0D0 , 8.0D0 , 1.0D1 , 1.5D1 , 2.0D1 ,      &
         & 1.0D2/
      DATA c1 , c3 , c4 , c5 , c6 , c9/1.0D4 , 2.0D2 , 2.02D1 , 1.98D1 ,&
         & 1.8D2 , 2.9D1/
      DFLOAT(ivar) = ivar
!
!     JACOBIAN ROUTINE SELECTOR.
!
      SELECT CASE (Nprob)
      CASE (2)
!
!     POWELL SINGULAR FUNCTION.
!
         DO k = 1 , 4
            DO j = 1 , 4
               Fjac(k,j) = zero
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
      CASE (3)
!
!     POWELL BADLY SCALED FUNCTION.
!
         Fjac(1,1) = c1*X(2)
         Fjac(1,2) = c1*X(1)
         Fjac(2,1) = -DEXP(-X(1))
         Fjac(2,2) = -DEXP(-X(2))
      CASE (4)
!
!     WOOD FUNCTION.
!
         DO k = 1 , 4
            DO j = 1 , 4
               Fjac(k,j) = zero
            ENDDO
         ENDDO
         temp1 = X(2) - three*X(1)**2
         temp2 = X(4) - three*X(3)**2
         Fjac(1,1) = -c3*temp1 + one
         Fjac(1,2) = -c3*X(1)
         Fjac(2,1) = -two*c3*X(1)
         Fjac(2,2) = c3 + c4
         Fjac(2,4) = c5
         Fjac(3,3) = -c6*temp2 + one
         Fjac(3,4) = -c6*X(3)
         Fjac(4,2) = c5
         Fjac(4,3) = -two*c6*X(3)
         Fjac(4,4) = c6 + c4
      CASE (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*DATAN(one)
         temp = X(1)**2 + X(2)**2
         temp1 = tpi*temp
         temp2 = DSQRT(temp)
         Fjac(1,1) = hundrd*X(2)/temp1
         Fjac(1,2) = -hundrd*X(1)/temp1
         Fjac(1,3) = ten
         Fjac(2,1) = ten*X(1)/temp2
         Fjac(2,2) = ten*X(2)/temp2
         Fjac(2,3) = zero
         Fjac(3,1) = zero
         Fjac(3,2) = zero
         Fjac(3,3) = one
      CASE (6)
!
!     WATSON FUNCTION.
!
         DO k = 1 , N
            DO j = k , N
               Fjac(k,j) = zero
            ENDDO
         ENDDO
         DO i = 1 , 29
            ti = DFLOAT(i)/c9
            sum1 = zero
            temp = one
            DO j = 2 , N
               sum1 = sum1 + DFLOAT(j-1)*temp*X(j)
               temp = ti*temp
            ENDDO
            sum2 = zero
            temp = one
            DO j = 1 , N
               sum2 = sum2 + temp*X(j)
               temp = ti*temp
            ENDDO
            temp1 = two*(sum1-sum2**2-one)
            temp2 = two*sum2
            temp = ti**2
            tk = one
            DO k = 1 , N
               tj = tk
               DO j = k , N
                  Fjac(k,j) = Fjac(k,j)                                 &
                            & + tj*((DFLOAT(k-1)/ti-temp2)*(DFLOAT(j-1) &
                            & /ti-temp2)-temp1)
                  tj = ti*tj
               ENDDO
               tk = temp*tk
            ENDDO
         ENDDO
         Fjac(1,1) = Fjac(1,1) + six*X(1)**2 - two*X(2) + three
         Fjac(1,2) = Fjac(1,2) - two*X(1)
         Fjac(2,2) = Fjac(2,2) + one
         DO k = 1 , N
            DO j = k , N
               Fjac(j,k) = Fjac(k,j)
            ENDDO
         ENDDO
      CASE (7)
!
!     CHEBYQUAD FUNCTION.
!
         tk = one/DFLOAT(N)
         DO j = 1 , N
            temp1 = one
            temp2 = two*X(j) - one
            temp = two*temp2
            temp3 = zero
            temp4 = two
            DO k = 1 , N
               Fjac(k,j) = tk*temp4
               ti = four*temp2 + temp*temp4 - temp3
               temp3 = temp4
               temp4 = ti
               ti = temp*temp2 - temp1
               temp1 = temp2
               temp2 = ti
            ENDDO
         ENDDO
      CASE (8)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         prod = one
         DO j = 1 , N
            prod = X(j)*prod
            DO k = 1 , N
               Fjac(k,j) = one
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
      CASE (9)
!
!     DISCRETE BOUNDARY VALUE FUNCTION.
!
         h = one/DFLOAT(N+1)
         DO k = 1 , N
            temp = three*(X(k)+DFLOAT(k)*h+one)**2
            DO j = 1 , N
               Fjac(k,j) = zero
            ENDDO
            Fjac(k,k) = two + temp*h**2/two
            IF ( k/=1 ) Fjac(k,k-1) = -one
            IF ( k/=N ) Fjac(k,k+1) = -one
         ENDDO
      CASE (10)
!
!     DISCRETE INTEGRAL EQUATION FUNCTION.
!
         h = one/DFLOAT(N+1)
         DO k = 1 , N
            tk = DFLOAT(k)*h
            DO j = 1 , N
               tj = DFLOAT(j)*h
               temp = three*(X(j)+tj+one)**2
               Fjac(k,j) = h*DMIN1(tj*(one-tk),tk*(one-tj))*temp/two
            ENDDO
            Fjac(k,k) = Fjac(k,k) + one
         ENDDO
      CASE (11)
!
!     TRIGONOMETRIC FUNCTION.
!
         DO j = 1 , N
            temp = DSIN(X(j))
            DO k = 1 , N
               Fjac(k,j) = temp
            ENDDO
            Fjac(j,j) = DFLOAT(j+1)*temp - DCOS(X(j))
         ENDDO
      CASE (12)
!
!     VARIABLY DIMENSIONED FUNCTION.
!
         sum = zero
         DO j = 1 , N
            sum = sum + DFLOAT(j)*(X(j)-one)
         ENDDO
         temp = one + six*sum**2
         DO k = 1 , N
            DO j = k , N
               Fjac(k,j) = DFLOAT(k*j)*temp
               Fjac(j,k) = Fjac(k,j)
            ENDDO
            Fjac(k,k) = Fjac(k,k) + one
         ENDDO
      CASE (13)
!
!     BROYDEN TRIDIAGONAL FUNCTION.
!
         DO k = 1 , N
            DO j = 1 , N
               Fjac(k,j) = zero
            ENDDO
            Fjac(k,k) = three - four*X(k)
            IF ( k/=1 ) Fjac(k,k-1) = -one
            IF ( k/=N ) Fjac(k,k+1) = -two
         ENDDO
      CASE (14)
!
!     BROYDEN BANDED FUNCTION.
!
         ml = 5
         mu = 1
         DO k = 1 , N
            DO j = 1 , N
               Fjac(k,j) = zero
            ENDDO
            k1 = MAX0(1,k-ml)
            k2 = MIN0(k+mu,N)
            DO j = k1 , k2
               IF ( j/=k ) Fjac(k,j) = -(one+two*X(j))
            ENDDO
            Fjac(k,k) = two + fiftn*X(k)**2
         ENDDO
      CASE DEFAULT
!
!     ROSENBROCK FUNCTION.
!
         Fjac(1,1) = -one
         Fjac(1,2) = zero
         Fjac(2,1) = -twenty*X(1)
         Fjac(2,2) = ten
      END SELECT
!
!     LAST CARD OF SUBROUTINE VECJAC.
!
      END
!*==INITPT.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE INITPT(N,X,Nprob,Factor)
      IMPLICIT NONE
!*--INITPT444
      INTEGER N , Nprob
      DOUBLE PRECISION Factor
      DOUBLE PRECISION X(N)
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
      INTEGER ivar , j
      DOUBLE PRECISION c1 , h , half , one , three , tj , zero
      DOUBLE PRECISION DFLOAT
      DATA zero , half , one , three , c1/0.0D0 , 5.0D-1 , 1.0D0 ,      &
         & 3.0D0 , 1.2D0/
      DFLOAT(ivar) = ivar
!
!     SELECTION OF INITIAL POINT.
!
      SELECT CASE (Nprob)
      CASE (2)
!
!     POWELL SINGULAR FUNCTION.
!
         X(1) = three
         X(2) = -one
         X(3) = zero
         X(4) = one
      CASE (3)
!
!     POWELL BADLY SCALED FUNCTION.
!
         X(1) = zero
         X(2) = one
      CASE (4)
!
!     WOOD FUNCTION.
!
         X(1) = -three
         X(2) = -one
         X(3) = -three
         X(4) = -one
      CASE (5)
!
!     HELICAL VALLEY FUNCTION.
!
         X(1) = -one
         X(2) = zero
         X(3) = zero
      CASE (6)
!
!     WATSON FUNCTION.
!
         DO j = 1 , N
            X(j) = zero
         ENDDO
      CASE (7)
!
!     CHEBYQUAD FUNCTION.
!
         h = one/DFLOAT(N+1)
         DO j = 1 , N
            X(j) = DFLOAT(j)*h
         ENDDO
      CASE (8)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         DO j = 1 , N
            X(j) = half
         ENDDO
      CASE (9,10)
!
!     DISCRETE BOUNDARY VALUE AND INTEGRAL EQUATION FUNCTIONS.
!
         h = one/DFLOAT(N+1)
         DO j = 1 , N
            tj = DFLOAT(j)*h
            X(j) = tj*(tj-one)
         ENDDO
      CASE (11)
!
!     TRIGONOMETRIC FUNCTION.
!
         h = one/DFLOAT(N)
         DO j = 1 , N
            X(j) = h
         ENDDO
      CASE (12)
!
!     VARIABLY DIMENSIONED FUNCTION.
!
         h = one/DFLOAT(N)
         DO j = 1 , N
            X(j) = one - DFLOAT(j)*h
         ENDDO
      CASE (13,14)
!
!     BROYDEN TRIDIAGONAL AND BANDED FUNCTIONS.
!
         DO j = 1 , N
            X(j) = -one
         ENDDO
      CASE DEFAULT
!
!     ROSENBROCK FUNCTION.
!
         X(1) = -c1
         X(2) = one
      END SELECT
!
!     COMPUTE MULTIPLE OF INITIAL POINT.
!
      IF ( Factor/=one ) THEN
         IF ( Nprob==6 ) THEN
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
!*==VECFCN.spg  processed by SPAG 6.72Dc at 04:41 on 19 Sep 2021
      SUBROUTINE VECFCN(N,X,Fvec,Nprob)
      IMPLICIT NONE
!*--VECFCN602
      INTEGER N , Nprob
      DOUBLE PRECISION X(N) , Fvec(N)
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
      INTEGER i , iev , ivar , j , k , k1 , k2 , kp1 , ml , mu
      DOUBLE PRECISION c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 ,     &
                     & eight , five , h , one , prod , sum , sum1 ,     &
                     & sum2 , temp , temp1 , temp2 , ten , three , ti , &
                     & tj , tk , tpi , two , zero
      DOUBLE PRECISION DFLOAT
      DATA zero , one , two , three , five , eight , ten/0.0D0 , 1.0D0 ,&
         & 2.0D0 , 3.0D0 , 5.0D0 , 8.0D0 , 1.0D1/
      DATA c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9/1.0D4 , 1.0001D0 ,&
         & 2.0D2 , 2.02D1 , 1.98D1 , 1.8D2 , 2.5D-1 , 5.0D-1 , 2.9D1/
      DFLOAT(ivar) = ivar
!
!     PROBLEM SELECTOR.
!
      SELECT CASE (Nprob)
      CASE (2)
!
!     POWELL SINGULAR FUNCTION.
!
         Fvec(1) = X(1) + ten*X(2)
         Fvec(2) = DSQRT(five)*(X(3)-X(4))
         Fvec(3) = (X(2)-two*X(3))**2
         Fvec(4) = DSQRT(ten)*(X(1)-X(4))**2
      CASE (3)
!
!     POWELL BADLY SCALED FUNCTION.
!
         Fvec(1) = c1*X(1)*X(2) - one
         Fvec(2) = DEXP(-X(1)) + DEXP(-X(2)) - c2
      CASE (4)
!
!     WOOD FUNCTION.
!
         temp1 = X(2) - X(1)**2
         temp2 = X(4) - X(3)**2
         Fvec(1) = -c3*X(1)*temp1 - (one-X(1))
         Fvec(2) = c3*temp1 + c4*(X(2)-one) + c5*(X(4)-one)
         Fvec(3) = -c6*X(3)*temp2 - (one-X(3))
         Fvec(4) = c6*temp2 + c4*(X(4)-one) + c5*(X(2)-one)
      CASE (5)
!
!     HELICAL VALLEY FUNCTION.
!
         tpi = eight*DATAN(one)
         temp1 = DSIGN(c7,X(2))
         IF ( X(1)>zero ) temp1 = DATAN(X(2)/X(1))/tpi
         IF ( X(1)<zero ) temp1 = DATAN(X(2)/X(1))/tpi + c8
         temp2 = DSQRT(X(1)**2+X(2)**2)
         Fvec(1) = ten*(X(3)-ten*temp1)
         Fvec(2) = ten*(temp2-one)
         Fvec(3) = X(3)
      CASE (6)
!
!     WATSON FUNCTION.
!
         DO k = 1 , N
            Fvec(k) = zero
         ENDDO
         DO i = 1 , 29
            ti = DFLOAT(i)/c9
            sum1 = zero
            temp = one
            DO j = 2 , N
               sum1 = sum1 + DFLOAT(j-1)*temp*X(j)
               temp = ti*temp
            ENDDO
            sum2 = zero
            temp = one
            DO j = 1 , N
               sum2 = sum2 + temp*X(j)
               temp = ti*temp
            ENDDO
            temp1 = sum1 - sum2**2 - one
            temp2 = two*ti*sum2
            temp = one/ti
            DO k = 1 , N
               Fvec(k) = Fvec(k) + temp*(DFLOAT(k-1)-temp2)*temp1
               temp = ti*temp
            ENDDO
         ENDDO
         temp = X(2) - X(1)**2 - one
         Fvec(1) = Fvec(1) + X(1)*(one-two*temp)
         Fvec(2) = Fvec(2) + temp
      CASE (7)
!
!     CHEBYQUAD FUNCTION.
!
         DO k = 1 , N
            Fvec(k) = zero
         ENDDO
         DO j = 1 , N
            temp1 = one
            temp2 = two*X(j) - one
            temp = two*temp2
            DO i = 1 , N
               Fvec(i) = Fvec(i) + temp2
               ti = temp*temp2 - temp1
               temp1 = temp2
               temp2 = ti
            ENDDO
         ENDDO
         tk = one/DFLOAT(N)
         iev = -1
         DO k = 1 , N
            Fvec(k) = tk*Fvec(k)
            IF ( iev>0 ) Fvec(k) = Fvec(k) + one/(DFLOAT(k)**2-one)
            iev = -iev
         ENDDO
      CASE (8)
!
!     BROWN ALMOST-LINEAR FUNCTION.
!
         sum = -DFLOAT(N+1)
         prod = one
         DO j = 1 , N
            sum = sum + X(j)
            prod = X(j)*prod
         ENDDO
         DO k = 1 , N
            Fvec(k) = X(k) + sum
         ENDDO
         Fvec(N) = prod - one
      CASE (9)
!
!     DISCRETE BOUNDARY VALUE FUNCTION.
!
         h = one/DFLOAT(N+1)
         DO k = 1 , N
            temp = (X(k)+DFLOAT(k)*h+one)**3
            temp1 = zero
            IF ( k/=1 ) temp1 = X(k-1)
            temp2 = zero
            IF ( k/=N ) temp2 = X(k+1)
            Fvec(k) = two*X(k) - temp1 - temp2 + temp*h**2/two
         ENDDO
      CASE (10)
!
!     DISCRETE INTEGRAL EQUATION FUNCTION.
!
         h = one/DFLOAT(N+1)
         DO k = 1 , N
            tk = DFLOAT(k)*h
            sum1 = zero
            DO j = 1 , k
               tj = DFLOAT(j)*h
               temp = (X(j)+tj+one)**3
               sum1 = sum1 + tj*temp
            ENDDO
            sum2 = zero
            kp1 = k + 1
            IF ( N>=kp1 ) THEN
               DO j = kp1 , N
                  tj = DFLOAT(j)*h
                  temp = (X(j)+tj+one)**3
                  sum2 = sum2 + (one-tj)*temp
               ENDDO
            ENDIF
            Fvec(k) = X(k) + h*((one-tk)*sum1+tk*sum2)/two
         ENDDO
      CASE (11)
!
!     TRIGONOMETRIC FUNCTION.
!
         sum = zero
         DO j = 1 , N
            Fvec(j) = DCOS(X(j))
            sum = sum + Fvec(j)
         ENDDO
         DO k = 1 , N
            Fvec(k) = DFLOAT(N+k) - DSIN(X(k)) - sum - DFLOAT(k)*Fvec(k)
         ENDDO
      CASE (12)
!
!     VARIABLY DIMENSIONED FUNCTION.
!
         sum = zero
         DO j = 1 , N
            sum = sum + DFLOAT(j)*(X(j)-one)
         ENDDO
         temp = sum*(one+two*sum**2)
         DO k = 1 , N
            Fvec(k) = X(k) - one + DFLOAT(k)*temp
         ENDDO
      CASE (13)
!
!     BROYDEN TRIDIAGONAL FUNCTION.
!
         DO k = 1 , N
            temp = (three-two*X(k))*X(k)
            temp1 = zero
            IF ( k/=1 ) temp1 = X(k-1)
            temp2 = zero
            IF ( k/=N ) temp2 = X(k+1)
            Fvec(k) = temp - temp1 - two*temp2 + one
         ENDDO
      CASE (14)
!
!     BROYDEN BANDED FUNCTION.
!
         ml = 5
         mu = 1
         DO k = 1 , N
            k1 = MAX0(1,k-ml)
            k2 = MIN0(k+mu,N)
            temp = zero
            DO j = k1 , k2
               IF ( j/=k ) temp = temp + X(j)*(one+X(j))
            ENDDO
            Fvec(k) = X(k)*(two+five*X(k)**2) + one - temp
         ENDDO
      CASE DEFAULT
!
!     ROSENBROCK FUNCTION.
!
         Fvec(1) = one - X(1)
         Fvec(2) = ten*(X(2)-X(1)**2)
      END SELECT
!
!     LAST CARD OF SUBROUTINE VECFCN.
!
      END
