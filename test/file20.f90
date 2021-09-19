!*==AA0001.spg  processed by SPAG 6.72Dc at 04:43 on 19 Sep 2021
      PROGRAM TEST
      IMPLICIT NONE
!*--AA00013
!     **********
!
!     THIS PROGRAM TESTS THE ABILITY OF CHKDER TO DETECT
!     INCONSISTENCIES BETWEEN FUNCTIONS AND THEIR FIRST DERIVATIVES.
!     FOURTEEN TEST FUNCTION VECTORS AND JACOBIANS ARE USED. ELEVEN OF
!     THE TESTS ARE FALSE(F), I.E. THERE ARE INCONSISTENCIES BETWEEN
!     THE FUNCTION VECTORS AND THE CORRESPONDING JACOBIANS. THREE OF
!     THE TESTS ARE TRUE(T), I.E. THERE ARE NO INCONSISTENCIES. THE
!     DRIVER READS IN DATA, CALLS CHKDER AND PRINTS OUT INFORMATION
!     REQUIRED BY AND RECEIVED FROM CHKDER.
!
!     SUBPROGRAMS CALLED
!
!       MINPACK SUPPLIED ... CHKDER,ERRJAC,INITPT,VECFCN
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER i , ldfjac , lnp , mode , n , nprob , nread , nwrite
      INTEGER na(14) , np(14)
      LOGICAL a(14)
      DOUBLE PRECISION cp , one
      DOUBLE PRECISION diff(10) , err(10) , errmax(14) , errmin(14) ,   &
                     & fjac(10,10) , fvec1(10) , fvec2(10) , x1(10) ,   &
                     & x2(10)
!
!     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.
!     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.
!
      DATA nread , nwrite/5 , 6/
!
      DATA a(1) , a(2) , a(3) , a(4) , a(5) , a(6) , a(7) , a(8) ,      &
         & a(9) , a(10) , a(11) , a(12) , a(13) , a(14)/.FALSE. ,       &
         & .FALSE. , .FALSE. , .TRUE. , .FALSE. , .FALSE. , .FALSE. ,   &
         & .TRUE. , .FALSE. , .FALSE. , .FALSE. , .FALSE. , .TRUE. ,    &
         & .FALSE./
      DATA cp , one/1.23D-1 , 1.0D0/
      ldfjac = 10
 100  READ (nread,99001) nprob , n
99001 FORMAT (2I5)
      IF ( nprob<=0 ) THEN
         WRITE (nwrite,99002) lnp
99002    FORMAT ('1SUMMARY OF ',I3,' TESTS OF CHKDER'/)
         WRITE (nwrite,99003)
99003    FORMAT (' NPROB   N    STATUS     ERRMIN         ERRMAX'/)
         DO i = 1 , lnp
            WRITE (nwrite,99004) np(i) , na(i) , a(i) , errmin(i) ,     &
                               & errmax(i)
99004       FORMAT (I4,I6,6X,L1,3X,2D15.7)
         ENDDO
         STOP
      ELSE
         CALL INITPT(n,x1,nprob,one)
         DO i = 1 , n
            x1(i) = x1(i) + cp
            cp = -cp
         ENDDO
         WRITE (nwrite,99005) nprob , n , a(nprob)
99005    FORMAT (///5X,' PROBLEM',I5,5X,' WITH DIMENSION',I5,2X,' IS  ',&
               & L1)
         mode = 1
         CALL CHKDER(n,n,x1,fvec1,fjac,ldfjac,x2,fvec2,mode,err)
         mode = 2
         CALL VECFCN(n,x1,fvec1,nprob)
         CALL ERRJAC(n,x1,fjac,ldfjac,nprob)
         CALL VECFCN(n,x2,fvec2,nprob)
         CALL CHKDER(n,n,x1,fvec1,fjac,ldfjac,x2,fvec2,mode,err)
         errmin(nprob) = err(1)
         errmax(nprob) = err(1)
         DO i = 1 , n
            diff(i) = fvec2(i) - fvec1(i)
            IF ( errmin(nprob)>err(i) ) errmin(nprob) = err(i)
            IF ( errmax(nprob)<err(i) ) errmax(nprob) = err(i)
         ENDDO
         np(nprob) = nprob
         lnp = nprob
         na(nprob) = n
         WRITE (nwrite,99006) (fvec1(i),i=1,n)
99006    FORMAT (//5X,' FIRST FUNCTION VECTOR   '//(5X,5D15.7))
         WRITE (nwrite,99007) (diff(i),i=1,n)
99007    FORMAT (//5X,' FUNCTION DIFFERENCE VECTOR'//(5X,5D15.7))
         WRITE (nwrite,99008) (err(i),i=1,n)
99008    FORMAT (//5X,' ERROR VECTOR'//(5X,5D15.7))
         GOTO 100
      ENDIF
!
!     LAST CARD OF DERIVATIVE CHECK TEST DRIVER.
!
      END
!*==ERRJAC.spg  processed by SPAG 6.72Dc at 04:43 on 19 Sep 2021
      SUBROUTINE ERRJAC(N,X,Fjac,Ldfjac,Nprob)
      IMPLICIT NONE
!*--ERRJAC97
      INTEGER N , Ldfjac , Nprob
      DOUBLE PRECISION X(N) , Fjac(Ldfjac,N)
!     **********
!
!     SUBROUTINE ERRJAC
!
!     THIS SUBROUTINE IS DERIVED FROM VECJAC WHICH DEFINES THE
!     JACOBIAN MATRICES OF FOURTEEN TEST FUNCTIONS. THE PROBLEM
!     DIMENSIONS ARE AS DESCRIBED IN THE PROLOGUE COMMENTS OF VECFCN.
!     VARIOUS ERRORS ARE DELIBERATELY INTRODUCED TO PROVIDE A TEST
!     FOR CHKDER.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE ERRJAC(N,X,FJAC,LDFJAC,NPROB)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER VARIABLE.
!
!       X IS AN ARRAY OF LENGTH N.
!
!       FJAC IS AN N BY N ARRAY. ON OUTPUT FJAC CONTAINS THE
!         JACOBIAN MATRIX, WITH VARIOUS ERRORS DELIBERATELY
!         INTRODUCED, OF THE NPROB FUNCTION EVALUATED AT X.
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
!     POWELL SINGULAR FUNCTION WITH SIGN REVERSAL AFFECTING ELEMENT
!     (3,3).
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
         Fjac(3,3) = two*Fjac(3,2)
         Fjac(4,1) = two*DSQRT(ten)*(X(1)-X(4))
         Fjac(4,4) = -Fjac(4,1)
      CASE (3)
!
!     POWELL BADLY SCALED FUNCTION WITH THE SIGN OF THE JACOBIAN
!     REVERSED.
!
         Fjac(1,1) = -c1*X(2)
         Fjac(1,2) = -c1*X(1)
         Fjac(2,1) = DEXP(-X(1))
         Fjac(2,2) = DEXP(-X(2))
      CASE (4)
!
!     WOOD FUNCTION WITHOUT ERROR.
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
!     HELICAL VALLEY FUNCTION WITH MULTIPLICATIVE ERROR AFFECTING
!     ELEMENTS (2,1) AND (2,2).
!
         tpi = eight*DATAN(one)
         temp = X(1)**2 + X(2)**2
         temp1 = tpi*temp
         temp2 = DSQRT(temp)
         Fjac(1,1) = hundrd*X(2)/temp1
         Fjac(1,2) = -hundrd*X(1)/temp1
         Fjac(1,3) = ten
         Fjac(2,1) = five*X(1)/temp2
         Fjac(2,2) = five*X(2)/temp2
         Fjac(2,3) = zero
         Fjac(3,1) = zero
         Fjac(3,2) = zero
         Fjac(3,3) = one
      CASE (6)
!
!     WATSON FUNCTION WITH SIGN REVERSALS AFFECTING THE COMPUTATION OF
!     TEMP1.
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
            temp1 = two*(sum1+sum2**2+one)
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
!     CHEBYQUAD FUNCTION WITH JACOBIAN TWICE CORRECT SIZE.
!
         tk = one/DFLOAT(N)
         DO j = 1 , N
            temp1 = one
            temp2 = two*X(j) - one
            temp = two*temp2
            temp3 = zero
            temp4 = two
            DO k = 1 , N
               Fjac(k,j) = two*tk*temp4
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
!     BROWN ALMOST-LINEAR FUNCTION WITHOUT ERROR.
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
!     DISCRETE BOUNDARY VALUE FUNCTION WITH MULTIPLICATIVE ERROR
!     AFFECTING THE JACOBIAN DIAGONAL.
!
         h = one/DFLOAT(N+1)
         DO k = 1 , N
            temp = three*(X(k)+DFLOAT(k)*h+one)**2
            DO j = 1 , N
               Fjac(k,j) = zero
            ENDDO
            Fjac(k,k) = four + temp*h**2
            IF ( k/=1 ) Fjac(k,k-1) = -one
            IF ( k/=N ) Fjac(k,k+1) = -one
         ENDDO
      CASE (10)
!
!     DISCRETE INTEGRAL EQUATION FUNCTION WITH SIGN ERROR AFFECTING
!     THE JACOBIAN DIAGONAL.
!
         h = one/DFLOAT(N+1)
         DO k = 1 , N
            tk = DFLOAT(k)*h
            DO j = 1 , N
               tj = DFLOAT(j)*h
               temp = three*(X(j)+tj+one)**2
               Fjac(k,j) = h*DMIN1(tj*(one-tk),tk*(one-tj))*temp/two
            ENDDO
            Fjac(k,k) = Fjac(k,k) - one
         ENDDO
      CASE (11)
!
!     TRIGONOMETRIC FUNCTION WITH SIGN ERRORS AFFECTING THE
!     OFFDIAGONAL ELEMENTS OF THE JACOBIAN.
!
         DO j = 1 , N
            temp = DSIN(X(j))
            DO k = 1 , N
               Fjac(k,j) = -temp
            ENDDO
            Fjac(j,j) = DFLOAT(j+1)*temp - DCOS(X(j))
         ENDDO
      CASE (12)
!
!     VARIABLY DIMENSIONED FUNCTION WITH OPERATION ERROR AFFECTING
!     THE UPPER TRIANGULAR ELEMENTS OF THE JACOBIAN.
!
         sum = zero
         DO j = 1 , N
            sum = sum + DFLOAT(j)*(X(j)-one)
         ENDDO
         temp = one + six*sum**2
         DO k = 1 , N
            DO j = k , N
               Fjac(k,j) = DFLOAT(k*j)/temp
               Fjac(j,k) = Fjac(k,j)
            ENDDO
            Fjac(k,k) = Fjac(k,k) + one
         ENDDO
      CASE (13)
!
!     BROYDEN TRIDIAGONAL FUNCTION WITHOUT ERROR.
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
!     BROYDEN BANDED FUNCTION WITH SIGN ERROR AFFECTING THE JACOBIAN
!     DIAGONAL.
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
            Fjac(k,k) = two - fiftn*X(k)**2
         ENDDO
      CASE DEFAULT
!
!     ROSENBROCK FUNCTION WITH SIGN REVERSAL AFFECTING ELEMENT (1,1).
!
         Fjac(1,1) = one
         Fjac(1,2) = zero
         Fjac(2,1) = -twenty*X(1)
         Fjac(2,2) = ten
      END SELECT
!
!     LAST CARD OF SUBROUTINE ERRJAC.
!
      END
!*==INITPT.spg  processed by SPAG 6.72Dc at 04:43 on 19 Sep 2021
      SUBROUTINE INITPT(N,X,Nprob,Factor)
      IMPLICIT NONE
!*--INITPT419
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
!*==VECFCN.spg  processed by SPAG 6.72Dc at 04:43 on 19 Sep 2021
      SUBROUTINE VECFCN(N,X,Fvec,Nprob)
      IMPLICIT NONE
!*--VECFCN577
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
