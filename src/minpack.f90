!*==CHKDER.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE CHKDER(M,N,X,Fvec,Fjac,Ldfjac,Xp,Fvecp,Mode,Err)
      IMPLICIT NONE
!*--CHKDER4
      INTEGER M , N , Ldfjac , Mode
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N) , Xp(N) ,        &
                     & Fvecp(M) , Err(M)
!     **********
!
!     subroutine chkder
!
!     this subroutine checks the gradients of m nonlinear functions
!     in n variables, evaluated at a point x, for consistency with
!     the functions themselves. the user must call chkder twice,
!     first with mode = 1 and then with mode = 2.
!
!     mode = 1. on input, x must contain the point of evaluation.
!               on output, xp is set to a neighboring point.
!
!     mode = 2. on input, fvec must contain the functions and the
!                         rows of fjac must contain the gradients
!                         of the respective functions each evaluated
!                         at x, and fvecp must contain the functions
!                         evaluated at xp.
!               on output, err contains measures of correctness of
!                          the respective gradients.
!
!     the subroutine does not perform reliably if cancellation or
!     rounding errors cause a severe loss of significance in the
!     evaluation of a function. therefore, none of the components
!     of x should be unusually small (in particular, zero) or any
!     other value which may cause loss of significance.
!
!     the subroutine statement is
!
!       subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables.
!
!       x is an input array of length n.
!
!       fvec is an array of length m. on input when mode = 2,
!         fvec must contain the functions evaluated at x.
!
!       fjac is an m by n array. on input when mode = 2,
!         the rows of fjac must contain the gradients of
!         the respective functions evaluated at x.
!
!       ldfjac is a positive integer input parameter not less than m
!         which specifies the leading dimension of the array fjac.
!
!       xp is an array of length n. on output when mode = 1,
!         xp is set to a neighboring point of x.
!
!       fvecp is an array of length m. on input when mode = 2,
!         fvecp must contain the functions evaluated at xp.
!
!       mode is an integer input variable set to 1 on the first call
!         and 2 on the second. other values of mode are equivalent
!         to mode = 1.
!
!       err is an array of length m. on output when mode = 2,
!         err contains measures of correctness of the respective
!         gradients. if there is no severe loss of significance,
!         then if err(i) is 1.0 the i-th gradient is correct,
!         while if err(i) is 0.0 the i-th gradient is incorrect.
!         for values of err between 0.0 and 1.0, the categorization
!         is less certain. in general, a value of err(i) greater
!         than 0.5 indicates that the i-th gradient is probably
!         correct, while a value of err(i) less than 0.5 indicates
!         that the i-th gradient is probably incorrect.
!
!     subprograms called
!
!       minpack supplied ... dpmpar
!
!       fortran supplied ... dabs,dlog10,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j
      DOUBLE PRECISION eps , epsf , epslog , epsmch , factor , one ,    &
                     & temp , zero
      DOUBLE PRECISION DPMPAR
      DATA factor , one , zero/1.0D2 , 1.0D0 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      eps = DSQRT(epsmch)
!
      IF ( Mode==2 ) THEN
!
!        mode = 2.
!
         epsf = factor*epsmch
         epslog = DLOG10(eps)
         DO i = 1 , M
            Err(i) = zero
         ENDDO
         DO j = 1 , N
            temp = DABS(X(j))
            IF ( temp==zero ) temp = one
            DO i = 1 , M
               Err(i) = Err(i) + temp*Fjac(i,j)
            ENDDO
         ENDDO
         DO i = 1 , M
            temp = one
            IF ( Fvec(i)/=zero .AND. Fvecp(i)/=zero .AND.               &
               & DABS(Fvecp(i)-Fvec(i))>=epsf*DABS(Fvec(i)) )           &
               & temp = eps*DABS((Fvecp(i)-Fvec(i))/eps-Err(i))         &
               & /(DABS(Fvec(i))+DABS(Fvecp(i)))
            Err(i) = one
            IF ( temp>epsmch .AND. temp<eps ) Err(i)                    &
               & = (DLOG10(temp)-epslog)/epslog
            IF ( temp>=eps ) Err(i) = zero
         ENDDO
      ELSE
!
!        mode = 1.
!
         DO j = 1 , N
            temp = eps*DABS(X(j))
            IF ( temp==zero ) temp = eps
            Xp(j) = X(j) + temp
         ENDDO
      ENDIF
!
!
!     last card of subroutine chkder.
!
      END
!*==DOGLEG.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE DOGLEG(N,R,Lr,Diag,Qtb,Delta,X,Wa1,Wa2)
      IMPLICIT NONE
!*--DOGLEG146
      INTEGER N , Lr
      DOUBLE PRECISION Delta
      DOUBLE PRECISION R(Lr) , Diag(N) , Qtb(N) , X(N) , Wa1(N) , Wa2(N)
!     **********
!
!     subroutine dogleg
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta, the
!     problem is to determine the convex combination x of the
!     gauss-newton and scaled gradient directions that minimizes
!     (a*x - b) in the least squares sense, subject to the
!     restriction that the euclidean norm of d*x be at most delta.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization of a. that is, if a = q*r, where q has
!     orthogonal columns and r is an upper triangular matrix,
!     then dogleg expects the full upper triangle of r and
!     the first n components of (q transpose)*b.
!
!     the subroutine statement is
!
!       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an input array of length lr which must contain the upper
!         triangular matrix r stored by rows.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       x is an output array of length n which contains the desired
!         convex combination of the gauss-newton direction and the
!         scaled gradient direction.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jj , jp1 , k , l
      DOUBLE PRECISION alpha , bnorm , epsmch , gnorm , one , qnorm ,   &
                     & sgnorm , sum , temp , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , zero/1.0D0 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
!     first, calculate the gauss-newton direction.
!
      jj = (N*(N+1))/2 + 1
      DO k = 1 , N
         j = N - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         IF ( N>=jp1 ) THEN
            DO i = jp1 , N
               sum = sum + R(l)*X(i)
               l = l + 1
            ENDDO
         ENDIF
         temp = R(jj)
         IF ( temp==zero ) THEN
            l = j
            DO i = 1 , j
               temp = DMAX1(temp,DABS(R(l)))
               l = l + N - i
            ENDDO
            temp = epsmch*temp
            IF ( temp==zero ) temp = epsmch
         ENDIF
         X(j) = (Qtb(j)-sum)/temp
      ENDDO
!
!     test whether the gauss-newton direction is acceptable.
!
      DO j = 1 , N
         Wa1(j) = zero
         Wa2(j) = Diag(j)*X(j)
      ENDDO
      qnorm = ENORM(N,Wa2)
      IF ( qnorm>Delta ) THEN
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
         l = 1
         DO j = 1 , N
            temp = Qtb(j)
            DO i = j , N
               Wa1(i) = Wa1(i) + R(l)*temp
               l = l + 1
            ENDDO
            Wa1(j) = Wa1(j)/Diag(j)
         ENDDO
!
!     calculate the norm of the scaled gradient and test for
!     the special case in which the scaled gradient is zero.
!
         gnorm = ENORM(N,Wa1)
         sgnorm = zero
         alpha = Delta/qnorm
         IF ( gnorm/=zero ) THEN
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
            DO j = 1 , N
               Wa1(j) = (Wa1(j)/gnorm)/Diag(j)
            ENDDO
            l = 1
            DO j = 1 , N
               sum = zero
               DO i = j , N
                  sum = sum + R(l)*Wa1(i)
                  l = l + 1
               ENDDO
               Wa2(j) = sum
            ENDDO
            temp = ENORM(N,Wa2)
            sgnorm = (gnorm/temp)/temp
!
!     test whether the scaled gradient direction is acceptable.
!
            alpha = zero
            IF ( sgnorm<Delta ) THEN
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
               bnorm = ENORM(N,Qtb)
               temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/Delta)
               temp = temp - (Delta/qnorm)*(sgnorm/Delta)               &
                    & **2 + DSQRT((temp-(Delta/qnorm))                  &
                    & **2+(one-(Delta/qnorm)**2)*(one-(sgnorm/Delta)**2)&
                    & )
               alpha = ((Delta/qnorm)*(one-(sgnorm/Delta)**2))/temp
            ENDIF
         ENDIF
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
         temp = (one-alpha)*DMIN1(sgnorm,Delta)
         DO j = 1 , N
            X(j) = temp*Wa1(j) + alpha*X(j)
         ENDDO
      ENDIF
!
!     last card of subroutine dogleg.
!
      END
!*==DPMPAR.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      DOUBLE PRECISION FUNCTION DPMPAR(I)
      IMPLICIT NONE
!*--DPMPAR327
      INTEGER I
!     **********
!
!     Function dpmpar
!
!     This function provides double precision machine parameters
!     when the appropriate set of data statements is activated (by
!     removing the c from column 1) and all other data statements are
!     rendered inactive. Most of the parameter values were obtained
!     from the corresponding Bell Laboratories Port Library function.
!
!     The function statement is
!
!       double precision function dpmpar(i)
!
!     where
!
!       i is an integer input variable set to 1, 2, or 3 which
!         selects the desired machine parameter. If the machine has
!         t base b digits and its smallest and largest exponents are
!         emin and emax, respectively, then these parameters are
!
!         dpmpar(1) = b**(1 - t), the machine precision,
!
!         dpmpar(2) = b**(emin - 1), the smallest magnitude,
!
!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
!
!     Argonne National Laboratory. MINPACK Project. November 1996.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
!
!     **********
      INTEGER mcheps(4)
      INTEGER minmag(4)
      INTEGER maxmag(4)
      DOUBLE PRECISION dmach(3)
      EQUIVALENCE (dmach(1),mcheps(1))
      EQUIVALENCE (dmach(2),minmag(1))
      EQUIVALENCE (dmach(3),maxmag(1))
!
!     Machine constants for the IBM 360/370 series,
!     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
!     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
!
!     data mcheps(1),mcheps(2) / z34100000, z00000000 /
!     data minmag(1),minmag(2) / z00100000, z00000000 /
!     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
!
!     Machine constants for the Honeywell 600/6000 series.
!
!     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
!     data minmag(1),minmag(2) / o402400000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
!
!     Machine constants for the CDC 6000/7000 series.
!
!     data mcheps(1) / 15614000000000000000b /
!     data mcheps(2) / 15010000000000000000b /
!
!     data minmag(1) / 00604000000000000000b /
!     data minmag(2) / 00000000000000000000b /
!
!     data maxmag(1) / 37767777777777777777b /
!     data maxmag(2) / 37167777777777777777b /
!
!     Machine constants for the PDP-10 (KA processor).
!
!     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
!     data minmag(1),minmag(2) / "033400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
!
!     Machine constants for the PDP-10 (KI processor).
!
!     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
!     data minmag(1),minmag(2) / "000400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
!
!     Machine constants for the PDP-11.
!
!     data mcheps(1),mcheps(2) /   9472,      0 /
!     data mcheps(3),mcheps(4) /      0,      0 /
!
!     data minmag(1),minmag(2) /    128,      0 /
!     data minmag(3),minmag(4) /      0,      0 /
!
!     data maxmag(1),maxmag(2) /  32767,     -1 /
!     data maxmag(3),maxmag(4) /     -1,     -1 /
!
!     Machine constants for the Burroughs 6700/7700 systems.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o7770000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o7777777777777777 /
!
!     Machine constants for the Burroughs 5700 system.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o0000000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o0007777777777777 /
!
!     Machine constants for the Burroughs 1700 system.
!
!     data mcheps(1) / zcc6800000 /
!     data mcheps(2) / z000000000 /
!
!     data minmag(1) / zc00800000 /
!     data minmag(2) / z000000000 /
!
!     data maxmag(1) / zdffffffff /
!     data maxmag(2) / zfffffffff /
!
!     Machine constants for the Univac 1100 series.
!
!     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
!     data minmag(1),minmag(2) / o000040000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
!
!     Machine constants for the Data General Eclipse S/200.
!
!     Note - it may be appropriate to include the following card -
!     static dmach(3)
!
!     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
!     data mcheps/32020k,3*0/
!
!     Machine constants for the Harris 220.
!
!     data mcheps(1),mcheps(2) / '20000000, '00000334 /
!     data minmag(1),minmag(2) / '20000000, '00000201 /
!     data maxmag(1),maxmag(2) / '37777777, '37777577 /
!
!     Machine constants for the Cray-1.
!
!     data mcheps(1) / 0376424000000000000000b /
!     data mcheps(2) / 0000000000000000000000b /
!
!     data minmag(1) / 0200034000000000000000b /
!     data minmag(2) / 0000000000000000000000b /
!
!     data maxmag(1) / 0577777777777777777777b /
!     data maxmag(2) / 0000007777777777777776b /
!
!     Machine constants for the Prime 400.
!
!     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
!     data minmag(1),minmag(2) / :10000000000, :00000100000 /
!     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
!
!     Machine constants for the VAX-11.
!
!     data mcheps(1),mcheps(2) /   9472,  0 /
!     data minmag(1),minmag(2) /    128,  0 /
!     data maxmag(1),maxmag(2) / -32769, -1 /
!
!     Machine constants for IEEE machines.
!
      DATA dmach(1)/2.22044604926D-16/
      DATA dmach(2)/2.22507385852D-308/
      DATA dmach(3)/1.79769313485D+308/
!
      DPMPAR = dmach(I)
!
!     Last card of function dpmpar.
!
      END
!*==ENORM.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      DOUBLE PRECISION FUNCTION ENORM(N,X)
      IMPLICIT NONE
!*--ENORM506
      INTEGER N
      DOUBLE PRECISION X(N)
!     **********
!
!     function enorm
!
!     given an n-vector x, this function calculates the
!     euclidean norm of x.
!
!     the euclidean norm is computed by accumulating the sum of
!     squares in three different sums. the sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. non-destructive underflows are permitted. underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     the definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. the main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. the constants
!     given here are suitable for every known computer.
!
!     the function statement is
!
!       double precision function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i
      DOUBLE PRECISION agiant , floatn , one , rdwarf , rgiant , s1 ,   &
                     & s2 , s3 , xabs , x1max , x3max , zero
      DATA one , zero , rdwarf , rgiant/1.0D0 , 0.0D0 , 3.834D-20 ,     &
         & 1.304D19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = N
      agiant = rgiant/floatn
      DO i = 1 , N
         xabs = DABS(X(i))
         IF ( xabs>rdwarf .AND. xabs<agiant ) THEN
!
!           sum for intermediate components.
!
            s2 = s2 + xabs**2
         ELSEIF ( xabs<=rdwarf ) THEN
!
!              sum for small components.
!
            IF ( xabs<=x3max ) THEN
               IF ( xabs/=zero ) s3 = s3 + (xabs/x3max)**2
            ELSE
               s3 = one + s3*(x3max/xabs)**2
               x3max = xabs
            ENDIF
!
!              sum for large components.
!
         ELSEIF ( xabs<=x1max ) THEN
            s1 = s1 + (xabs/x1max)**2
         ELSE
            s1 = one + s1*(x1max/xabs)**2
            x1max = xabs
         ENDIF
      ENDDO
!
!     calculation of norm.
!
      IF ( s1/=zero ) THEN
         ENORM = x1max*DSQRT(s1+(s2/x1max)/x1max)
      ELSEIF ( s2==zero ) THEN
         ENORM = x3max*DSQRT(s3)
      ELSE
         IF ( s2>=x3max ) ENORM = DSQRT(s2*(one+(x3max/s2)*(x3max*s3)))
         IF ( s2<x3max ) ENORM = DSQRT(x3max*((s2/x3max)+(x3max*s3)))
      ENDIF
!
!     last card of function enorm.
!
      END
!*==FDJAC1.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE FDJAC1(FCN,N,X,Fvec,Fjac,Ldfjac,Iflag,Ml,Mu,Epsfcn,Wa1,&
                      & Wa2)
      IMPLICIT NONE
!*--FDJAC1604
      INTEGER N , Ldfjac , Iflag , Ml , Mu
      DOUBLE PRECISION Epsfcn
      DOUBLE PRECISION X(N) , Fvec(N) , Fjac(Ldfjac,N) , Wa1(N) , Wa2(N)
!     **********
!
!     subroutine fdjac1
!
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     the subroutine statement is
!
!       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
!                         wa1,wa2)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!         functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac1. see description of fcn.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
!         least n, then the jacobian is considered dense, and wa2 is
!         not referenced.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , k , msum
      DOUBLE PRECISION eps , epsmch , h , temp , zero
      DOUBLE PRECISION DPMPAR
      DATA zero/0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      eps = DSQRT(DMAX1(Epsfcn,epsmch))
      msum = Ml + Mu + 1
      IF ( msum<N ) THEN
!
!        computation of banded approximate jacobian.
!
         DO k = 1 , msum
            DO j = k , N , msum
               Wa2(j) = X(j)
               h = eps*DABS(Wa2(j))
               IF ( h==zero ) h = eps
               X(j) = Wa2(j) + h
            ENDDO
            CALL FCN(N,X,Wa1,Iflag)
            IF ( Iflag<0 ) GOTO 99999
            DO j = k , N , msum
               X(j) = Wa2(j)
               h = eps*DABS(Wa2(j))
               IF ( h==zero ) h = eps
               DO i = 1 , N
                  Fjac(i,j) = zero
                  IF ( i>=j-Mu .AND. i<=j+Ml ) Fjac(i,j)                &
                     & = (Wa1(i)-Fvec(i))/h
               ENDDO
            ENDDO
         ENDDO
      ELSE
!
!        computation of dense approximate jacobian.
!
         DO j = 1 , N
            temp = X(j)
            h = eps*DABS(temp)
            IF ( h==zero ) h = eps
            X(j) = temp + h
            CALL FCN(N,X,Wa1,Iflag)
            IF ( Iflag<0 ) GOTO 99999
            X(j) = temp
            DO i = 1 , N
               Fjac(i,j) = (Wa1(i)-Fvec(i))/h
            ENDDO
         ENDDO
      ENDIF
!
!     last card of subroutine fdjac1.
!
99999 END
!*==FDJAC2.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
 
      SUBROUTINE FDJAC2(FCN,M,N,X,Fvec,Fjac,Ldfjac,Iflag,Epsfcn,Wa)
      IMPLICIT NONE
!*--FDJAC2753
      INTEGER M , N , Ldfjac , Iflag
      DOUBLE PRECISION Epsfcn
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N) , Wa(M)
!     **********
!
!     subroutine fdjac2
!
!     this subroutine computes a forward-difference approximation
!     to the m by n jacobian matrix associated with a specified
!     problem of m functions in n variables.
!
!     the subroutine statement is
!
!       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac2.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an input array of length n.
!
!       fvec is an input array of length m which must contain the
!         functions evaluated at x.
!
!       fjac is an output m by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac2. see description of fcn.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j
      DOUBLE PRECISION eps , epsmch , h , temp , zero
      DOUBLE PRECISION DPMPAR
      DATA zero/0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      eps = DSQRT(DMAX1(Epsfcn,epsmch))
      DO j = 1 , N
         temp = X(j)
         h = eps*DABS(temp)
         IF ( h==zero ) h = eps
         X(j) = temp + h
         CALL FCN(M,N,X,Wa,Iflag)
         IF ( Iflag<0 ) GOTO 99999
         X(j) = temp
         DO i = 1 , M
            Fjac(i,j) = (Wa(i)-Fvec(i))/h
         ENDDO
      ENDDO
!
!     last card of subroutine fdjac2.
!
99999 END
!*==HYBRD.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE HYBRD(FCN,N,X,Fvec,Xtol,Maxfev,Ml,Mu,Epsfcn,Diag,Mode, &
                     & Factor,Nprint,Info,Nfev,Fjac,Ldfjac,R,Lr,Qtf,Wa1,&
                     & Wa2,Wa3,Wa4)
      IMPLICIT NONE
!*--HYBRD863
      INTEGER N , Maxfev , Ml , Mu , Mode , Nprint , Info , Nfev ,      &
            & Ldfjac , Lr
      DOUBLE PRECISION Xtol , Epsfcn , Factor
      DOUBLE PRECISION X(N) , Fvec(N) , Diag(N) , Fjac(Ldfjac,N) ,      &
                     & R(Lr) , Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) ,      &
                     & Wa4(N)
      EXTERNAL FCN
!     **********
!
!     subroutine hybrd
!
!     the purpose of hybrd is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least maxfev
!         by the end of an iteration.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
!                            qform,qrfac,r1mpyq,r1updt
!
!       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , iflag , iter , j , jm1 , l , msum , ncfail , ncsuc ,  &
            & nslow1 , nslow2
      INTEGER iwa(1)
      LOGICAL jeval , sing
      DOUBLE PRECISION actred , delta , epsmch , fnorm , fnorm1 , one , &
                     & pnorm , prered , p1 , p5 , p001 , p0001 , ratio ,&
                     & sum , temp , xnorm , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p1 , p5 , p001 , p0001 , zero/1.0D0 , 1.0D-1 , 5.0D-1 ,&
         & 1.0D-3 , 1.0D-4 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
!
!     check the input parameters for errors.
!
      IF ( N<=0 .OR. Xtol<zero .OR. Maxfev<=0 .OR. Ml<0 .OR. Mu<0 .OR.  &
         & Factor<=zero .OR. Ldfjac<N .OR. Lr<(N*(N+1))/2 ) GOTO 300
      IF ( Mode==2 ) THEN
         DO j = 1 , N
            IF ( Diag(j)<=zero ) GOTO 300
         ENDDO
      ENDIF
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      CALL FCN(N,X,Fvec,iflag)
      Nfev = 1
      IF ( iflag<0 ) GOTO 300
      fnorm = ENORM(N,Fvec)
!
!     determine the number of calls to fcn needed to compute
!     the jacobian matrix.
!
      msum = MIN0(Ml+Mu+1,N)
!
!     initialize iteration counter and monitors.
!
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!
!     beginning of the outer loop.
!
 100  jeval = .TRUE.
!
!        calculate the jacobian matrix.
!
      iflag = 2
      CALL FDJAC1(FCN,N,X,Fvec,Fjac,Ldfjac,iflag,Ml,Mu,Epsfcn,Wa1,Wa2)
      Nfev = Nfev + msum
      IF ( iflag<0 ) GOTO 300
!
!        compute the qr factorization of the jacobian.
!
      CALL QRFAC(N,N,Fjac,Ldfjac,.FALSE.,iwa,1,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      IF ( iter==1 ) THEN
         IF ( Mode/=2 ) THEN
            DO j = 1 , N
               Diag(j) = Wa2(j)
               IF ( Wa2(j)==zero ) Diag(j) = one
            ENDDO
         ENDIF
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         DO j = 1 , N
            Wa3(j) = Diag(j)*X(j)
         ENDDO
         xnorm = ENORM(N,Wa3)
         delta = Factor*xnorm
         IF ( delta==zero ) delta = Factor
      ENDIF
!
!        form (q transpose)*fvec and store in qtf.
!
      DO i = 1 , N
         Qtf(i) = Fvec(i)
      ENDDO
      DO j = 1 , N
         IF ( Fjac(j,j)/=zero ) THEN
            sum = zero
            DO i = j , N
               sum = sum + Fjac(i,j)*Qtf(i)
            ENDDO
            temp = -sum/Fjac(j,j)
            DO i = j , N
               Qtf(i) = Qtf(i) + Fjac(i,j)*temp
            ENDDO
         ENDIF
      ENDDO
!
!        copy the triangular factor of the qr factorization into r.
!
      sing = .FALSE.
      DO j = 1 , N
         l = j
         jm1 = j - 1
         IF ( jm1>=1 ) THEN
            DO i = 1 , jm1
               R(l) = Fjac(i,j)
               l = l + N - i
            ENDDO
         ENDIF
         R(l) = Wa1(j)
         IF ( Wa1(j)==zero ) sing = .TRUE.
      ENDDO
!
!        accumulate the orthogonal factor in fjac.
!
      CALL QFORM(N,N,Fjac,Ldfjac,Wa1)
!
!        rescale if necessary.
!
      IF ( Mode/=2 ) THEN
         DO j = 1 , N
            Diag(j) = DMAX1(Diag(j),Wa2(j))
         ENDDO
      ENDIF
!
!        beginning of the inner loop.
!
!
!           if requested, call fcn to enable printing of iterates.
!
 200  IF ( Nprint>0 ) THEN
         iflag = 0
         IF ( MOD(iter-1,Nprint)==0 ) CALL FCN(N,X,Fvec,iflag)
         IF ( iflag<0 ) GOTO 300
      ENDIF
!
!           determine the direction p.
!
      CALL DOGLEG(N,R,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
      DO j = 1 , N
         Wa1(j) = -Wa1(j)
         Wa2(j) = X(j) + Wa1(j)
         Wa3(j) = Diag(j)*Wa1(j)
      ENDDO
      pnorm = ENORM(N,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      IF ( iter==1 ) delta = DMIN1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      CALL FCN(N,Wa2,Wa4,iflag)
      Nfev = Nfev + 1
      IF ( iflag>=0 ) THEN
         fnorm1 = ENORM(N,Wa4)
!
!           compute the scaled actual reduction.
!
         actred = -one
         IF ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
         l = 1
         DO i = 1 , N
            sum = zero
            DO j = i , N
               sum = sum + R(l)*Wa1(j)
               l = l + 1
            ENDDO
            Wa3(i) = Qtf(i) + sum
         ENDDO
         temp = ENORM(N,Wa3)
         prered = zero
         IF ( temp<fnorm ) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
         ratio = zero
         IF ( prered>zero ) ratio = actred/prered
!
!           update the step bound.
!
         IF ( ratio>=p1 ) THEN
            ncfail = 0
            ncsuc = ncsuc + 1
            IF ( ratio>=p5 .OR. ncsuc>1 ) delta = DMAX1(delta,pnorm/p5)
            IF ( DABS(ratio-one)<=p1 ) delta = pnorm/p5
         ELSE
            ncsuc = 0
            ncfail = ncfail + 1
            delta = p5*delta
         ENDIF
!
!           test for successful iteration.
!
         IF ( ratio>=p0001 ) THEN
!
!           successful iteration. update x, fvec, and their norms.
!
            DO j = 1 , N
               X(j) = Wa2(j)
               Wa2(j) = Diag(j)*X(j)
               Fvec(j) = Wa4(j)
            ENDDO
            xnorm = ENORM(N,Wa2)
            fnorm = fnorm1
            iter = iter + 1
         ENDIF
!
!           determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         IF ( actred>=p001 ) nslow1 = 0
         IF ( jeval ) nslow2 = nslow2 + 1
         IF ( actred>=p1 ) nslow2 = 0
!
!           test for convergence.
!
         IF ( delta<=Xtol*xnorm .OR. fnorm==zero ) Info = 1
         IF ( Info==0 ) THEN
!
!           tests for termination and stringent tolerances.
!
            IF ( Nfev>=Maxfev ) Info = 2
            IF ( p1*DMAX1(p1*delta,pnorm)<=epsmch*xnorm ) Info = 3
            IF ( nslow2==5 ) Info = 4
            IF ( nslow1==10 ) Info = 5
            IF ( Info==0 ) THEN
!
!           criterion for recalculating jacobian approximation
!           by forward differences.
!
               IF ( ncfail==2 ) GOTO 100
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
               DO j = 1 , N
                  sum = zero
                  DO i = 1 , N
                     sum = sum + Fjac(i,j)*Wa4(i)
                  ENDDO
                  Wa2(j) = (sum-Wa3(j))/pnorm
                  Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
                  IF ( ratio>=p0001 ) Qtf(j) = sum
               ENDDO
!
!           compute the qr factorization of the updated jacobian.
!
               CALL R1UPDT(N,N,R,Lr,Wa1,Wa2,Wa3,sing)
               CALL R1MPYQ(N,N,Fjac,Ldfjac,Wa2,Wa3)
               CALL R1MPYQ(1,N,Qtf,1,Wa2,Wa3)
!
!           end of the inner loop.
!
               jeval = .FALSE.
!
!        end of the outer loop.
!
               GOTO 200
            ENDIF
         ENDIF
      ENDIF
!
!     termination, either normal or user imposed.
!
 300  IF ( iflag<0 ) Info = iflag
      iflag = 0
      IF ( Nprint>0 ) CALL FCN(N,X,Fvec,iflag)
!
!     last card of subroutine hybrd.
!
      END
!*==HYBRD1.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE HYBRD1(FCN,N,X,Fvec,Tol,Info,Wa,Lwa)
      IMPLICIT NONE
!*--HYBRD11319
!*** Start of declarations inserted by SPAG
      REAL FCN
!*** End of declarations inserted by SPAG
      INTEGER N , Info , Lwa
      DOUBLE PRECISION Tol
      DOUBLE PRECISION X(N) , Fvec(N) , Wa(Lwa)
      EXTERNAL FCN
!     **********
!
!     subroutine hybrd1
!
!     the purpose of hybrd1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrd. the user
!     must provide a subroutine which calculates the functions.
!     the jacobian is then calculated by a forward-difference
!     approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates that the relative error
!         between x and the solution is at most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   algorithm estimates that the relative error
!                    between x and the solution is at most tol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    200*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (n*(3*n+13))/2.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... hybrd
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER index , j , lr , maxfev , ml , mode , mu , nfev , nprint
      DOUBLE PRECISION epsfcn , factor , one , xtol , zero
      DATA factor , one , zero/1.0D2 , 1.0D0 , 0.0D0/
      Info = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. Tol>=zero .AND. Lwa>=(N*(3*N+13))/2 ) THEN
!
!     call hybrd.
!
         maxfev = 200*(N+1)
         xtol = Tol
         ml = N - 1
         mu = N - 1
         epsfcn = zero
         mode = 2
         DO j = 1 , N
            Wa(j) = one
         ENDDO
         nprint = 0
         lr = (N*(N+1))/2
         index = 6*N + lr
         CALL HYBRD(FCN,N,X,Fvec,xtol,maxfev,ml,mu,epsfcn,Wa(1),mode,   &
                  & factor,nprint,Info,nfev,Wa(index+1),N,Wa(6*N+1),lr, &
                  & Wa(N+1),Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
         IF ( Info==5 ) Info = 4
      ENDIF
!
!     last card of subroutine hybrd1.
!
      END
!*==HYBRJ.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE HYBRJ(FCN,N,X,Fvec,Fjac,Ldfjac,Xtol,Maxfev,Diag,Mode,  &
                     & Factor,Nprint,Info,Nfev,Njev,R,Lr,Qtf,Wa1,Wa2,   &
                     & Wa3,Wa4)
      IMPLICIT NONE
!*--HYBRJ1448
      INTEGER N , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev , Njev ,&
            & Lr
      DOUBLE PRECISION Xtol , Factor
      DOUBLE PRECISION X(N) , Fvec(N) , Fjac(Ldfjac,N) , Diag(N) ,      &
                     & R(Lr) , Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) ,      &
                     & Wa4(N)
!     **********
!
!     subroutine hybrj
!
!     the purpose of hybrj is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
!                        mode,factor,nprint,info,nfev,njev,r,lr,qtf,
!                        wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         double precision x(n),fvec(n),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrj.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. fvec and fjac should not be altered.
!         if nprint is not positive, no special calls of fcn
!         with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn with iflag = 1 has
!                    reached maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dogleg,dpmpar,enorm,
!                            qform,qrfac,r1mpyq,r1updt
!
!       fortran-supplied ... dabs,dmax1,dmin1,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , iflag , iter , j , jm1 , l , ncfail , ncsuc , nslow1 ,&
            & nslow2
      INTEGER iwa(1)
      LOGICAL jeval , sing
      DOUBLE PRECISION actred , delta , epsmch , fnorm , fnorm1 , one , &
                     & pnorm , prered , p1 , p5 , p001 , p0001 , ratio ,&
                     & sum , temp , xnorm , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p1 , p5 , p001 , p0001 , zero/1.0D0 , 1.0D-1 , 5.0D-1 ,&
         & 1.0D-3 , 1.0D-4 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
      Njev = 0
!
!     check the input parameters for errors.
!
      IF ( N<=0 .OR. Ldfjac<N .OR. Xtol<zero .OR. Maxfev<=0 .OR.        &
         & Factor<=zero .OR. Lr<(N*(N+1))/2 ) GOTO 300
      IF ( Mode==2 ) THEN
         DO j = 1 , N
            IF ( Diag(j)<=zero ) GOTO 300
         ENDDO
      ENDIF
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      CALL FCN(N,X,Fvec,Fjac,Ldfjac,iflag)
      Nfev = 1
      IF ( iflag<0 ) GOTO 300
      fnorm = ENORM(N,Fvec)
!
!     initialize iteration counter and monitors.
!
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!
!     beginning of the outer loop.
!
 100  jeval = .TRUE.
!
!        calculate the jacobian matrix.
!
      iflag = 2
      CALL FCN(N,X,Fvec,Fjac,Ldfjac,iflag)
      Njev = Njev + 1
      IF ( iflag<0 ) GOTO 300
!
!        compute the qr factorization of the jacobian.
!
      CALL QRFAC(N,N,Fjac,Ldfjac,.FALSE.,iwa,1,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      IF ( iter==1 ) THEN
         IF ( Mode/=2 ) THEN
            DO j = 1 , N
               Diag(j) = Wa2(j)
               IF ( Wa2(j)==zero ) Diag(j) = one
            ENDDO
         ENDIF
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         DO j = 1 , N
            Wa3(j) = Diag(j)*X(j)
         ENDDO
         xnorm = ENORM(N,Wa3)
         delta = Factor*xnorm
         IF ( delta==zero ) delta = Factor
      ENDIF
!
!        form (q transpose)*fvec and store in qtf.
!
      DO i = 1 , N
         Qtf(i) = Fvec(i)
      ENDDO
      DO j = 1 , N
         IF ( Fjac(j,j)/=zero ) THEN
            sum = zero
            DO i = j , N
               sum = sum + Fjac(i,j)*Qtf(i)
            ENDDO
            temp = -sum/Fjac(j,j)
            DO i = j , N
               Qtf(i) = Qtf(i) + Fjac(i,j)*temp
            ENDDO
         ENDIF
      ENDDO
!
!        copy the triangular factor of the qr factorization into r.
!
      sing = .FALSE.
      DO j = 1 , N
         l = j
         jm1 = j - 1
         IF ( jm1>=1 ) THEN
            DO i = 1 , jm1
               R(l) = Fjac(i,j)
               l = l + N - i
            ENDDO
         ENDIF
         R(l) = Wa1(j)
         IF ( Wa1(j)==zero ) sing = .TRUE.
      ENDDO
!
!        accumulate the orthogonal factor in fjac.
!
      CALL QFORM(N,N,Fjac,Ldfjac,Wa1)
!
!        rescale if necessary.
!
      IF ( Mode/=2 ) THEN
         DO j = 1 , N
            Diag(j) = DMAX1(Diag(j),Wa2(j))
         ENDDO
      ENDIF
!
!        beginning of the inner loop.
!
!
!           if requested, call fcn to enable printing of iterates.
!
 200  IF ( Nprint>0 ) THEN
         iflag = 0
         IF ( MOD(iter-1,Nprint)==0 )                                   &
            & CALL FCN(N,X,Fvec,Fjac,Ldfjac,iflag)
         IF ( iflag<0 ) GOTO 300
      ENDIF
!
!           determine the direction p.
!
      CALL DOGLEG(N,R,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
      DO j = 1 , N
         Wa1(j) = -Wa1(j)
         Wa2(j) = X(j) + Wa1(j)
         Wa3(j) = Diag(j)*Wa1(j)
      ENDDO
      pnorm = ENORM(N,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      IF ( iter==1 ) delta = DMIN1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      CALL FCN(N,Wa2,Wa4,Fjac,Ldfjac,iflag)
      Nfev = Nfev + 1
      IF ( iflag>=0 ) THEN
         fnorm1 = ENORM(N,Wa4)
!
!           compute the scaled actual reduction.
!
         actred = -one
         IF ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
         l = 1
         DO i = 1 , N
            sum = zero
            DO j = i , N
               sum = sum + R(l)*Wa1(j)
               l = l + 1
            ENDDO
            Wa3(i) = Qtf(i) + sum
         ENDDO
         temp = ENORM(N,Wa3)
         prered = zero
         IF ( temp<fnorm ) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
         ratio = zero
         IF ( prered>zero ) ratio = actred/prered
!
!           update the step bound.
!
         IF ( ratio>=p1 ) THEN
            ncfail = 0
            ncsuc = ncsuc + 1
            IF ( ratio>=p5 .OR. ncsuc>1 ) delta = DMAX1(delta,pnorm/p5)
            IF ( DABS(ratio-one)<=p1 ) delta = pnorm/p5
         ELSE
            ncsuc = 0
            ncfail = ncfail + 1
            delta = p5*delta
         ENDIF
!
!           test for successful iteration.
!
         IF ( ratio>=p0001 ) THEN
!
!           successful iteration. update x, fvec, and their norms.
!
            DO j = 1 , N
               X(j) = Wa2(j)
               Wa2(j) = Diag(j)*X(j)
               Fvec(j) = Wa4(j)
            ENDDO
            xnorm = ENORM(N,Wa2)
            fnorm = fnorm1
            iter = iter + 1
         ENDIF
!
!           determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         IF ( actred>=p001 ) nslow1 = 0
         IF ( jeval ) nslow2 = nslow2 + 1
         IF ( actred>=p1 ) nslow2 = 0
!
!           test for convergence.
!
         IF ( delta<=Xtol*xnorm .OR. fnorm==zero ) Info = 1
         IF ( Info==0 ) THEN
!
!           tests for termination and stringent tolerances.
!
            IF ( Nfev>=Maxfev ) Info = 2
            IF ( p1*DMAX1(p1*delta,pnorm)<=epsmch*xnorm ) Info = 3
            IF ( nslow2==5 ) Info = 4
            IF ( nslow1==10 ) Info = 5
            IF ( Info==0 ) THEN
!
!           criterion for recalculating jacobian.
!
               IF ( ncfail==2 ) GOTO 100
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
               DO j = 1 , N
                  sum = zero
                  DO i = 1 , N
                     sum = sum + Fjac(i,j)*Wa4(i)
                  ENDDO
                  Wa2(j) = (sum-Wa3(j))/pnorm
                  Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
                  IF ( ratio>=p0001 ) Qtf(j) = sum
               ENDDO
!
!           compute the qr factorization of the updated jacobian.
!
               CALL R1UPDT(N,N,R,Lr,Wa1,Wa2,Wa3,sing)
               CALL R1MPYQ(N,N,Fjac,Ldfjac,Wa2,Wa3)
               CALL R1MPYQ(1,N,Qtf,1,Wa2,Wa3)
!
!           end of the inner loop.
!
               jeval = .FALSE.
!
!        end of the outer loop.
!
               GOTO 200
            ENDIF
         ENDIF
      ENDIF
!
!     termination, either normal or user imposed.
!
 300  IF ( iflag<0 ) Info = iflag
      iflag = 0
      IF ( Nprint>0 ) CALL FCN(N,X,Fvec,Fjac,Ldfjac,iflag)
!
!     last card of subroutine hybrj.
!
      END
!*==HYBRJ1.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE HYBRJ1(FCN,N,X,Fvec,Fjac,Ldfjac,Tol,Info,Wa,Lwa)
      IMPLICIT NONE
!*--HYBRJ11886
!*** Start of declarations inserted by SPAG
      REAL FCN
!*** End of declarations inserted by SPAG
      INTEGER N , Ldfjac , Info , Lwa
      DOUBLE PRECISION Tol
      DOUBLE PRECISION X(N) , Fvec(N) , Fjac(Ldfjac,N) , Wa(Lwa)
      EXTERNAL FCN
!     **********
!
!     subroutine hybrj1
!
!     the purpose of hybrj1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrj. the user
!     must provide a subroutine which calculates the functions
!     and the jacobian.
!
!     the subroutine statement is
!
!       subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         double precision x(n),fvec(n),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrj1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates that the relative error
!         between x and the solution is at most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   algorithm estimates that the relative error
!                    between x and the solution is at most tol.
!
!         info = 2   number of calls to fcn with iflag = 1 has
!                    reached 100*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (n*(n+13))/2.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... hybrj
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER j , lr , maxfev , mode , nfev , njev , nprint
      DOUBLE PRECISION factor , one , xtol , zero
      DATA factor , one , zero/1.0D2 , 1.0D0 , 0.0D0/
      Info = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. Ldfjac>=N .AND. Tol>=zero .AND. Lwa>=(N*(N+13))/2 )&
         & THEN
!
!     call hybrj.
!
         maxfev = 100*(N+1)
         xtol = Tol
         mode = 2
         DO j = 1 , N
            Wa(j) = one
         ENDDO
         nprint = 0
         lr = (N*(N+1))/2
         CALL HYBRJ(FCN,N,X,Fvec,Fjac,Ldfjac,xtol,maxfev,Wa(1),mode,    &
                  & factor,nprint,Info,nfev,njev,Wa(6*N+1),lr,Wa(N+1),  &
                  & Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
         IF ( Info==5 ) Info = 4
      ENDIF
!
!     last card of subroutine hybrj1.
!
      END
!*==LMDER.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMDER(FCN,M,N,X,Fvec,Fjac,Ldfjac,Ftol,Xtol,Gtol,Maxfev,&
                     & Diag,Mode,Factor,Nprint,Info,Nfev,Njev,Ipvt,Qtf, &
                     & Wa1,Wa2,Wa3,Wa4)
      IMPLICIT NONE
!*--LMDER2020
      INTEGER M , N , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev ,   &
            & Njev
      INTEGER Ipvt(N)
      DOUBLE PRECISION Ftol , Xtol , Gtol , Factor
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N) , Diag(N) ,      &
                     & Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) , Wa4(M)
!     **********
!
!     subroutine lmder
!
!     the purpose of lmder is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                        maxfev,diag,mode,factor,nprint,info,nfev,
!                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         double precision x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmder.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.).100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x, fvec, and fjac
!         available for printing. fvec and fjac should not be
!         altered. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , iflag , iter , j , l
      DOUBLE PRECISION actred , delta , dirder , epsmch , fnorm ,       &
                     & fnorm1 , gnorm , one , par , pnorm , prered ,    &
                     & p1 , p5 , p25 , p75 , p0001 , ratio , sum ,      &
                     & temp , temp1 , temp2 , xnorm , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p1 , p5 , p25 , p75 , p0001 , zero/1.0D0 , 1.0D-1 ,    &
         & 5.0D-1 , 2.5D-1 , 7.5D-1 , 1.0D-4 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
      Njev = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. M>=N .AND. Ldfjac>=M .AND. Ftol>=zero .AND.        &
         & Xtol>=zero .AND. Gtol>=zero .AND. Maxfev>0 .AND.             &
         & Factor>zero ) THEN
         IF ( Mode==2 ) THEN
            DO j = 1 , N
               IF ( Diag(j)<=zero ) GOTO 100
            ENDDO
         ENDIF
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
         iflag = 1
         CALL FCN(M,N,X,Fvec,Fjac,Ldfjac,iflag)
         Nfev = 1
         IF ( iflag>=0 ) THEN
            fnorm = ENORM(M,Fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
            par = zero
            iter = 1
!
!     beginning of the outer loop.
!
!
!        calculate the jacobian matrix.
!
 20         iflag = 2
            CALL FCN(M,N,X,Fvec,Fjac,Ldfjac,iflag)
            Njev = Njev + 1
            IF ( iflag>=0 ) THEN
!
!        if requested, call fcn to enable printing of iterates.
!
               IF ( Nprint>0 ) THEN
                  iflag = 0
                  IF ( MOD(iter-1,Nprint)==0 )                          &
                     & CALL FCN(M,N,X,Fvec,Fjac,Ldfjac,iflag)
                  IF ( iflag<0 ) GOTO 100
               ENDIF
!
!        compute the qr factorization of the jacobian.
!
               CALL QRFAC(M,N,Fjac,Ldfjac,.TRUE.,Ipvt,N,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
               IF ( iter==1 ) THEN
                  IF ( Mode/=2 ) THEN
                     DO j = 1 , N
                        Diag(j) = Wa2(j)
                        IF ( Wa2(j)==zero ) Diag(j) = one
                     ENDDO
                  ENDIF
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
                  DO j = 1 , N
                     Wa3(j) = Diag(j)*X(j)
                  ENDDO
                  xnorm = ENORM(N,Wa3)
                  delta = Factor*xnorm
                  IF ( delta==zero ) delta = Factor
               ENDIF
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
               DO i = 1 , M
                  Wa4(i) = Fvec(i)
               ENDDO
               DO j = 1 , N
                  IF ( Fjac(j,j)/=zero ) THEN
                     sum = zero
                     DO i = j , M
                        sum = sum + Fjac(i,j)*Wa4(i)
                     ENDDO
                     temp = -sum/Fjac(j,j)
                     DO i = j , M
                        Wa4(i) = Wa4(i) + Fjac(i,j)*temp
                     ENDDO
                  ENDIF
                  Fjac(j,j) = Wa1(j)
                  Qtf(j) = Wa4(j)
               ENDDO
!
!        compute the norm of the scaled gradient.
!
               gnorm = zero
               IF ( fnorm/=zero ) THEN
                  DO j = 1 , N
                     l = Ipvt(j)
                     IF ( Wa2(l)/=zero ) THEN
                        sum = zero
                        DO i = 1 , j
                           sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
                        ENDDO
                        gnorm = DMAX1(gnorm,DABS(sum/Wa2(l)))
                     ENDIF
                  ENDDO
               ENDIF
!
!        test for convergence of the gradient norm.
!
               IF ( gnorm<=Gtol ) Info = 4
               IF ( Info==0 ) THEN
!
!        rescale if necessary.
!
                  IF ( Mode/=2 ) THEN
                     DO j = 1 , N
                        Diag(j) = DMAX1(Diag(j),Wa2(j))
                     ENDDO
                  ENDIF
!
!        beginning of the inner loop.
!
!
!           determine the levenberg-marquardt parameter.
!
 25               CALL LMPAR(N,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1, &
                           & Wa2,Wa3,Wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
                  DO j = 1 , N
                     Wa1(j) = -Wa1(j)
                     Wa2(j) = X(j) + Wa1(j)
                     Wa3(j) = Diag(j)*Wa1(j)
                  ENDDO
                  pnorm = ENORM(N,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
                  IF ( iter==1 ) delta = DMIN1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
                  iflag = 1
                  CALL FCN(M,N,Wa2,Wa4,Fjac,Ldfjac,iflag)
                  Nfev = Nfev + 1
                  IF ( iflag>=0 ) THEN
                     fnorm1 = ENORM(M,Wa4)
!
!           compute the scaled actual reduction.
!
                     actred = -one
                     IF ( p1*fnorm1<fnorm ) actred = one -              &
                        & (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
                     DO j = 1 , N
                        Wa3(j) = zero
                        l = Ipvt(j)
                        temp = Wa1(l)
                        DO i = 1 , j
                           Wa3(i) = Wa3(i) + Fjac(i,j)*temp
                        ENDDO
                     ENDDO
                     temp1 = ENORM(N,Wa3)/fnorm
                     temp2 = (DSQRT(par)*pnorm)/fnorm
                     prered = temp1**2 + temp2**2/p5
                     dirder = -(temp1**2+temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
                     ratio = zero
                     IF ( prered/=zero ) ratio = actred/prered
!
!           update the step bound.
!
                     IF ( ratio<=p25 ) THEN
                        IF ( actred>=zero ) temp = p5
                        IF ( actred<zero )                              &
                           & temp = p5*dirder/(dirder+p5*actred)
                        IF ( p1*fnorm1>=fnorm .OR. temp<p1 ) temp = p1
                        delta = temp*DMIN1(delta,pnorm/p1)
                        par = par/temp
                     ELSEIF ( par==zero .OR. ratio>=p75 ) THEN
                        delta = pnorm/p5
                        par = p5*par
                     ENDIF
!
!           test for successful iteration.
!
                     IF ( ratio>=p0001 ) THEN
!
!           successful iteration. update x, fvec, and their norms.
!
                        DO j = 1 , N
                           X(j) = Wa2(j)
                           Wa2(j) = Diag(j)*X(j)
                        ENDDO
                        DO i = 1 , M
                           Fvec(i) = Wa4(i)
                        ENDDO
                        xnorm = ENORM(N,Wa2)
                        fnorm = fnorm1
                        iter = iter + 1
                     ENDIF
!
!           tests for convergence.
!
                     IF ( DABS(actred)<=Ftol .AND. prered<=Ftol .AND.   &
                        & p5*ratio<=one ) Info = 1
                     IF ( delta<=Xtol*xnorm ) Info = 2
                     IF ( DABS(actred)<=Ftol .AND. prered<=Ftol .AND.   &
                        & p5*ratio<=one .AND. Info==2 ) Info = 3
                     IF ( Info==0 ) THEN
!
!           tests for termination and stringent tolerances.
!
                        IF ( Nfev>=Maxfev ) Info = 5
                        IF ( DABS(actred)<=epsmch .AND.                 &
                           & prered<=epsmch .AND. p5*ratio<=one )       &
                           & Info = 6
                        IF ( delta<=epsmch*xnorm ) Info = 7
                        IF ( gnorm<=epsmch ) Info = 8
                        IF ( Info==0 ) THEN
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
!
!        end of the outer loop.
!
                           IF ( ratio>=p0001 ) GOTO 20
                           GOTO 25
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
!     termination, either normal or user imposed.
!
 100  IF ( iflag<0 ) Info = iflag
      iflag = 0
      IF ( Nprint>0 ) CALL FCN(M,N,X,Fvec,Fjac,Ldfjac,iflag)
!
!     last card of subroutine lmder.
!
      END
!*==LMDER1.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMDER1(FCN,M,N,X,Fvec,Fjac,Ldfjac,Tol,Info,Ipvt,Wa,Lwa)
      IMPLICIT NONE
!*--LMDER12477
!*** Start of declarations inserted by SPAG
      REAL FCN
!*** End of declarations inserted by SPAG
      INTEGER M , N , Ldfjac , Info , Lwa
      INTEGER Ipvt(N)
      DOUBLE PRECISION Tol
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N) , Wa(Lwa)
      EXTERNAL FCN
!     **********
!
!     subroutine lmder1
!
!     the purpose of lmder1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmder. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
!                         ipvt,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         double precision x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmder1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached 100*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than 5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmder
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER maxfev , mode , nfev , njev , nprint
      DOUBLE PRECISION factor , ftol , gtol , xtol , zero
      DATA factor , zero/1.0D2 , 0.0D0/
      Info = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. M>=N .AND. Ldfjac>=M .AND. Tol>=zero .AND.         &
         & Lwa>=5*N+M ) THEN
!
!     call lmder.
!
         maxfev = 100*(N+1)
         ftol = Tol
         xtol = Tol
         gtol = zero
         mode = 1
         nprint = 0
         CALL LMDER(FCN,M,N,X,Fvec,Fjac,Ldfjac,ftol,xtol,gtol,maxfev,   &
                  & Wa(1),mode,factor,nprint,Info,nfev,njev,Ipvt,Wa(N+1)&
                  & ,Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
         IF ( Info==8 ) Info = 4
      ENDIF
!
!     last card of subroutine lmder1.
!
      END
!*==LMDIF.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMDIF(FCN,M,N,X,Fvec,Ftol,Xtol,Gtol,Maxfev,Epsfcn,Diag,&
                     & Mode,Factor,Nprint,Info,Nfev,Fjac,Ldfjac,Ipvt,   &
                     & Qtf,Wa1,Wa2,Wa3,Wa4)
      IMPLICIT NONE
!*--LMDIF2639
      INTEGER M , N , Maxfev , Mode , Nprint , Info , Nfev , Ldfjac
      INTEGER Ipvt(N)
      DOUBLE PRECISION Ftol , Xtol , Gtol , Epsfcn , Factor
      DOUBLE PRECISION X(N) , Fvec(M) , Diag(N) , Fjac(Ldfjac,N) ,      &
                     & Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) , Wa4(M)
      EXTERNAL FCN
!     **********
!
!     subroutine lmdif
!
!     the purpose of lmdif is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , iflag , iter , j , l
      DOUBLE PRECISION actred , delta , dirder , epsmch , fnorm ,       &
                     & fnorm1 , gnorm , one , par , pnorm , prered ,    &
                     & p1 , p5 , p25 , p75 , p0001 , ratio , sum ,      &
                     & temp , temp1 , temp2 , xnorm , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p1 , p5 , p25 , p75 , p0001 , zero/1.0D0 , 1.0D-1 ,    &
         & 5.0D-1 , 2.5D-1 , 7.5D-1 , 1.0D-4 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. M>=N .AND. Ldfjac>=M .AND. Ftol>=zero .AND.        &
         & Xtol>=zero .AND. Gtol>=zero .AND. Maxfev>0 .AND.             &
         & Factor>zero ) THEN
         IF ( Mode==2 ) THEN
            DO j = 1 , N
               IF ( Diag(j)<=zero ) GOTO 100
            ENDDO
         ENDIF
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
         iflag = 1
         CALL FCN(M,N,X,Fvec,iflag)
         Nfev = 1
         IF ( iflag>=0 ) THEN
            fnorm = ENORM(M,Fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
            par = zero
            iter = 1
!
!     beginning of the outer loop.
!
!
!        calculate the jacobian matrix.
!
 20         iflag = 2
            CALL FDJAC2(FCN,M,N,X,Fvec,Fjac,Ldfjac,iflag,Epsfcn,Wa4)
            Nfev = Nfev + N
            IF ( iflag>=0 ) THEN
!
!        if requested, call fcn to enable printing of iterates.
!
               IF ( Nprint>0 ) THEN
                  iflag = 0
                  IF ( MOD(iter-1,Nprint)==0 )                          &
                     & CALL FCN(M,N,X,Fvec,iflag)
                  IF ( iflag<0 ) GOTO 100
               ENDIF
!
!        compute the qr factorization of the jacobian.
!
               CALL QRFAC(M,N,Fjac,Ldfjac,.TRUE.,Ipvt,N,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
               IF ( iter==1 ) THEN
                  IF ( Mode/=2 ) THEN
                     DO j = 1 , N
                        Diag(j) = Wa2(j)
                        IF ( Wa2(j)==zero ) Diag(j) = one
                     ENDDO
                  ENDIF
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
                  DO j = 1 , N
                     Wa3(j) = Diag(j)*X(j)
                  ENDDO
                  xnorm = ENORM(N,Wa3)
                  delta = Factor*xnorm
                  IF ( delta==zero ) delta = Factor
               ENDIF
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
               DO i = 1 , M
                  Wa4(i) = Fvec(i)
               ENDDO
               DO j = 1 , N
                  IF ( Fjac(j,j)/=zero ) THEN
                     sum = zero
                     DO i = j , M
                        sum = sum + Fjac(i,j)*Wa4(i)
                     ENDDO
                     temp = -sum/Fjac(j,j)
                     DO i = j , M
                        Wa4(i) = Wa4(i) + Fjac(i,j)*temp
                     ENDDO
                  ENDIF
                  Fjac(j,j) = Wa1(j)
                  Qtf(j) = Wa4(j)
               ENDDO
!
!        compute the norm of the scaled gradient.
!
               gnorm = zero
               IF ( fnorm/=zero ) THEN
                  DO j = 1 , N
                     l = Ipvt(j)
                     IF ( Wa2(l)/=zero ) THEN
                        sum = zero
                        DO i = 1 , j
                           sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
                        ENDDO
                        gnorm = DMAX1(gnorm,DABS(sum/Wa2(l)))
                     ENDIF
                  ENDDO
               ENDIF
!
!        test for convergence of the gradient norm.
!
               IF ( gnorm<=Gtol ) Info = 4
               IF ( Info==0 ) THEN
!
!        rescale if necessary.
!
                  IF ( Mode/=2 ) THEN
                     DO j = 1 , N
                        Diag(j) = DMAX1(Diag(j),Wa2(j))
                     ENDDO
                  ENDIF
!
!        beginning of the inner loop.
!
!
!           determine the levenberg-marquardt parameter.
!
 25               CALL LMPAR(N,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1, &
                           & Wa2,Wa3,Wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
                  DO j = 1 , N
                     Wa1(j) = -Wa1(j)
                     Wa2(j) = X(j) + Wa1(j)
                     Wa3(j) = Diag(j)*Wa1(j)
                  ENDDO
                  pnorm = ENORM(N,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
                  IF ( iter==1 ) delta = DMIN1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
                  iflag = 1
                  CALL FCN(M,N,Wa2,Wa4,iflag)
                  Nfev = Nfev + 1
                  IF ( iflag>=0 ) THEN
                     fnorm1 = ENORM(M,Wa4)
!
!           compute the scaled actual reduction.
!
                     actred = -one
                     IF ( p1*fnorm1<fnorm ) actred = one -              &
                        & (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
                     DO j = 1 , N
                        Wa3(j) = zero
                        l = Ipvt(j)
                        temp = Wa1(l)
                        DO i = 1 , j
                           Wa3(i) = Wa3(i) + Fjac(i,j)*temp
                        ENDDO
                     ENDDO
                     temp1 = ENORM(N,Wa3)/fnorm
                     temp2 = (DSQRT(par)*pnorm)/fnorm
                     prered = temp1**2 + temp2**2/p5
                     dirder = -(temp1**2+temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
                     ratio = zero
                     IF ( prered/=zero ) ratio = actred/prered
!
!           update the step bound.
!
                     IF ( ratio<=p25 ) THEN
                        IF ( actred>=zero ) temp = p5
                        IF ( actred<zero )                              &
                           & temp = p5*dirder/(dirder+p5*actred)
                        IF ( p1*fnorm1>=fnorm .OR. temp<p1 ) temp = p1
                        delta = temp*DMIN1(delta,pnorm/p1)
                        par = par/temp
                     ELSEIF ( par==zero .OR. ratio>=p75 ) THEN
                        delta = pnorm/p5
                        par = p5*par
                     ENDIF
!
!           test for successful iteration.
!
                     IF ( ratio>=p0001 ) THEN
!
!           successful iteration. update x, fvec, and their norms.
!
                        DO j = 1 , N
                           X(j) = Wa2(j)
                           Wa2(j) = Diag(j)*X(j)
                        ENDDO
                        DO i = 1 , M
                           Fvec(i) = Wa4(i)
                        ENDDO
                        xnorm = ENORM(N,Wa2)
                        fnorm = fnorm1
                        iter = iter + 1
                     ENDIF
!
!           tests for convergence.
!
                     IF ( DABS(actred)<=Ftol .AND. prered<=Ftol .AND.   &
                        & p5*ratio<=one ) Info = 1
                     IF ( delta<=Xtol*xnorm ) Info = 2
                     IF ( DABS(actred)<=Ftol .AND. prered<=Ftol .AND.   &
                        & p5*ratio<=one .AND. Info==2 ) Info = 3
                     IF ( Info==0 ) THEN
!
!           tests for termination and stringent tolerances.
!
                        IF ( Nfev>=Maxfev ) Info = 5
                        IF ( DABS(actred)<=epsmch .AND.                 &
                           & prered<=epsmch .AND. p5*ratio<=one )       &
                           & Info = 6
                        IF ( delta<=epsmch*xnorm ) Info = 7
                        IF ( gnorm<=epsmch ) Info = 8
                        IF ( Info==0 ) THEN
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
!
!        end of the outer loop.
!
                           IF ( ratio>=p0001 ) GOTO 20
                           GOTO 25
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
!     termination, either normal or user imposed.
!
 100  IF ( iflag<0 ) Info = iflag
      iflag = 0
      IF ( Nprint>0 ) CALL FCN(M,N,X,Fvec,iflag)
!
!     last card of subroutine lmdif.
!
      END
!*==LMDIF1.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMDIF1(FCN,M,N,X,Fvec,Tol,Info,Iwa,Wa,Lwa)
      IMPLICIT NONE
!*--LMDIF13098
!*** Start of declarations inserted by SPAG
      REAL FCN
!*** End of declarations inserted by SPAG
      INTEGER M , N , Info , Lwa
      INTEGER Iwa(N)
      DOUBLE PRECISION Tol
      DOUBLE PRECISION X(N) , Fvec(M) , Wa(Lwa)
      EXTERNAL FCN
!     **********
!
!     subroutine lmdif1
!
!     the purpose of lmdif1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmdif. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded 200*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       iwa is an integer work array of length n.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         m*n+5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmdif
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER maxfev , mode , mp5n , nfev , nprint
      DOUBLE PRECISION epsfcn , factor , ftol , gtol , xtol , zero
      DATA factor , zero/1.0D2 , 0.0D0/
      Info = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. M>=N .AND. Tol>=zero .AND. Lwa>=M*N+5*N+M ) THEN
!
!     call lmdif.
!
         maxfev = 200*(N+1)
         ftol = Tol
         xtol = Tol
         gtol = zero
         epsfcn = zero
         mode = 1
         nprint = 0
         mp5n = M + 5*N
         CALL LMDIF(FCN,M,N,X,Fvec,ftol,xtol,gtol,maxfev,epsfcn,Wa(1),  &
                  & mode,factor,nprint,Info,nfev,Wa(mp5n+1),M,Iwa,      &
                  & Wa(N+1),Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
         IF ( Info==8 ) Info = 4
      ENDIF
!
!     last card of subroutine lmdif1.
!
      END
!*==LMPAR.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMPAR(N,R,Ldr,Ipvt,Diag,Qtb,Delta,Par,X,Sdiag,Wa1,Wa2)
      IMPLICIT NONE
!*--LMPAR3237
      INTEGER N , Ldr
      INTEGER Ipvt(N)
      DOUBLE PRECISION Delta , Par
      DOUBLE PRECISION R(Ldr,N) , Diag(N) , Qtb(N) , X(N) , Sdiag(N) ,  &
                     & Wa1(N) , Wa2(N)
!     **********
!
!     subroutine lmpar
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system
!
!           a*x = b ,     sqrt(par)*d*x = 0 ,
!
!     in the least squares sense, and dxnorm is the euclidean
!     norm of d*x, then either par is zero and
!
!           (dxnorm-delta) .le. 0.1*delta ,
!
!     or par is positive and
!
!           abs(dxnorm-delta) .le. 0.1*delta .
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then lmpar expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. on output
!     lmpar also provides an upper triangular matrix s such that
!
!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .
!
!     s is employed within lmpar and may be of separate interest.
!
!     only a few iterations are generally needed for convergence
!     of the algorithm. if, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.
!
!     the subroutine statement is
!
!       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
!                        wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       par is a nonnegative variable. on input par contains an
!         initial estimate of the levenberg-marquardt parameter.
!         on output par contains the final estimate.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!         for the output par.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm,qrsolv
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , iter , j , jm1 , jp1 , k , l , nsing
      DOUBLE PRECISION dxnorm , dwarf , fp , gnorm , parc , parl ,      &
                     & paru , p1 , p001 , sum , temp , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA p1 , p001 , zero/1.0D-1 , 1.0D-3 , 0.0D0/
!
!     dwarf is the smallest positive magnitude.
!
      dwarf = DPMPAR(2)
!
!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.
!
      nsing = N
      DO j = 1 , N
         Wa1(j) = Qtb(j)
         IF ( R(j,j)==zero .AND. nsing==N ) nsing = j - 1
         IF ( nsing<N ) Wa1(j) = zero
      ENDDO
      IF ( nsing>=1 ) THEN
         DO k = 1 , nsing
            j = nsing - k + 1
            Wa1(j) = Wa1(j)/R(j,j)
            temp = Wa1(j)
            jm1 = j - 1
            IF ( jm1>=1 ) THEN
               DO i = 1 , jm1
                  Wa1(i) = Wa1(i) - R(i,j)*temp
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      DO j = 1 , N
         l = Ipvt(j)
         X(l) = Wa1(j)
      ENDDO
!
!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.
!
      iter = 0
      DO j = 1 , N
         Wa2(j) = Diag(j)*X(j)
      ENDDO
      dxnorm = ENORM(N,Wa2)
      fp = dxnorm - Delta
      IF ( fp<=p1*Delta ) THEN
!
!     termination.
!
         IF ( iter==0 ) Par = zero
      ELSE
!
!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function. otherwise set this bound to zero.
!
         parl = zero
         IF ( nsing>=N ) THEN
            DO j = 1 , N
               l = Ipvt(j)
               Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
            ENDDO
            DO j = 1 , N
               sum = zero
               jm1 = j - 1
               IF ( jm1>=1 ) THEN
                  DO i = 1 , jm1
                     sum = sum + R(i,j)*Wa1(i)
                  ENDDO
               ENDIF
               Wa1(j) = (Wa1(j)-sum)/R(j,j)
            ENDDO
            temp = ENORM(N,Wa1)
            parl = ((fp/Delta)/temp)/temp
         ENDIF
!
!     calculate an upper bound, paru, for the zero of the function.
!
         DO j = 1 , N
            sum = zero
            DO i = 1 , j
               sum = sum + R(i,j)*Qtb(i)
            ENDDO
            l = Ipvt(j)
            Wa1(j) = sum/Diag(l)
         ENDDO
         gnorm = ENORM(N,Wa1)
         paru = gnorm/Delta
         IF ( paru==zero ) paru = dwarf/DMIN1(Delta,p1)
!
!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.
!
         Par = DMAX1(Par,parl)
         Par = DMIN1(Par,paru)
         IF ( Par==zero ) Par = gnorm/dxnorm
!
!     beginning of an iteration.
!
 50      iter = iter + 1
!
!        evaluate the function at the current value of par.
!
         IF ( Par==zero ) Par = DMAX1(dwarf,p001*paru)
         temp = DSQRT(Par)
         DO j = 1 , N
            Wa1(j) = temp*Diag(j)
         ENDDO
         CALL QRSOLV(N,R,Ldr,Ipvt,Wa1,Qtb,X,Sdiag,Wa2)
         DO j = 1 , N
            Wa2(j) = Diag(j)*X(j)
         ENDDO
         dxnorm = ENORM(N,Wa2)
         temp = fp
         fp = dxnorm - Delta
!
!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.
!
         IF ( DABS(fp)<=p1*Delta .OR. parl==zero .AND. fp<=temp .AND.   &
            & temp<zero .OR. iter==10 ) THEN
            IF ( iter==0 ) Par = zero
         ELSE
!
!        compute the newton correction.
!
            DO j = 1 , N
               l = Ipvt(j)
               Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
            ENDDO
            DO j = 1 , N
               Wa1(j) = Wa1(j)/Sdiag(j)
               temp = Wa1(j)
               jp1 = j + 1
               IF ( N>=jp1 ) THEN
                  DO i = jp1 , N
                     Wa1(i) = Wa1(i) - R(i,j)*temp
                  ENDDO
               ENDIF
            ENDDO
            temp = ENORM(N,Wa1)
            parc = ((fp/Delta)/temp)/temp
!
!        depending on the sign of the function, update parl or paru.
!
            IF ( fp>zero ) parl = DMAX1(parl,Par)
            IF ( fp<zero ) paru = DMIN1(paru,Par)
!
!        compute an improved estimate for par.
!
            Par = DMAX1(parl,Par+parc)
!
!        end of an iteration.
!
            GOTO 50
         ENDIF
      ENDIF
!
!     last card of subroutine lmpar.
!
      END
!*==LMSTR.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMSTR(FCN,M,N,X,Fvec,Fjac,Ldfjac,Ftol,Xtol,Gtol,Maxfev,&
                     & Diag,Mode,Factor,Nprint,Info,Nfev,Njev,Ipvt,Qtf, &
                     & Wa1,Wa2,Wa3,Wa4)
      IMPLICIT NONE
!*--LMSTR3506
      INTEGER M , N , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev ,   &
            & Njev
      INTEGER Ipvt(N)
      LOGICAL sing
      DOUBLE PRECISION Ftol , Xtol , Gtol , Factor
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N) , Diag(N) ,      &
                     & Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) , Wa4(M)
!     **********
!
!     subroutine lmstr
!
!     the purpose of lmstr is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm which uses minimal storage.
!     the user must provide a subroutine which calculates the
!     functions and the rows of the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                        maxfev,diag,mode,factor,nprint,info,nfev,
!                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the rows of the jacobian.
!         fcn must be declared in an external statement in the
!         user calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjrow,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m),fjrow(n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec.
!         if iflag = i calculate the (i-1)-st row of the
!         jacobian at x and return this vector in fjrow.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmstr.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array. the upper triangle of fjac
!         contains an upper triangular matrix r such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower triangular
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,lmpar,qrfac,rwupdt
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
!     jorge j. more
!
!     **********
      INTEGER i , iflag , iter , j , l
      DOUBLE PRECISION actred , delta , dirder , epsmch , fnorm ,       &
                     & fnorm1 , gnorm , one , par , pnorm , prered ,    &
                     & p1 , p5 , p25 , p75 , p0001 , ratio , sum ,      &
                     & temp , temp1 , temp2 , xnorm , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p1 , p5 , p25 , p75 , p0001 , zero/1.0D0 , 1.0D-1 ,    &
         & 5.0D-1 , 2.5D-1 , 7.5D-1 , 1.0D-4 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
      Njev = 0
!
!     check the input parameters for errors.
!
      IF ( N<=0 .OR. M<N .OR. Ldfjac<N .OR. Ftol<zero .OR.              &
         & Xtol<zero .OR. Gtol<zero .OR. Maxfev<=0 .OR. Factor<=zero )  &
         & GOTO 200
      IF ( Mode==2 ) THEN
         DO j = 1 , N
            IF ( Diag(j)<=zero ) GOTO 200
         ENDDO
      ENDIF
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      CALL FCN(M,N,X,Fvec,Wa3,iflag)
      Nfev = 1
      IF ( iflag<0 ) GOTO 200
      fnorm = ENORM(M,Fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      par = zero
      iter = 1
!
!     beginning of the outer loop.
!
!
!        if requested, call fcn to enable printing of iterates.
!
 100  IF ( Nprint>0 ) THEN
         iflag = 0
         IF ( MOD(iter-1,Nprint)==0 ) CALL FCN(M,N,X,Fvec,Wa3,iflag)
         IF ( iflag<0 ) GOTO 200
      ENDIF
!
!        compute the qr factorization of the jacobian matrix
!        calculated one row at a time, while simultaneously
!        forming (q transpose)*fvec and storing the first
!        n components in qtf.
!
      DO j = 1 , N
         Qtf(j) = zero
         DO i = 1 , N
            Fjac(i,j) = zero
         ENDDO
      ENDDO
      iflag = 2
      DO i = 1 , M
         CALL FCN(M,N,X,Fvec,Wa3,iflag)
         IF ( iflag<0 ) GOTO 200
         temp = Fvec(i)
         CALL RWUPDT(N,Fjac,Ldfjac,Wa3,Qtf,temp,Wa1,Wa2)
         iflag = iflag + 1
      ENDDO
      Njev = Njev + 1
!
!        if the jacobian is rank deficient, call qrfac to
!        reorder its columns and update the components of qtf.
!
      sing = .FALSE.
      DO j = 1 , N
         IF ( Fjac(j,j)==zero ) sing = .TRUE.
         Ipvt(j) = j
         Wa2(j) = ENORM(j,Fjac(1,j))
      ENDDO
      IF ( sing ) THEN
         CALL QRFAC(N,N,Fjac,Ldfjac,.TRUE.,Ipvt,N,Wa1,Wa2,Wa3)
         DO j = 1 , N
            IF ( Fjac(j,j)/=zero ) THEN
               sum = zero
               DO i = j , N
                  sum = sum + Fjac(i,j)*Qtf(i)
               ENDDO
               temp = -sum/Fjac(j,j)
               DO i = j , N
                  Qtf(i) = Qtf(i) + Fjac(i,j)*temp
               ENDDO
            ENDIF
            Fjac(j,j) = Wa1(j)
         ENDDO
      ENDIF
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      IF ( iter==1 ) THEN
         IF ( Mode/=2 ) THEN
            DO j = 1 , N
               Diag(j) = Wa2(j)
               IF ( Wa2(j)==zero ) Diag(j) = one
            ENDDO
         ENDIF
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         DO j = 1 , N
            Wa3(j) = Diag(j)*X(j)
         ENDDO
         xnorm = ENORM(N,Wa3)
         delta = Factor*xnorm
         IF ( delta==zero ) delta = Factor
      ENDIF
!
!        compute the norm of the scaled gradient.
!
      gnorm = zero
      IF ( fnorm/=zero ) THEN
         DO j = 1 , N
            l = Ipvt(j)
            IF ( Wa2(l)/=zero ) THEN
               sum = zero
               DO i = 1 , j
                  sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
               ENDDO
               gnorm = DMAX1(gnorm,DABS(sum/Wa2(l)))
            ENDIF
         ENDDO
      ENDIF
!
!        test for convergence of the gradient norm.
!
      IF ( gnorm<=Gtol ) Info = 4
      IF ( Info==0 ) THEN
!
!        rescale if necessary.
!
         IF ( Mode/=2 ) THEN
            DO j = 1 , N
               Diag(j) = DMAX1(Diag(j),Wa2(j))
            ENDDO
         ENDIF
!
!        beginning of the inner loop.
!
!
!           determine the levenberg-marquardt parameter.
!
 150     CALL LMPAR(N,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1,Wa2,Wa3,  &
                  & Wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
         DO j = 1 , N
            Wa1(j) = -Wa1(j)
            Wa2(j) = X(j) + Wa1(j)
            Wa3(j) = Diag(j)*Wa1(j)
         ENDDO
         pnorm = ENORM(N,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
         IF ( iter==1 ) delta = DMIN1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
         iflag = 1
         CALL FCN(M,N,Wa2,Wa4,Wa3,iflag)
         Nfev = Nfev + 1
         IF ( iflag>=0 ) THEN
            fnorm1 = ENORM(M,Wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            IF ( p1*fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            DO j = 1 , N
               Wa3(j) = zero
               l = Ipvt(j)
               temp = Wa1(l)
               DO i = 1 , j
                  Wa3(i) = Wa3(i) + Fjac(i,j)*temp
               ENDDO
            ENDDO
            temp1 = ENORM(N,Wa3)/fnorm
            temp2 = (DSQRT(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2+temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            IF ( prered/=zero ) ratio = actred/prered
!
!           update the step bound.
!
            IF ( ratio<=p25 ) THEN
               IF ( actred>=zero ) temp = p5
               IF ( actred<zero ) temp = p5*dirder/(dirder+p5*actred)
               IF ( p1*fnorm1>=fnorm .OR. temp<p1 ) temp = p1
               delta = temp*DMIN1(delta,pnorm/p1)
               par = par/temp
            ELSEIF ( par==zero .OR. ratio>=p75 ) THEN
               delta = pnorm/p5
               par = p5*par
            ENDIF
!
!           test for successful iteration.
!
            IF ( ratio>=p0001 ) THEN
!
!           successful iteration. update x, fvec, and their norms.
!
               DO j = 1 , N
                  X(j) = Wa2(j)
                  Wa2(j) = Diag(j)*X(j)
               ENDDO
               DO i = 1 , M
                  Fvec(i) = Wa4(i)
               ENDDO
               xnorm = ENORM(N,Wa2)
               fnorm = fnorm1
               iter = iter + 1
            ENDIF
!
!           tests for convergence.
!
            IF ( DABS(actred)<=Ftol .AND. prered<=Ftol .AND.            &
               & p5*ratio<=one ) Info = 1
            IF ( delta<=Xtol*xnorm ) Info = 2
            IF ( DABS(actred)<=Ftol .AND. prered<=Ftol .AND.            &
               & p5*ratio<=one .AND. Info==2 ) Info = 3
            IF ( Info==0 ) THEN
!
!           tests for termination and stringent tolerances.
!
               IF ( Nfev>=Maxfev ) Info = 5
               IF ( DABS(actred)<=epsmch .AND. prered<=epsmch .AND.     &
                  & p5*ratio<=one ) Info = 6
               IF ( delta<=epsmch*xnorm ) Info = 7
               IF ( gnorm<=epsmch ) Info = 8
               IF ( Info==0 ) THEN
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
!
!        end of the outer loop.
!
                  IF ( ratio>=p0001 ) GOTO 100
                  GOTO 150
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
!     termination, either normal or user imposed.
!
 200  IF ( iflag<0 ) Info = iflag
      iflag = 0
      IF ( Nprint>0 ) CALL FCN(M,N,X,Fvec,Wa3,iflag)
!
!     last card of subroutine lmstr.
!
      END
!*==LMSTR1.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE LMSTR1(FCN,M,N,X,Fvec,Fjac,Ldfjac,Tol,Info,Ipvt,Wa,Lwa)
      IMPLICIT NONE
!*--LMSTR13971
!*** Start of declarations inserted by SPAG
      REAL FCN
!*** End of declarations inserted by SPAG
      INTEGER M , N , Ldfjac , Info , Lwa
      INTEGER Ipvt(N)
      DOUBLE PRECISION Tol
      DOUBLE PRECISION X(N) , Fvec(M) , Fjac(Ldfjac,N) , Wa(Lwa)
      EXTERNAL FCN
!     **********
!
!     subroutine lmstr1
!
!     the purpose of lmstr1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm which uses minimal storage.
!     this is done by using the more general least-squares solver
!     lmstr. the user must provide a subroutine which calculates
!     the functions and the rows of the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
!                         ipvt,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the rows of the jacobian.
!         fcn must be declared in an external statement in the
!         user calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjrow,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m),fjrow(n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec.
!         if iflag = i calculate the (i-1)-st row of the
!         jacobian at x and return this vector in fjrow.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmstr1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array. the upper triangle of fjac
!         contains an upper triangular matrix r such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower triangular
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached 100*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than 5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmstr
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
!     jorge j. more
!
!     **********
      INTEGER maxfev , mode , nfev , njev , nprint
      DOUBLE PRECISION factor , ftol , gtol , xtol , zero
      DATA factor , zero/1.0D2 , 0.0D0/
      Info = 0
!
!     check the input parameters for errors.
!
      IF ( N>0 .AND. M>=N .AND. Ldfjac>=N .AND. Tol>=zero .AND.         &
         & Lwa>=5*N+M ) THEN
!
!     call lmstr.
!
         maxfev = 100*(N+1)
         ftol = Tol
         xtol = Tol
         gtol = zero
         mode = 1
         nprint = 0
         CALL LMSTR(FCN,M,N,X,Fvec,Fjac,Ldfjac,ftol,xtol,gtol,maxfev,   &
                  & Wa(1),mode,factor,nprint,Info,nfev,njev,Ipvt,Wa(N+1)&
                  & ,Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
         IF ( Info==8 ) Info = 4
      ENDIF
!
!     last card of subroutine lmstr1.
!
      END
!*==QFORM.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE QFORM(M,N,Q,Ldq,Wa)
      IMPLICIT NONE
!*--QFORM4131
      INTEGER M , N , Ldq
      DOUBLE PRECISION Q(Ldq,M) , Wa(M)
!     **********
!
!     subroutine qform
!
!     this subroutine proceeds from the computed qr factorization of
!     an m by n matrix a to accumulate the m by m orthogonal matrix
!     q from its factored form.
!
!     the subroutine statement is
!
!       subroutine qform(m,n,q,ldq,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a and the order of q.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       q is an m by m array. on input the full lower trapezoid in
!         the first min(m,n) columns of q contains the factored form.
!         on output q has been accumulated into a square matrix.
!
!       ldq is a positive integer input variable not less than m
!         which specifies the leading dimension of the array q.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       fortran-supplied ... min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jm1 , k , l , minmn , np1
      DOUBLE PRECISION one , sum , temp , zero
      DATA one , zero/1.0D0 , 0.0D0/
!
!     zero out upper triangle of q in the first min(m,n) columns.
!
      minmn = MIN0(M,N)
      IF ( minmn>=2 ) THEN
         DO j = 2 , minmn
            jm1 = j - 1
            DO i = 1 , jm1
               Q(i,j) = zero
            ENDDO
         ENDDO
      ENDIF
!
!     initialize remaining columns to those of the identity matrix.
!
      np1 = N + 1
      IF ( M>=np1 ) THEN
         DO j = np1 , M
            DO i = 1 , M
               Q(i,j) = zero
            ENDDO
            Q(j,j) = one
         ENDDO
      ENDIF
!
!     accumulate q from its factored form.
!
      DO l = 1 , minmn
         k = minmn - l + 1
         DO i = k , M
            Wa(i) = Q(i,k)
            Q(i,k) = zero
         ENDDO
         Q(k,k) = one
         IF ( Wa(k)/=zero ) THEN
            DO j = k , M
               sum = zero
               DO i = k , M
                  sum = sum + Q(i,j)*Wa(i)
               ENDDO
               temp = sum/Wa(k)
               DO i = k , M
                  Q(i,j) = Q(i,j) - temp*Wa(i)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
!     last card of subroutine qform.
!
      END
!*==QRFAC.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE QRFAC(M,N,A,Lda,Pivot,Ipvt,Lipvt,Rdiag,Acnorm,Wa)
      IMPLICIT NONE
!*--QRFAC4228
      INTEGER M , N , Lda , Lipvt
      INTEGER Ipvt(Lipvt)
      LOGICAL Pivot
      DOUBLE PRECISION A(Lda,N) , Rdiag(N) , Acnorm(N) , Wa(N)
!     **********
!
!     subroutine qrfac
!
!     this subroutine uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!     the subroutine statement is
!
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed. on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. if pivot is set true,
!         then column pivoting is enforced. if pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dmax1,dsqrt,min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jp1 , k , kmax , minmn
      DOUBLE PRECISION ajnorm , epsmch , one , p05 , sum , temp , zero
      DOUBLE PRECISION DPMPAR , ENORM
      DATA one , p05 , zero/1.0D0 , 5.0D-2 , 0.0D0/
!
!     epsmch is the machine precision.
!
      epsmch = DPMPAR(1)
!
!     compute the initial column norms and initialize several arrays.
!
      DO j = 1 , N
         Acnorm(j) = ENORM(M,A(1,j))
         Rdiag(j) = Acnorm(j)
         Wa(j) = Rdiag(j)
         IF ( Pivot ) Ipvt(j) = j
      ENDDO
!
!     reduce a to r with householder transformations.
!
      minmn = MIN0(M,N)
      DO j = 1 , minmn
         IF ( Pivot ) THEN
!
!        bring the column of largest norm into the pivot position.
!
            kmax = j
            DO k = j , N
               IF ( Rdiag(k)>Rdiag(kmax) ) kmax = k
            ENDDO
            IF ( kmax/=j ) THEN
               DO i = 1 , M
                  temp = A(i,j)
                  A(i,j) = A(i,kmax)
                  A(i,kmax) = temp
               ENDDO
               Rdiag(kmax) = Rdiag(j)
               Wa(kmax) = Wa(j)
               k = Ipvt(j)
               Ipvt(j) = Ipvt(kmax)
               Ipvt(kmax) = k
            ENDIF
         ENDIF
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = ENORM(M-j+1,A(j,j))
         IF ( ajnorm/=zero ) THEN
            IF ( A(j,j)<zero ) ajnorm = -ajnorm
            DO i = j , M
               A(i,j) = A(i,j)/ajnorm
            ENDDO
            A(j,j) = A(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
            jp1 = j + 1
            IF ( N>=jp1 ) THEN
               DO k = jp1 , N
                  sum = zero
                  DO i = j , M
                     sum = sum + A(i,j)*A(i,k)
                  ENDDO
                  temp = sum/A(j,j)
                  DO i = j , M
                     A(i,k) = A(i,k) - temp*A(i,j)
                  ENDDO
                  IF ( .NOT.(.NOT.Pivot .OR. Rdiag(k)==zero) ) THEN
                     temp = A(j,k)/Rdiag(k)
                     Rdiag(k) = Rdiag(k)*DSQRT(DMAX1(zero,one-temp**2))
                     IF ( p05*(Rdiag(k)/Wa(k))**2<=epsmch ) THEN
                        Rdiag(k) = ENORM(M-j,A(jp1,k))
                        Wa(k) = Rdiag(k)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         Rdiag(j) = -ajnorm
      ENDDO
!
!     last card of subroutine qrfac.
!
      END
!*==QRSOLV.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE QRSOLV(N,R,Ldr,Ipvt,Diag,Qtb,X,Sdiag,Wa)
      IMPLICIT NONE
!*--QRSOLV4397
      INTEGER N , Ldr
      INTEGER Ipvt(N)
      DOUBLE PRECISION R(Ldr,N) , Diag(N) , Qtb(N) , X(N) , Sdiag(N) ,  &
                     & Wa(N)
!     **********
!
!     subroutine qrsolv
!
!     given an m by n matrix a, an n by n diagonal matrix d,
!     and an m-vector b, the problem is to determine an x which
!     solves the system
!
!           a*x = b ,     d*x = 0 ,
!
!     in the least squares sense.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then qrsolv expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. the system
!     a*x = b, d*x = 0, is then equivalent to
!
!                  t       t
!           r*z = q *b ,  p *d*p*z = 0 ,
!
!     where x = p*z. if this system does not have full rank,
!     then a least squares solution is obtained. on output qrsolv
!     also provides an upper triangular matrix s such that
!
!            t   t               t
!           p *(a *a + d*d)*p = s *s .
!
!     s is computed within qrsolv and may be of separate interest.
!
!     the subroutine statement is
!
!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, d*x = 0.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa is a work array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , jp1 , k , kp1 , l , nsing
      DOUBLE PRECISION cos , cotan , p5 , p25 , qtbpj , sin , sum ,     &
                     & tan , temp , zero
      DATA p5 , p25 , zero/5.0D-1 , 2.5D-1 , 0.0D0/
!
!     copy r and (q transpose)*b to preserve input and initialize s.
!     in particular, save the diagonal elements of r in x.
!
      DO j = 1 , N
         DO i = j , N
            R(i,j) = R(j,i)
         ENDDO
         X(j) = R(j,j)
         Wa(j) = Qtb(j)
      ENDDO
!
!     eliminate the diagonal matrix d using a givens rotation.
!
      DO j = 1 , N
!
!        prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
!
         l = Ipvt(j)
         IF ( Diag(l)/=zero ) THEN
            DO k = j , N
               Sdiag(k) = zero
            ENDDO
            Sdiag(j) = Diag(l)
!
!        the transformations to eliminate the row of d
!        modify only a single element of (q transpose)*b
!        beyond the first n, which is initially zero.
!
            qtbpj = zero
            DO k = j , N
!
!           determine a givens rotation which eliminates the
!           appropriate element in the current row of d.
!
               IF ( Sdiag(k)/=zero ) THEN
                  IF ( DABS(R(k,k))>=DABS(Sdiag(k)) ) THEN
                     tan = Sdiag(k)/R(k,k)
                     cos = p5/DSQRT(p25+p25*tan**2)
                     sin = cos*tan
                  ELSE
                     cotan = R(k,k)/Sdiag(k)
                     sin = p5/DSQRT(p25+p25*cotan**2)
                     cos = sin*cotan
                  ENDIF
!
!           compute the modified diagonal element of r and
!           the modified element of ((q transpose)*b,0).
!
                  R(k,k) = cos*R(k,k) + sin*Sdiag(k)
                  temp = cos*Wa(k) + sin*qtbpj
                  qtbpj = -sin*Wa(k) + cos*qtbpj
                  Wa(k) = temp
!
!           accumulate the tranformation in the row of s.
!
                  kp1 = k + 1
                  IF ( N>=kp1 ) THEN
                     DO i = kp1 , N
                        temp = cos*R(i,k) + sin*Sdiag(i)
                        Sdiag(i) = -sin*R(i,k) + cos*Sdiag(i)
                        R(i,k) = temp
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
!
!        store the diagonal element of s and restore
!        the corresponding diagonal element of r.
!
         Sdiag(j) = R(j,j)
         R(j,j) = X(j)
      ENDDO
!
!     solve the triangular system for z. if the system is
!     singular, then obtain a least squares solution.
!
      nsing = N
      DO j = 1 , N
         IF ( Sdiag(j)==zero .AND. nsing==N ) nsing = j - 1
         IF ( nsing<N ) Wa(j) = zero
      ENDDO
      IF ( nsing>=1 ) THEN
         DO k = 1 , nsing
            j = nsing - k + 1
            sum = zero
            jp1 = j + 1
            IF ( nsing>=jp1 ) THEN
               DO i = jp1 , nsing
                  sum = sum + R(i,j)*Wa(i)
               ENDDO
            ENDIF
            Wa(j) = (Wa(j)-sum)/Sdiag(j)
         ENDDO
      ENDIF
!
!     permute the components of z back to components of x.
!
      DO j = 1 , N
         l = Ipvt(j)
         X(l) = Wa(j)
      ENDDO
!
!     last card of subroutine qrsolv.
!
      END
!*==R1MPYQ.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE R1MPYQ(M,N,A,Lda,V,W)
      IMPLICIT NONE
!*--R1MPYQ4594
      INTEGER M , N , Lda
      DOUBLE PRECISION A(Lda,N) , V(N) , W(N)
!     **********
!
!     subroutine r1mpyq
!
!     given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     and gv(i), gw(i) are givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.
!
!     the subroutine statement is
!
!       subroutine r1mpyq(m,n,a,lda,v,w)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a must contain the matrix
!         to be postmultiplied by the orthogonal matrix q
!         described above. on output a*q has replaced a.
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       v is an input array of length n. v(i) must contain the
!         information necessary to recover the givens rotation gv(i)
!         described above.
!
!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the givens rotation gw(i)
!         described above.
!
!     subroutines called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i , j , nmj , nm1
      DOUBLE PRECISION cos , one , sin , temp
      DATA one/1.0D0/
!
!     apply the first set of givens rotations to a.
!
      nm1 = N - 1
      IF ( nm1>=1 ) THEN
         DO nmj = 1 , nm1
            j = N - nmj
            IF ( DABS(V(j))>one ) cos = one/V(j)
            IF ( DABS(V(j))>one ) sin = DSQRT(one-cos**2)
            IF ( DABS(V(j))<=one ) sin = V(j)
            IF ( DABS(V(j))<=one ) cos = DSQRT(one-sin**2)
            DO i = 1 , M
               temp = cos*A(i,j) - sin*A(i,N)
               A(i,N) = sin*A(i,j) + cos*A(i,N)
               A(i,j) = temp
            ENDDO
         ENDDO
!
!     apply the second set of givens rotations to a.
!
         DO j = 1 , nm1
            IF ( DABS(W(j))>one ) cos = one/W(j)
            IF ( DABS(W(j))>one ) sin = DSQRT(one-cos**2)
            IF ( DABS(W(j))<=one ) sin = W(j)
            IF ( DABS(W(j))<=one ) cos = DSQRT(one-sin**2)
            DO i = 1 , M
               temp = cos*A(i,j) + sin*A(i,N)
               A(i,N) = -sin*A(i,j) + cos*A(i,N)
               A(i,j) = temp
            ENDDO
         ENDDO
      ENDIF
!
!     last card of subroutine r1mpyq.
!
      END
!*==R1UPDT.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE R1UPDT(M,N,S,Ls,U,V,W,Sing)
      IMPLICIT NONE
!*--R1UPDT4688
      INTEGER M , N , Ls
      LOGICAL Sing
      DOUBLE PRECISION S(Ls) , U(M) , V(N) , W(M)
!     **********
!
!     subroutine r1updt
!
!     given an m by n lower trapezoidal matrix s, an m-vector u,
!     and an n-vector v, the problem is to determine an
!     orthogonal matrix q such that
!
!                   t
!           (s + u*v )*q
!
!     is again lower trapezoidal.
!
!     this subroutine determines q as the product of 2*(n - 1)
!     transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     where gv(i), gw(i) are givens rotations in the (i,n) plane
!     which eliminate elements in the i-th and n-th planes,
!     respectively. q itself is not accumulated, rather the
!     information to recover the gv, gw rotations is returned.
!
!     the subroutine statement is
!
!       subroutine r1updt(m,n,s,ls,u,v,w,sing)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of s.
!
!       n is a positive integer input variable set to the number
!         of columns of s. n must not exceed m.
!
!       s is an array of length ls. on input s must contain the lower
!         trapezoidal matrix s stored by columns. on output s contains
!         the lower trapezoidal matrix produced as described above.
!
!       ls is a positive integer input variable not less than
!         (n*(2*m-n+1))/2.
!
!       u is an input array of length m which must contain the
!         vector u.
!
!       v is an array of length n. on input v must contain the vector
!         v. on output v(i) contains the information necessary to
!         recover the givens rotation gv(i) described above.
!
!       w is an output array of length m. w(i) contains information
!         necessary to recover the givens rotation gw(i) described
!         above.
!
!       sing is a logical output variable. sing is set true if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set false.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more,
!     john l. nazareth
!
!     **********
      INTEGER i , j , jj , l , nmj , nm1
      DOUBLE PRECISION cos , cotan , giant , one , p5 , p25 , sin ,     &
                     & tan , tau , temp , zero
      DOUBLE PRECISION DPMPAR
      DATA one , p5 , p25 , zero/1.0D0 , 5.0D-1 , 2.5D-1 , 0.0D0/
!
!     giant is the largest magnitude.
!
      giant = DPMPAR(3)
!
!     initialize the diagonal element pointer.
!
      jj = (N*(2*M-N+1))/2 - (M-N)
!
!     move the nontrivial part of the last column of s into w.
!
      l = jj
      DO i = N , M
         W(i) = S(l)
         l = l + 1
      ENDDO
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
      nm1 = N - 1
      IF ( nm1>=1 ) THEN
         DO nmj = 1 , nm1
            j = N - nmj
            jj = jj - (M-j+1)
            W(j) = zero
            IF ( V(j)/=zero ) THEN
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
               IF ( DABS(V(N))>=DABS(V(j)) ) THEN
                  tan = V(j)/V(N)
                  cos = p5/DSQRT(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               ELSE
                  cotan = V(N)/V(j)
                  sin = p5/DSQRT(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  IF ( DABS(cos)*giant>one ) tau = one/cos
               ENDIF
!
!        apply the transformation to v and store the information
!        necessary to recover the givens rotation.
!
               V(N) = sin*V(j) + cos*V(N)
               V(j) = tau
!
!        apply the transformation to s and extend the spike in w.
!
               l = jj
               DO i = j , M
                  temp = cos*S(l) - sin*W(i)
                  W(i) = sin*S(l) + cos*W(i)
                  S(l) = temp
                  l = l + 1
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!     add the spike from the rank 1 update to w.
!
      DO i = 1 , M
         W(i) = W(i) + V(N)*U(i)
      ENDDO
!
!     eliminate the spike.
!
      Sing = .FALSE.
      IF ( nm1>=1 ) THEN
         DO j = 1 , nm1
            IF ( W(j)/=zero ) THEN
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
               IF ( DABS(S(jj))>=DABS(W(j)) ) THEN
                  tan = W(j)/S(jj)
                  cos = p5/DSQRT(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               ELSE
                  cotan = S(jj)/W(j)
                  sin = p5/DSQRT(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  IF ( DABS(cos)*giant>one ) tau = one/cos
               ENDIF
!
!        apply the transformation to s and reduce the spike in w.
!
               l = jj
               DO i = j , M
                  temp = cos*S(l) + sin*W(i)
                  W(i) = -sin*S(l) + cos*W(i)
                  S(l) = temp
                  l = l + 1
               ENDDO
!
!        store the information necessary to recover the
!        givens rotation.
!
               W(j) = tau
            ENDIF
!
!        test for zero diagonal elements in the output s.
!
            IF ( S(jj)==zero ) Sing = .TRUE.
            jj = jj + (M-j+1)
         ENDDO
      ENDIF
!
!     move w back into the last column of the output s.
!
      l = jj
      DO i = N , M
         S(l) = W(i)
         l = l + 1
      ENDDO
      IF ( S(jj)==zero ) Sing = .TRUE.
!
!     last card of subroutine r1updt.
!
      END
!*==RWUPDT.spg  processed by SPAG 6.72Dc at 03:59 on 19 Sep 2021
      SUBROUTINE RWUPDT(N,R,Ldr,W,B,Alpha,Cos,Sin)
      IMPLICIT NONE
!*--RWUPDT4895
      INTEGER N , Ldr
      DOUBLE PRECISION Alpha
      DOUBLE PRECISION R(Ldr,N) , W(N) , B(N) , Cos(N) , Sin(N)
!     **********
!
!     subroutine rwupdt
!
!     given an n by n upper triangular matrix r, this subroutine
!     computes the qr decomposition of the matrix formed when a row
!     is added to r. if the row is specified by the vector w, then
!     rwupdt determines an orthogonal matrix q such that when the
!     n+1 by n matrix composed of r augmented by w is premultiplied
!     by (q transpose), the resulting matrix is upper trapezoidal.
!     the matrix (q transpose) is the product of n transformations
!
!           g(n)*g(n-1)* ... *g(1)
!
!     where g(i) is a givens rotation in the (i,n+1) plane which
!     eliminates elements in the (n+1)-st plane. rwupdt also
!     computes the product (q transpose)*c where c is the
!     (n+1)-vector (b,alpha). q itself is not accumulated, rather
!     the information to recover the g rotations is supplied.
!
!     the subroutine statement is
!
!       subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the upper triangular part of
!         r must contain the matrix to be updated. on output r
!         contains the updated triangular matrix.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       w is an input array of length n which must contain the row
!         vector to be added to r.
!
!       b is an array of length n. on input b must contain the
!         first n elements of the vector c. on output b contains
!         the first n elements of the vector (q transpose)*c.
!
!       alpha is a variable. on input alpha must contain the
!         (n+1)-st element of the vector c. on output alpha contains
!         the (n+1)-st element of the vector (q transpose)*c.
!
!       cos is an output array of length n which contains the
!         cosines of the transforming givens rotations.
!
!       sin is an output array of length n which contains the
!         sines of the transforming givens rotations.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
!     jorge j. more
!
!     **********
      INTEGER i , j , jm1
      DOUBLE PRECISION cotan , one , p5 , p25 , rowj , tan , temp , zero
      DATA one , p5 , p25 , zero/1.0D0 , 5.0D-1 , 2.5D-1 , 0.0D0/
!
      DO j = 1 , N
         rowj = W(j)
         jm1 = j - 1
!
!        apply the previous transformations to
!        r(i,j), i=1,2,...,j-1, and to w(j).
!
         IF ( jm1>=1 ) THEN
            DO i = 1 , jm1
               temp = Cos(i)*R(i,j) + Sin(i)*rowj
               rowj = -Sin(i)*R(i,j) + Cos(i)*rowj
               R(i,j) = temp
            ENDDO
         ENDIF
!
!        determine a givens rotation which eliminates w(j).
!
         Cos(j) = one
         Sin(j) = zero
         IF ( rowj/=zero ) THEN
            IF ( DABS(R(j,j))>=DABS(rowj) ) THEN
               tan = rowj/R(j,j)
               Cos(j) = p5/DSQRT(p25+p25*tan**2)
               Sin(j) = Cos(j)*tan
            ELSE
               cotan = R(j,j)/rowj
               Sin(j) = p5/DSQRT(p25+p25*cotan**2)
               Cos(j) = Sin(j)*cotan
            ENDIF
!
!        apply the current transformation to r(j,j), b(j), and alpha.
!
            R(j,j) = Cos(j)*R(j,j) + Sin(j)*rowj
            temp = Cos(j)*B(j) + Sin(j)*Alpha
            Alpha = -Sin(j)*B(j) + Cos(j)*Alpha
            B(j) = temp
         ENDIF
      ENDDO
!
!     last card of subroutine rwupdt.
!
      END
