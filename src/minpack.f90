
      subroutine chkder(m,n,x,Fvec,Fjac,Ldfjac,Xp,Fvecp,Mode,Err)
      implicit none

      integer m , n , Ldfjac , Mode
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n) , Xp(n) ,        &
                     & Fvecp(m) , Err(m)
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
      integer i , j
      double precision eps , epsf , epslog , epsmch , factor , one ,    &
                     & temp , zero
      double precision dpmpar
      data factor , one , zero/1.0d2 , 1.0d0 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(epsmch)
!
      if ( Mode==2 ) then
!
!        mode = 2.
!
         epsf = factor*epsmch
         epslog = dlog10(eps)
         do i = 1 , m
            Err(i) = zero
         enddo
         do j = 1 , n
            temp = dabs(x(j))
            if ( temp==zero ) temp = one
            do i = 1 , m
               Err(i) = Err(i) + temp*Fjac(i,j)
            enddo
         enddo
         do i = 1 , m
            temp = one
            if ( Fvec(i)/=zero .and. Fvecp(i)/=zero .and.               &
               & dabs(Fvecp(i)-Fvec(i))>=epsf*dabs(Fvec(i)) )           &
               & temp = eps*dabs((Fvecp(i)-Fvec(i))/eps-Err(i))         &
               & /(dabs(Fvec(i))+dabs(Fvecp(i)))
            Err(i) = one
            if ( temp>epsmch .and. temp<eps ) Err(i)                    &
               & = (dlog10(temp)-epslog)/epslog
            if ( temp>=eps ) Err(i) = zero
         enddo
      else
!
!        mode = 1.
!
         do j = 1 , n
            temp = eps*dabs(x(j))
            if ( temp==zero ) temp = eps
            Xp(j) = x(j) + temp
         enddo
      endif
!
!
!     last card of subroutine chkder.
!
      end

      subroutine dogleg(n,r,Lr,Diag,Qtb,Delta,x,Wa1,Wa2)
      implicit none

      integer n , Lr
      double precision Delta
      double precision r(Lr) , Diag(n) , Qtb(n) , x(n) , Wa1(n) , Wa2(n)
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
      integer i , j , jj , jp1 , k , l
      double precision alpha , bnorm , epsmch , gnorm , one , qnorm ,   &
                     & sgnorm , sum , temp , zero
      double precision dpmpar , enorm
      data one , zero/1.0d0 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
!     first, calculate the gauss-newton direction.
!
      jj = (n*(n+1))/2 + 1
      do k = 1 , n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if ( n>=jp1 ) then
            do i = jp1 , n
               sum = sum + r(l)*x(i)
               l = l + 1
            enddo
         endif
         temp = r(jj)
         if ( temp==zero ) then
            l = j
            do i = 1 , j
               temp = dmax1(temp,dabs(r(l)))
               l = l + n - i
            enddo
            temp = epsmch*temp
            if ( temp==zero ) temp = epsmch
         endif
         x(j) = (Qtb(j)-sum)/temp
      enddo
!
!     test whether the gauss-newton direction is acceptable.
!
      do j = 1 , n
         Wa1(j) = zero
         Wa2(j) = Diag(j)*x(j)
      enddo
      qnorm = enorm(n,Wa2)
      if ( qnorm>Delta ) then
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
         l = 1
         do j = 1 , n
            temp = Qtb(j)
            do i = j , n
               Wa1(i) = Wa1(i) + r(l)*temp
               l = l + 1
            enddo
            Wa1(j) = Wa1(j)/Diag(j)
         enddo
!
!     calculate the norm of the scaled gradient and test for
!     the special case in which the scaled gradient is zero.
!
         gnorm = enorm(n,Wa1)
         sgnorm = zero
         alpha = Delta/qnorm
         if ( gnorm/=zero ) then
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
            do j = 1 , n
               Wa1(j) = (Wa1(j)/gnorm)/Diag(j)
            enddo
            l = 1
            do j = 1 , n
               sum = zero
               do i = j , n
                  sum = sum + r(l)*Wa1(i)
                  l = l + 1
               enddo
               Wa2(j) = sum
            enddo
            temp = enorm(n,Wa2)
            sgnorm = (gnorm/temp)/temp
!
!     test whether the scaled gradient direction is acceptable.
!
            alpha = zero
            if ( sgnorm<Delta ) then
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
               bnorm = enorm(n,Qtb)
               temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/Delta)
               temp = temp - (Delta/qnorm)*(sgnorm/Delta)               &
                    & **2 + dsqrt((temp-(Delta/qnorm))                  &
                    & **2+(one-(Delta/qnorm)**2)*(one-(sgnorm/Delta)**2)&
                    & )
               alpha = ((Delta/qnorm)*(one-(sgnorm/Delta)**2))/temp
            endif
         endif
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
         temp = (one-alpha)*dmin1(sgnorm,Delta)
         do j = 1 , n
            x(j) = temp*Wa1(j) + alpha*x(j)
         enddo
      endif
!
!     last card of subroutine dogleg.
!
      end

      double precision function dpmpar(i)
      implicit none

      integer i
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
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      double precision dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
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
      data dmach(1)/2.22044604926d-16/
      data dmach(2)/2.22507385852d-308/
      data dmach(3)/1.79769313485d+308/
!
      dpmpar = dmach(i)
!
!     Last card of function dpmpar.
!
      end

      double precision function enorm(n,x)
      implicit none

      integer n
      double precision x(n)
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
      integer i
      double precision agiant , floatn , one , rdwarf , rgiant , s1 ,   &
                     & s2 , s3 , xabs , x1max , x3max , zero
      data one , zero , rdwarf , rgiant/1.0d0 , 0.0d0 , 3.834d-20 ,     &
         & 1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do i = 1 , n
         xabs = dabs(x(i))
         if ( xabs>rdwarf .and. xabs<agiant ) then
!
!           sum for intermediate components.
!
            s2 = s2 + xabs**2
         elseif ( xabs<=rdwarf ) then
!
!              sum for small components.
!
            if ( xabs<=x3max ) then
               if ( xabs/=zero ) s3 = s3 + (xabs/x3max)**2
            else
               s3 = one + s3*(x3max/xabs)**2
               x3max = xabs
            endif
!
!              sum for large components.
!
         elseif ( xabs<=x1max ) then
            s1 = s1 + (xabs/x1max)**2
         else
            s1 = one + s1*(x1max/xabs)**2
            x1max = xabs
         endif
      enddo
!
!     calculation of norm.
!
      if ( s1/=zero ) then
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
      elseif ( s2==zero ) then
         enorm = x3max*dsqrt(s3)
      else
         if ( s2>=x3max ) enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
         if ( s2<x3max ) enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
      endif
!
!     last card of function enorm.
!
      end

      subroutine fdjac1(fcn,n,x,Fvec,Fjac,Ldfjac,Iflag,Ml,Mu,Epsfcn,Wa1,&
                      & Wa2)
      implicit none

      integer n , Ldfjac , Iflag , Ml , Mu
      double precision Epsfcn
      double precision x(n) , Fvec(n) , Fjac(Ldfjac,n) , Wa1(n) , Wa2(n)
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
      integer i , j , k , msum
      double precision eps , epsmch , h , temp , zero
      double precision dpmpar
      data zero/0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(dmax1(Epsfcn,epsmch))
      msum = Ml + Mu + 1
      if ( msum<n ) then
!
!        computation of banded approximate jacobian.
!
         do k = 1 , msum
            do j = k , n , msum
               Wa2(j) = x(j)
               h = eps*dabs(Wa2(j))
               if ( h==zero ) h = eps
               x(j) = Wa2(j) + h
            enddo
            call fcn(n,x,Wa1,Iflag)
            if ( Iflag<0 ) goto 99999
            do j = k , n , msum
               x(j) = Wa2(j)
               h = eps*dabs(Wa2(j))
               if ( h==zero ) h = eps
               do i = 1 , n
                  Fjac(i,j) = zero
                  if ( i>=j-Mu .and. i<=j+Ml ) Fjac(i,j)                &
                     & = (Wa1(i)-Fvec(i))/h
               enddo
            enddo
         enddo
      else
!
!        computation of dense approximate jacobian.
!
         do j = 1 , n
            temp = x(j)
            h = eps*dabs(temp)
            if ( h==zero ) h = eps
            x(j) = temp + h
            call fcn(n,x,Wa1,Iflag)
            if ( Iflag<0 ) goto 99999
            x(j) = temp
            do i = 1 , n
               Fjac(i,j) = (Wa1(i)-Fvec(i))/h
            enddo
         enddo
      endif
!
!     last card of subroutine fdjac1.
!
99999 end

 
      subroutine fdjac2(fcn,m,n,x,Fvec,Fjac,Ldfjac,Iflag,Epsfcn,Wa)
      implicit none

      integer m , n , Ldfjac , Iflag
      double precision Epsfcn
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n) , Wa(m)
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
      integer i , j
      double precision eps , epsmch , h , temp , zero
      double precision dpmpar
      data zero/0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(dmax1(Epsfcn,epsmch))
      do j = 1 , n
         temp = x(j)
         h = eps*dabs(temp)
         if ( h==zero ) h = eps
         x(j) = temp + h
         call fcn(m,n,x,Wa,Iflag)
         if ( Iflag<0 ) goto 99999
         x(j) = temp
         do i = 1 , m
            Fjac(i,j) = (Wa(i)-Fvec(i))/h
         enddo
      enddo
!
!     last card of subroutine fdjac2.
!
99999 end

      subroutine hybrd(fcn,n,x,Fvec,Xtol,Maxfev,Ml,Mu,Epsfcn,Diag,Mode, &
                     & Factor,Nprint,Info,Nfev,Fjac,Ldfjac,r,Lr,Qtf,Wa1,&
                     & Wa2,Wa3,Wa4)
      implicit none

      integer n , Maxfev , Ml , Mu , Mode , Nprint , Info , Nfev ,      &
            & Ldfjac , Lr
      double precision Xtol , Epsfcn , Factor
      double precision x(n) , Fvec(n) , Diag(n) , Fjac(Ldfjac,n) ,      &
                     & r(Lr) , Qtf(n) , Wa1(n) , Wa2(n) , Wa3(n) ,      &
                     & Wa4(n)
      external fcn
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
      integer i , iflag , iter , j , jm1 , l , msum , ncfail , ncsuc ,  &
            & nslow1 , nslow2
      integer iwa(1)
      logical jeval , sing
      double precision actred , delta , epsmch , fnorm , fnorm1 , one , &
                     & pnorm , prered , p1 , p5 , p001 , p0001 , ratio ,&
                     & sum , temp , xnorm , zero
      double precision dpmpar , enorm
      data one , p1 , p5 , p001 , p0001 , zero/1.0d0 , 1.0d-1 , 5.0d-1 ,&
         & 1.0d-3 , 1.0d-4 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
!
!     check the input parameters for errors.
!
      if ( n<=0 .or. Xtol<zero .or. Maxfev<=0 .or. Ml<0 .or. Mu<0 .or.  &
         & Factor<=zero .or. Ldfjac<n .or. Lr<(n*(n+1))/2 ) goto 300
      if ( Mode==2 ) then
         do j = 1 , n
            if ( Diag(j)<=zero ) goto 300
         enddo
      endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,Fvec,iflag)
      Nfev = 1
      if ( iflag<0 ) goto 300
      fnorm = enorm(n,Fvec)
!
!     determine the number of calls to fcn needed to compute
!     the jacobian matrix.
!
      msum = min0(Ml+Mu+1,n)
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
 100  jeval = .true.
!
!        calculate the jacobian matrix.
!
      iflag = 2
      call fdjac1(fcn,n,x,Fvec,Fjac,Ldfjac,iflag,Ml,Mu,Epsfcn,Wa1,Wa2)
      Nfev = Nfev + msum
      if ( iflag<0 ) goto 300
!
!        compute the qr factorization of the jacobian.
!
      call qrfac(n,n,Fjac,Ldfjac,.false.,iwa,1,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      if ( iter==1 ) then
         if ( Mode/=2 ) then
            do j = 1 , n
               Diag(j) = Wa2(j)
               if ( Wa2(j)==zero ) Diag(j) = one
            enddo
         endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do j = 1 , n
            Wa3(j) = Diag(j)*x(j)
         enddo
         xnorm = enorm(n,Wa3)
         delta = Factor*xnorm
         if ( delta==zero ) delta = Factor
      endif
!
!        form (q transpose)*fvec and store in qtf.
!
      do i = 1 , n
         Qtf(i) = Fvec(i)
      enddo
      do j = 1 , n
         if ( Fjac(j,j)/=zero ) then
            sum = zero
            do i = j , n
               sum = sum + Fjac(i,j)*Qtf(i)
            enddo
            temp = -sum/Fjac(j,j)
            do i = j , n
               Qtf(i) = Qtf(i) + Fjac(i,j)*temp
            enddo
         endif
      enddo
!
!        copy the triangular factor of the qr factorization into r.
!
      sing = .false.
      do j = 1 , n
         l = j
         jm1 = j - 1
         if ( jm1>=1 ) then
            do i = 1 , jm1
               r(l) = Fjac(i,j)
               l = l + n - i
            enddo
         endif
         r(l) = Wa1(j)
         if ( Wa1(j)==zero ) sing = .true.
      enddo
!
!        accumulate the orthogonal factor in fjac.
!
      call qform(n,n,Fjac,Ldfjac,Wa1)
!
!        rescale if necessary.
!
      if ( Mode/=2 ) then
         do j = 1 , n
            Diag(j) = dmax1(Diag(j),Wa2(j))
         enddo
      endif
!
!        beginning of the inner loop.
!
!
!           if requested, call fcn to enable printing of iterates.
!
 200  if ( Nprint>0 ) then
         iflag = 0
         if ( mod(iter-1,Nprint)==0 ) call fcn(n,x,Fvec,iflag)
         if ( iflag<0 ) goto 300
      endif
!
!           determine the direction p.
!
      call dogleg(n,r,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
      do j = 1 , n
         Wa1(j) = -Wa1(j)
         Wa2(j) = x(j) + Wa1(j)
         Wa3(j) = Diag(j)*Wa1(j)
      enddo
      pnorm = enorm(n,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      if ( iter==1 ) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      call fcn(n,Wa2,Wa4,iflag)
      Nfev = Nfev + 1
      if ( iflag>=0 ) then
         fnorm1 = enorm(n,Wa4)
!
!           compute the scaled actual reduction.
!
         actred = -one
         if ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
         l = 1
         do i = 1 , n
            sum = zero
            do j = i , n
               sum = sum + r(l)*Wa1(j)
               l = l + 1
            enddo
            Wa3(i) = Qtf(i) + sum
         enddo
         temp = enorm(n,Wa3)
         prered = zero
         if ( temp<fnorm ) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
         ratio = zero
         if ( prered>zero ) ratio = actred/prered
!
!           update the step bound.
!
         if ( ratio>=p1 ) then
            ncfail = 0
            ncsuc = ncsuc + 1
            if ( ratio>=p5 .or. ncsuc>1 ) delta = dmax1(delta,pnorm/p5)
            if ( dabs(ratio-one)<=p1 ) delta = pnorm/p5
         else
            ncsuc = 0
            ncfail = ncfail + 1
            delta = p5*delta
         endif
!
!           test for successful iteration.
!
         if ( ratio>=p0001 ) then
!
!           successful iteration. update x, fvec, and their norms.
!
            do j = 1 , n
               x(j) = Wa2(j)
               Wa2(j) = Diag(j)*x(j)
               Fvec(j) = Wa4(j)
            enddo
            xnorm = enorm(n,Wa2)
            fnorm = fnorm1
            iter = iter + 1
         endif
!
!           determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         if ( actred>=p001 ) nslow1 = 0
         if ( jeval ) nslow2 = nslow2 + 1
         if ( actred>=p1 ) nslow2 = 0
!
!           test for convergence.
!
         if ( delta<=Xtol*xnorm .or. fnorm==zero ) Info = 1
         if ( Info==0 ) then
!
!           tests for termination and stringent tolerances.
!
            if ( Nfev>=Maxfev ) Info = 2
            if ( p1*dmax1(p1*delta,pnorm)<=epsmch*xnorm ) Info = 3
            if ( nslow2==5 ) Info = 4
            if ( nslow1==10 ) Info = 5
            if ( Info==0 ) then
!
!           criterion for recalculating jacobian approximation
!           by forward differences.
!
               if ( ncfail==2 ) goto 100
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
               do j = 1 , n
                  sum = zero
                  do i = 1 , n
                     sum = sum + Fjac(i,j)*Wa4(i)
                  enddo
                  Wa2(j) = (sum-Wa3(j))/pnorm
                  Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
                  if ( ratio>=p0001 ) Qtf(j) = sum
               enddo
!
!           compute the qr factorization of the updated jacobian.
!
               call r1updt(n,n,r,Lr,Wa1,Wa2,Wa3,sing)
               call r1mpyq(n,n,Fjac,Ldfjac,Wa2,Wa3)
               call r1mpyq(1,n,Qtf,1,Wa2,Wa3)
!
!           end of the inner loop.
!
               jeval = .false.
!
!        end of the outer loop.
!
               goto 200
            endif
         endif
      endif
!
!     termination, either normal or user imposed.
!
 300  if ( iflag<0 ) Info = iflag
      iflag = 0
      if ( Nprint>0 ) call fcn(n,x,Fvec,iflag)
!
!     last card of subroutine hybrd.
!
      end

      subroutine hybrd1(fcn,n,x,Fvec,Tol,Info,Wa,Lwa)
      implicit none

      integer n , Info , Lwa
      double precision Tol
      double precision x(n) , Fvec(n) , Wa(Lwa)
      external fcn
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
      integer index , j , lr , maxfev , ml , mode , mu , nfev , nprint
      double precision epsfcn , factor , one , xtol , zero
      data factor , one , zero/1.0d2 , 1.0d0 , 0.0d0/
      Info = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. Tol>=zero .and. Lwa>=(n*(3*n+13))/2 ) then
!
!     call hybrd.
!
         maxfev = 200*(n+1)
         xtol = Tol
         ml = n - 1
         mu = n - 1
         epsfcn = zero
         mode = 2
         do j = 1 , n
            Wa(j) = one
         enddo
         nprint = 0
         lr = (n*(n+1))/2
         index = 6*n + lr
         call hybrd(fcn,n,x,Fvec,xtol,maxfev,ml,mu,epsfcn,Wa(1),mode,   &
                  & factor,nprint,Info,nfev,Wa(index+1),n,Wa(6*n+1),lr, &
                  & Wa(n+1),Wa(2*n+1),Wa(3*n+1),Wa(4*n+1),Wa(5*n+1))
         if ( Info==5 ) Info = 4
      endif
!
!     last card of subroutine hybrd1.
!
      end

      subroutine hybrj(fcn,n,x,Fvec,Fjac,Ldfjac,Xtol,Maxfev,Diag,Mode,  &
                     & Factor,Nprint,Info,Nfev,Njev,r,Lr,Qtf,Wa1,Wa2,   &
                     & Wa3,Wa4)
      implicit none

      integer n , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev , Njev ,&
            & Lr
      double precision Xtol , Factor
      double precision x(n) , Fvec(n) , Fjac(Ldfjac,n) , Diag(n) ,      &
                     & r(Lr) , Qtf(n) , Wa1(n) , Wa2(n) , Wa3(n) ,      &
                     & Wa4(n)
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
      integer i , iflag , iter , j , jm1 , l , ncfail , ncsuc , nslow1 ,&
            & nslow2
      integer iwa(1)
      logical jeval , sing
      double precision actred , delta , epsmch , fnorm , fnorm1 , one , &
                     & pnorm , prered , p1 , p5 , p001 , p0001 , ratio ,&
                     & sum , temp , xnorm , zero
      double precision dpmpar , enorm
      data one , p1 , p5 , p001 , p0001 , zero/1.0d0 , 1.0d-1 , 5.0d-1 ,&
         & 1.0d-3 , 1.0d-4 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
      Njev = 0
!
!     check the input parameters for errors.
!
      if ( n<=0 .or. Ldfjac<n .or. Xtol<zero .or. Maxfev<=0 .or.        &
         & Factor<=zero .or. Lr<(n*(n+1))/2 ) goto 300
      if ( Mode==2 ) then
         do j = 1 , n
            if ( Diag(j)<=zero ) goto 300
         enddo
      endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,Fvec,Fjac,Ldfjac,iflag)
      Nfev = 1
      if ( iflag<0 ) goto 300
      fnorm = enorm(n,Fvec)
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
 100  jeval = .true.
!
!        calculate the jacobian matrix.
!
      iflag = 2
      call fcn(n,x,Fvec,Fjac,Ldfjac,iflag)
      Njev = Njev + 1
      if ( iflag<0 ) goto 300
!
!        compute the qr factorization of the jacobian.
!
      call qrfac(n,n,Fjac,Ldfjac,.false.,iwa,1,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      if ( iter==1 ) then
         if ( Mode/=2 ) then
            do j = 1 , n
               Diag(j) = Wa2(j)
               if ( Wa2(j)==zero ) Diag(j) = one
            enddo
         endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do j = 1 , n
            Wa3(j) = Diag(j)*x(j)
         enddo
         xnorm = enorm(n,Wa3)
         delta = Factor*xnorm
         if ( delta==zero ) delta = Factor
      endif
!
!        form (q transpose)*fvec and store in qtf.
!
      do i = 1 , n
         Qtf(i) = Fvec(i)
      enddo
      do j = 1 , n
         if ( Fjac(j,j)/=zero ) then
            sum = zero
            do i = j , n
               sum = sum + Fjac(i,j)*Qtf(i)
            enddo
            temp = -sum/Fjac(j,j)
            do i = j , n
               Qtf(i) = Qtf(i) + Fjac(i,j)*temp
            enddo
         endif
      enddo
!
!        copy the triangular factor of the qr factorization into r.
!
      sing = .false.
      do j = 1 , n
         l = j
         jm1 = j - 1
         if ( jm1>=1 ) then
            do i = 1 , jm1
               r(l) = Fjac(i,j)
               l = l + n - i
            enddo
         endif
         r(l) = Wa1(j)
         if ( Wa1(j)==zero ) sing = .true.
      enddo
!
!        accumulate the orthogonal factor in fjac.
!
      call qform(n,n,Fjac,Ldfjac,Wa1)
!
!        rescale if necessary.
!
      if ( Mode/=2 ) then
         do j = 1 , n
            Diag(j) = dmax1(Diag(j),Wa2(j))
         enddo
      endif
!
!        beginning of the inner loop.
!
!
!           if requested, call fcn to enable printing of iterates.
!
 200  if ( Nprint>0 ) then
         iflag = 0
         if ( mod(iter-1,Nprint)==0 )                                   &
            & call fcn(n,x,Fvec,Fjac,Ldfjac,iflag)
         if ( iflag<0 ) goto 300
      endif
!
!           determine the direction p.
!
      call dogleg(n,r,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
      do j = 1 , n
         Wa1(j) = -Wa1(j)
         Wa2(j) = x(j) + Wa1(j)
         Wa3(j) = Diag(j)*Wa1(j)
      enddo
      pnorm = enorm(n,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      if ( iter==1 ) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      call fcn(n,Wa2,Wa4,Fjac,Ldfjac,iflag)
      Nfev = Nfev + 1
      if ( iflag>=0 ) then
         fnorm1 = enorm(n,Wa4)
!
!           compute the scaled actual reduction.
!
         actred = -one
         if ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
         l = 1
         do i = 1 , n
            sum = zero
            do j = i , n
               sum = sum + r(l)*Wa1(j)
               l = l + 1
            enddo
            Wa3(i) = Qtf(i) + sum
         enddo
         temp = enorm(n,Wa3)
         prered = zero
         if ( temp<fnorm ) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
         ratio = zero
         if ( prered>zero ) ratio = actred/prered
!
!           update the step bound.
!
         if ( ratio>=p1 ) then
            ncfail = 0
            ncsuc = ncsuc + 1
            if ( ratio>=p5 .or. ncsuc>1 ) delta = dmax1(delta,pnorm/p5)
            if ( dabs(ratio-one)<=p1 ) delta = pnorm/p5
         else
            ncsuc = 0
            ncfail = ncfail + 1
            delta = p5*delta
         endif
!
!           test for successful iteration.
!
         if ( ratio>=p0001 ) then
!
!           successful iteration. update x, fvec, and their norms.
!
            do j = 1 , n
               x(j) = Wa2(j)
               Wa2(j) = Diag(j)*x(j)
               Fvec(j) = Wa4(j)
            enddo
            xnorm = enorm(n,Wa2)
            fnorm = fnorm1
            iter = iter + 1
         endif
!
!           determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         if ( actred>=p001 ) nslow1 = 0
         if ( jeval ) nslow2 = nslow2 + 1
         if ( actred>=p1 ) nslow2 = 0
!
!           test for convergence.
!
         if ( delta<=Xtol*xnorm .or. fnorm==zero ) Info = 1
         if ( Info==0 ) then
!
!           tests for termination and stringent tolerances.
!
            if ( Nfev>=Maxfev ) Info = 2
            if ( p1*dmax1(p1*delta,pnorm)<=epsmch*xnorm ) Info = 3
            if ( nslow2==5 ) Info = 4
            if ( nslow1==10 ) Info = 5
            if ( Info==0 ) then
!
!           criterion for recalculating jacobian.
!
               if ( ncfail==2 ) goto 100
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
               do j = 1 , n
                  sum = zero
                  do i = 1 , n
                     sum = sum + Fjac(i,j)*Wa4(i)
                  enddo
                  Wa2(j) = (sum-Wa3(j))/pnorm
                  Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
                  if ( ratio>=p0001 ) Qtf(j) = sum
               enddo
!
!           compute the qr factorization of the updated jacobian.
!
               call r1updt(n,n,r,Lr,Wa1,Wa2,Wa3,sing)
               call r1mpyq(n,n,Fjac,Ldfjac,Wa2,Wa3)
               call r1mpyq(1,n,Qtf,1,Wa2,Wa3)
!
!           end of the inner loop.
!
               jeval = .false.
!
!        end of the outer loop.
!
               goto 200
            endif
         endif
      endif
!
!     termination, either normal or user imposed.
!
 300  if ( iflag<0 ) Info = iflag
      iflag = 0
      if ( Nprint>0 ) call fcn(n,x,Fvec,Fjac,Ldfjac,iflag)
!
!     last card of subroutine hybrj.
!
      end

      subroutine hybrj1(fcn,n,x,Fvec,Fjac,Ldfjac,Tol,Info,Wa,Lwa)
      implicit none

      integer n , Ldfjac , Info , Lwa
      double precision Tol
      double precision x(n) , Fvec(n) , Fjac(Ldfjac,n) , Wa(Lwa)
      external fcn
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
      integer j , lr , maxfev , mode , nfev , njev , nprint
      double precision factor , one , xtol , zero
      data factor , one , zero/1.0d2 , 1.0d0 , 0.0d0/
      Info = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. Ldfjac>=n .and. Tol>=zero .and. Lwa>=(n*(n+13))/2 )&
         & then
!
!     call hybrj.
!
         maxfev = 100*(n+1)
         xtol = Tol
         mode = 2
         do j = 1 , n
            Wa(j) = one
         enddo
         nprint = 0
         lr = (n*(n+1))/2
         call hybrj(fcn,n,x,Fvec,Fjac,Ldfjac,xtol,maxfev,Wa(1),mode,    &
                  & factor,nprint,Info,nfev,njev,Wa(6*n+1),lr,Wa(n+1),  &
                  & Wa(2*n+1),Wa(3*n+1),Wa(4*n+1),Wa(5*n+1))
         if ( Info==5 ) Info = 4
      endif
!
!     last card of subroutine hybrj1.
!
      end

      subroutine lmder(fcn,m,n,x,Fvec,Fjac,Ldfjac,Ftol,Xtol,Gtol,Maxfev,&
                     & Diag,Mode,Factor,Nprint,Info,Nfev,Njev,Ipvt,Qtf, &
                     & Wa1,Wa2,Wa3,Wa4)
      implicit none

      integer m , n , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev ,   &
            & Njev
      integer Ipvt(n)
      double precision Ftol , Xtol , Gtol , Factor
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n) , Diag(n) ,      &
                     & Qtf(n) , Wa1(n) , Wa2(n) , Wa3(n) , Wa4(m)
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
      integer i , iflag , iter , j , l
      double precision actred , delta , dirder , epsmch , fnorm ,       &
                     & fnorm1 , gnorm , one , par , pnorm , prered ,    &
                     & p1 , p5 , p25 , p75 , p0001 , ratio , sum ,      &
                     & temp , temp1 , temp2 , xnorm , zero
      double precision dpmpar , enorm
      data one , p1 , p5 , p25 , p75 , p0001 , zero/1.0d0 , 1.0d-1 ,    &
         & 5.0d-1 , 2.5d-1 , 7.5d-1 , 1.0d-4 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
      Njev = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. m>=n .and. Ldfjac>=m .and. Ftol>=zero .and.        &
         & Xtol>=zero .and. Gtol>=zero .and. Maxfev>0 .and.             &
         & Factor>zero ) then
         if ( Mode==2 ) then
            do j = 1 , n
               if ( Diag(j)<=zero ) goto 100
            enddo
         endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
         iflag = 1
         call fcn(m,n,x,Fvec,Fjac,Ldfjac,iflag)
         Nfev = 1
         if ( iflag>=0 ) then
            fnorm = enorm(m,Fvec)
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
            call fcn(m,n,x,Fvec,Fjac,Ldfjac,iflag)
            Njev = Njev + 1
            if ( iflag>=0 ) then
!
!        if requested, call fcn to enable printing of iterates.
!
               if ( Nprint>0 ) then
                  iflag = 0
                  if ( mod(iter-1,Nprint)==0 )                          &
                     & call fcn(m,n,x,Fvec,Fjac,Ldfjac,iflag)
                  if ( iflag<0 ) goto 100
               endif
!
!        compute the qr factorization of the jacobian.
!
               call qrfac(m,n,Fjac,Ldfjac,.true.,Ipvt,n,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
               if ( iter==1 ) then
                  if ( Mode/=2 ) then
                     do j = 1 , n
                        Diag(j) = Wa2(j)
                        if ( Wa2(j)==zero ) Diag(j) = one
                     enddo
                  endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
                  do j = 1 , n
                     Wa3(j) = Diag(j)*x(j)
                  enddo
                  xnorm = enorm(n,Wa3)
                  delta = Factor*xnorm
                  if ( delta==zero ) delta = Factor
               endif
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
               do i = 1 , m
                  Wa4(i) = Fvec(i)
               enddo
               do j = 1 , n
                  if ( Fjac(j,j)/=zero ) then
                     sum = zero
                     do i = j , m
                        sum = sum + Fjac(i,j)*Wa4(i)
                     enddo
                     temp = -sum/Fjac(j,j)
                     do i = j , m
                        Wa4(i) = Wa4(i) + Fjac(i,j)*temp
                     enddo
                  endif
                  Fjac(j,j) = Wa1(j)
                  Qtf(j) = Wa4(j)
               enddo
!
!        compute the norm of the scaled gradient.
!
               gnorm = zero
               if ( fnorm/=zero ) then
                  do j = 1 , n
                     l = Ipvt(j)
                     if ( Wa2(l)/=zero ) then
                        sum = zero
                        do i = 1 , j
                           sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
                        enddo
                        gnorm = dmax1(gnorm,dabs(sum/Wa2(l)))
                     endif
                  enddo
               endif
!
!        test for convergence of the gradient norm.
!
               if ( gnorm<=Gtol ) Info = 4
               if ( Info==0 ) then
!
!        rescale if necessary.
!
                  if ( Mode/=2 ) then
                     do j = 1 , n
                        Diag(j) = dmax1(Diag(j),Wa2(j))
                     enddo
                  endif
!
!        beginning of the inner loop.
!
!
!           determine the levenberg-marquardt parameter.
!
 25               call lmpar(n,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1, &
                           & Wa2,Wa3,Wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
                  do j = 1 , n
                     Wa1(j) = -Wa1(j)
                     Wa2(j) = x(j) + Wa1(j)
                     Wa3(j) = Diag(j)*Wa1(j)
                  enddo
                  pnorm = enorm(n,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
                  if ( iter==1 ) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
                  iflag = 1
                  call fcn(m,n,Wa2,Wa4,Fjac,Ldfjac,iflag)
                  Nfev = Nfev + 1
                  if ( iflag>=0 ) then
                     fnorm1 = enorm(m,Wa4)
!
!           compute the scaled actual reduction.
!
                     actred = -one
                     if ( p1*fnorm1<fnorm ) actred = one -              &
                        & (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
                     do j = 1 , n
                        Wa3(j) = zero
                        l = Ipvt(j)
                        temp = Wa1(l)
                        do i = 1 , j
                           Wa3(i) = Wa3(i) + Fjac(i,j)*temp
                        enddo
                     enddo
                     temp1 = enorm(n,Wa3)/fnorm
                     temp2 = (dsqrt(par)*pnorm)/fnorm
                     prered = temp1**2 + temp2**2/p5
                     dirder = -(temp1**2+temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
                     ratio = zero
                     if ( prered/=zero ) ratio = actred/prered
!
!           update the step bound.
!
                     if ( ratio<=p25 ) then
                        if ( actred>=zero ) temp = p5
                        if ( actred<zero )                              &
                           & temp = p5*dirder/(dirder+p5*actred)
                        if ( p1*fnorm1>=fnorm .or. temp<p1 ) temp = p1
                        delta = temp*dmin1(delta,pnorm/p1)
                        par = par/temp
                     elseif ( par==zero .or. ratio>=p75 ) then
                        delta = pnorm/p5
                        par = p5*par
                     endif
!
!           test for successful iteration.
!
                     if ( ratio>=p0001 ) then
!
!           successful iteration. update x, fvec, and their norms.
!
                        do j = 1 , n
                           x(j) = Wa2(j)
                           Wa2(j) = Diag(j)*x(j)
                        enddo
                        do i = 1 , m
                           Fvec(i) = Wa4(i)
                        enddo
                        xnorm = enorm(n,Wa2)
                        fnorm = fnorm1
                        iter = iter + 1
                     endif
!
!           tests for convergence.
!
                     if ( dabs(actred)<=Ftol .and. prered<=Ftol .and.   &
                        & p5*ratio<=one ) Info = 1
                     if ( delta<=Xtol*xnorm ) Info = 2
                     if ( dabs(actred)<=Ftol .and. prered<=Ftol .and.   &
                        & p5*ratio<=one .and. Info==2 ) Info = 3
                     if ( Info==0 ) then
!
!           tests for termination and stringent tolerances.
!
                        if ( Nfev>=Maxfev ) Info = 5
                        if ( dabs(actred)<=epsmch .and.                 &
                           & prered<=epsmch .and. p5*ratio<=one )       &
                           & Info = 6
                        if ( delta<=epsmch*xnorm ) Info = 7
                        if ( gnorm<=epsmch ) Info = 8
                        if ( Info==0 ) then
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
!
!        end of the outer loop.
!
                           if ( ratio>=p0001 ) goto 20
                           goto 25
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!
!     termination, either normal or user imposed.
!
 100  if ( iflag<0 ) Info = iflag
      iflag = 0
      if ( Nprint>0 ) call fcn(m,n,x,Fvec,Fjac,Ldfjac,iflag)
!
!     last card of subroutine lmder.
!
      end

      subroutine lmder1(fcn,m,n,x,Fvec,Fjac,Ldfjac,Tol,Info,Ipvt,Wa,Lwa)
      implicit none

      integer m , n , Ldfjac , Info , Lwa
      integer Ipvt(n)
      double precision Tol
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n) , Wa(Lwa)
      external fcn
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
      integer maxfev , mode , nfev , njev , nprint
      double precision factor , ftol , gtol , xtol , zero
      data factor , zero/1.0d2 , 0.0d0/
      Info = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. m>=n .and. Ldfjac>=m .and. Tol>=zero .and.         &
         & Lwa>=5*n+m ) then
!
!     call lmder.
!
         maxfev = 100*(n+1)
         ftol = Tol
         xtol = Tol
         gtol = zero
         mode = 1
         nprint = 0
         call lmder(fcn,m,n,x,Fvec,Fjac,Ldfjac,ftol,xtol,gtol,maxfev,   &
                  & Wa(1),mode,factor,nprint,Info,nfev,njev,Ipvt,Wa(n+1)&
                  & ,Wa(2*n+1),Wa(3*n+1),Wa(4*n+1),Wa(5*n+1))
         if ( Info==8 ) Info = 4
      endif
!
!     last card of subroutine lmder1.
!
      end

      subroutine lmdif(fcn,m,n,x,Fvec,Ftol,Xtol,Gtol,Maxfev,Epsfcn,Diag,&
                     & Mode,Factor,Nprint,Info,Nfev,Fjac,Ldfjac,Ipvt,   &
                     & Qtf,Wa1,Wa2,Wa3,Wa4)
      implicit none

      integer m , n , Maxfev , Mode , Nprint , Info , Nfev , Ldfjac
      integer Ipvt(n)
      double precision Ftol , Xtol , Gtol , Epsfcn , Factor
      double precision x(n) , Fvec(m) , Diag(n) , Fjac(Ldfjac,n) ,      &
                     & Qtf(n) , Wa1(n) , Wa2(n) , Wa3(n) , Wa4(m)
      external fcn
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
      integer i , iflag , iter , j , l
      double precision actred , delta , dirder , epsmch , fnorm ,       &
                     & fnorm1 , gnorm , one , par , pnorm , prered ,    &
                     & p1 , p5 , p25 , p75 , p0001 , ratio , sum ,      &
                     & temp , temp1 , temp2 , xnorm , zero
      double precision dpmpar , enorm
      data one , p1 , p5 , p25 , p75 , p0001 , zero/1.0d0 , 1.0d-1 ,    &
         & 5.0d-1 , 2.5d-1 , 7.5d-1 , 1.0d-4 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. m>=n .and. Ldfjac>=m .and. Ftol>=zero .and.        &
         & Xtol>=zero .and. Gtol>=zero .and. Maxfev>0 .and.             &
         & Factor>zero ) then
         if ( Mode==2 ) then
            do j = 1 , n
               if ( Diag(j)<=zero ) goto 100
            enddo
         endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
         iflag = 1
         call fcn(m,n,x,Fvec,iflag)
         Nfev = 1
         if ( iflag>=0 ) then
            fnorm = enorm(m,Fvec)
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
            call fdjac2(fcn,m,n,x,Fvec,Fjac,Ldfjac,iflag,Epsfcn,Wa4)
            Nfev = Nfev + n
            if ( iflag>=0 ) then
!
!        if requested, call fcn to enable printing of iterates.
!
               if ( Nprint>0 ) then
                  iflag = 0
                  if ( mod(iter-1,Nprint)==0 )                          &
                     & call fcn(m,n,x,Fvec,iflag)
                  if ( iflag<0 ) goto 100
               endif
!
!        compute the qr factorization of the jacobian.
!
               call qrfac(m,n,Fjac,Ldfjac,.true.,Ipvt,n,Wa1,Wa2,Wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
               if ( iter==1 ) then
                  if ( Mode/=2 ) then
                     do j = 1 , n
                        Diag(j) = Wa2(j)
                        if ( Wa2(j)==zero ) Diag(j) = one
                     enddo
                  endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
                  do j = 1 , n
                     Wa3(j) = Diag(j)*x(j)
                  enddo
                  xnorm = enorm(n,Wa3)
                  delta = Factor*xnorm
                  if ( delta==zero ) delta = Factor
               endif
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
               do i = 1 , m
                  Wa4(i) = Fvec(i)
               enddo
               do j = 1 , n
                  if ( Fjac(j,j)/=zero ) then
                     sum = zero
                     do i = j , m
                        sum = sum + Fjac(i,j)*Wa4(i)
                     enddo
                     temp = -sum/Fjac(j,j)
                     do i = j , m
                        Wa4(i) = Wa4(i) + Fjac(i,j)*temp
                     enddo
                  endif
                  Fjac(j,j) = Wa1(j)
                  Qtf(j) = Wa4(j)
               enddo
!
!        compute the norm of the scaled gradient.
!
               gnorm = zero
               if ( fnorm/=zero ) then
                  do j = 1 , n
                     l = Ipvt(j)
                     if ( Wa2(l)/=zero ) then
                        sum = zero
                        do i = 1 , j
                           sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
                        enddo
                        gnorm = dmax1(gnorm,dabs(sum/Wa2(l)))
                     endif
                  enddo
               endif
!
!        test for convergence of the gradient norm.
!
               if ( gnorm<=Gtol ) Info = 4
               if ( Info==0 ) then
!
!        rescale if necessary.
!
                  if ( Mode/=2 ) then
                     do j = 1 , n
                        Diag(j) = dmax1(Diag(j),Wa2(j))
                     enddo
                  endif
!
!        beginning of the inner loop.
!
!
!           determine the levenberg-marquardt parameter.
!
 25               call lmpar(n,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1, &
                           & Wa2,Wa3,Wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
                  do j = 1 , n
                     Wa1(j) = -Wa1(j)
                     Wa2(j) = x(j) + Wa1(j)
                     Wa3(j) = Diag(j)*Wa1(j)
                  enddo
                  pnorm = enorm(n,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
                  if ( iter==1 ) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
                  iflag = 1
                  call fcn(m,n,Wa2,Wa4,iflag)
                  Nfev = Nfev + 1
                  if ( iflag>=0 ) then
                     fnorm1 = enorm(m,Wa4)
!
!           compute the scaled actual reduction.
!
                     actred = -one
                     if ( p1*fnorm1<fnorm ) actred = one -              &
                        & (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
                     do j = 1 , n
                        Wa3(j) = zero
                        l = Ipvt(j)
                        temp = Wa1(l)
                        do i = 1 , j
                           Wa3(i) = Wa3(i) + Fjac(i,j)*temp
                        enddo
                     enddo
                     temp1 = enorm(n,Wa3)/fnorm
                     temp2 = (dsqrt(par)*pnorm)/fnorm
                     prered = temp1**2 + temp2**2/p5
                     dirder = -(temp1**2+temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
                     ratio = zero
                     if ( prered/=zero ) ratio = actred/prered
!
!           update the step bound.
!
                     if ( ratio<=p25 ) then
                        if ( actred>=zero ) temp = p5
                        if ( actred<zero )                              &
                           & temp = p5*dirder/(dirder+p5*actred)
                        if ( p1*fnorm1>=fnorm .or. temp<p1 ) temp = p1
                        delta = temp*dmin1(delta,pnorm/p1)
                        par = par/temp
                     elseif ( par==zero .or. ratio>=p75 ) then
                        delta = pnorm/p5
                        par = p5*par
                     endif
!
!           test for successful iteration.
!
                     if ( ratio>=p0001 ) then
!
!           successful iteration. update x, fvec, and their norms.
!
                        do j = 1 , n
                           x(j) = Wa2(j)
                           Wa2(j) = Diag(j)*x(j)
                        enddo
                        do i = 1 , m
                           Fvec(i) = Wa4(i)
                        enddo
                        xnorm = enorm(n,Wa2)
                        fnorm = fnorm1
                        iter = iter + 1
                     endif
!
!           tests for convergence.
!
                     if ( dabs(actred)<=Ftol .and. prered<=Ftol .and.   &
                        & p5*ratio<=one ) Info = 1
                     if ( delta<=Xtol*xnorm ) Info = 2
                     if ( dabs(actred)<=Ftol .and. prered<=Ftol .and.   &
                        & p5*ratio<=one .and. Info==2 ) Info = 3
                     if ( Info==0 ) then
!
!           tests for termination and stringent tolerances.
!
                        if ( Nfev>=Maxfev ) Info = 5
                        if ( dabs(actred)<=epsmch .and.                 &
                           & prered<=epsmch .and. p5*ratio<=one )       &
                           & Info = 6
                        if ( delta<=epsmch*xnorm ) Info = 7
                        if ( gnorm<=epsmch ) Info = 8
                        if ( Info==0 ) then
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
!
!        end of the outer loop.
!
                           if ( ratio>=p0001 ) goto 20
                           goto 25
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!
!     termination, either normal or user imposed.
!
 100  if ( iflag<0 ) Info = iflag
      iflag = 0
      if ( Nprint>0 ) call fcn(m,n,x,Fvec,iflag)
!
!     last card of subroutine lmdif.
!
      end

      subroutine lmdif1(fcn,m,n,x,Fvec,Tol,Info,Iwa,Wa,Lwa)
      implicit none

      integer m , n , Info , Lwa
      integer Iwa(n)
      double precision Tol
      double precision x(n) , Fvec(m) , Wa(Lwa)
      external fcn
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
      integer maxfev , mode , mp5n , nfev , nprint
      double precision epsfcn , factor , ftol , gtol , xtol , zero
      data factor , zero/1.0d2 , 0.0d0/
      Info = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. m>=n .and. Tol>=zero .and. Lwa>=m*n+5*n+m ) then
!
!     call lmdif.
!
         maxfev = 200*(n+1)
         ftol = Tol
         xtol = Tol
         gtol = zero
         epsfcn = zero
         mode = 1
         nprint = 0
         mp5n = m + 5*n
         call lmdif(fcn,m,n,x,Fvec,ftol,xtol,gtol,maxfev,epsfcn,Wa(1),  &
                  & mode,factor,nprint,Info,nfev,Wa(mp5n+1),m,Iwa,      &
                  & Wa(n+1),Wa(2*n+1),Wa(3*n+1),Wa(4*n+1),Wa(5*n+1))
         if ( Info==8 ) Info = 4
      endif
!
!     last card of subroutine lmdif1.
!
      end

      subroutine lmpar(n,r,Ldr,Ipvt,Diag,Qtb,Delta,Par,x,Sdiag,Wa1,Wa2)
      implicit none

      integer n , Ldr
      integer Ipvt(n)
      double precision Delta , Par
      double precision r(Ldr,n) , Diag(n) , Qtb(n) , x(n) , Sdiag(n) ,  &
                     & Wa1(n) , Wa2(n)
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
      integer i , iter , j , jm1 , jp1 , k , l , nsing
      double precision dxnorm , dwarf , fp , gnorm , parc , parl ,      &
                     & paru , p1 , p001 , sum , temp , zero
      double precision dpmpar , enorm
      data p1 , p001 , zero/1.0d-1 , 1.0d-3 , 0.0d0/
!
!     dwarf is the smallest positive magnitude.
!
      dwarf = dpmpar(2)
!
!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.
!
      nsing = n
      do j = 1 , n
         Wa1(j) = Qtb(j)
         if ( r(j,j)==zero .and. nsing==n ) nsing = j - 1
         if ( nsing<n ) Wa1(j) = zero
      enddo
      if ( nsing>=1 ) then
         do k = 1 , nsing
            j = nsing - k + 1
            Wa1(j) = Wa1(j)/r(j,j)
            temp = Wa1(j)
            jm1 = j - 1
            if ( jm1>=1 ) then
               do i = 1 , jm1
                  Wa1(i) = Wa1(i) - r(i,j)*temp
               enddo
            endif
         enddo
      endif
      do j = 1 , n
         l = Ipvt(j)
         x(l) = Wa1(j)
      enddo
!
!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.
!
      iter = 0
      do j = 1 , n
         Wa2(j) = Diag(j)*x(j)
      enddo
      dxnorm = enorm(n,Wa2)
      fp = dxnorm - Delta
      if ( fp<=p1*Delta ) then
!
!     termination.
!
         if ( iter==0 ) Par = zero
      else
!
!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function. otherwise set this bound to zero.
!
         parl = zero
         if ( nsing>=n ) then
            do j = 1 , n
               l = Ipvt(j)
               Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
            enddo
            do j = 1 , n
               sum = zero
               jm1 = j - 1
               if ( jm1>=1 ) then
                  do i = 1 , jm1
                     sum = sum + r(i,j)*Wa1(i)
                  enddo
               endif
               Wa1(j) = (Wa1(j)-sum)/r(j,j)
            enddo
            temp = enorm(n,Wa1)
            parl = ((fp/Delta)/temp)/temp
         endif
!
!     calculate an upper bound, paru, for the zero of the function.
!
         do j = 1 , n
            sum = zero
            do i = 1 , j
               sum = sum + r(i,j)*Qtb(i)
            enddo
            l = Ipvt(j)
            Wa1(j) = sum/Diag(l)
         enddo
         gnorm = enorm(n,Wa1)
         paru = gnorm/Delta
         if ( paru==zero ) paru = dwarf/dmin1(Delta,p1)
!
!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.
!
         Par = dmax1(Par,parl)
         Par = dmin1(Par,paru)
         if ( Par==zero ) Par = gnorm/dxnorm
!
!     beginning of an iteration.
!
 50      iter = iter + 1
!
!        evaluate the function at the current value of par.
!
         if ( Par==zero ) Par = dmax1(dwarf,p001*paru)
         temp = dsqrt(Par)
         do j = 1 , n
            Wa1(j) = temp*Diag(j)
         enddo
         call qrsolv(n,r,Ldr,Ipvt,Wa1,Qtb,x,Sdiag,Wa2)
         do j = 1 , n
            Wa2(j) = Diag(j)*x(j)
         enddo
         dxnorm = enorm(n,Wa2)
         temp = fp
         fp = dxnorm - Delta
!
!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.
!
         if ( dabs(fp)<=p1*Delta .or. parl==zero .and. fp<=temp .and.   &
            & temp<zero .or. iter==10 ) then
            if ( iter==0 ) Par = zero
         else
!
!        compute the newton correction.
!
            do j = 1 , n
               l = Ipvt(j)
               Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
            enddo
            do j = 1 , n
               Wa1(j) = Wa1(j)/Sdiag(j)
               temp = Wa1(j)
               jp1 = j + 1
               if ( n>=jp1 ) then
                  do i = jp1 , n
                     Wa1(i) = Wa1(i) - r(i,j)*temp
                  enddo
               endif
            enddo
            temp = enorm(n,Wa1)
            parc = ((fp/Delta)/temp)/temp
!
!        depending on the sign of the function, update parl or paru.
!
            if ( fp>zero ) parl = dmax1(parl,Par)
            if ( fp<zero ) paru = dmin1(paru,Par)
!
!        compute an improved estimate for par.
!
            Par = dmax1(parl,Par+parc)
!
!        end of an iteration.
!
            goto 50
         endif
      endif
!
!     last card of subroutine lmpar.
!
      end

      subroutine lmstr(fcn,m,n,x,Fvec,Fjac,Ldfjac,Ftol,Xtol,Gtol,Maxfev,&
                     & Diag,Mode,Factor,Nprint,Info,Nfev,Njev,Ipvt,Qtf, &
                     & Wa1,Wa2,Wa3,Wa4)
      implicit none

      integer m , n , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev ,   &
            & Njev
      integer Ipvt(n)
      logical sing
      double precision Ftol , Xtol , Gtol , Factor
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n) , Diag(n) ,      &
                     & Qtf(n) , Wa1(n) , Wa2(n) , Wa3(n) , Wa4(m)
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
      integer i , iflag , iter , j , l
      double precision actred , delta , dirder , epsmch , fnorm ,       &
                     & fnorm1 , gnorm , one , par , pnorm , prered ,    &
                     & p1 , p5 , p25 , p75 , p0001 , ratio , sum ,      &
                     & temp , temp1 , temp2 , xnorm , zero
      double precision dpmpar , enorm
      data one , p1 , p5 , p25 , p75 , p0001 , zero/1.0d0 , 1.0d-1 ,    &
         & 5.0d-1 , 2.5d-1 , 7.5d-1 , 1.0d-4 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      Info = 0
      iflag = 0
      Nfev = 0
      Njev = 0
!
!     check the input parameters for errors.
!
      if ( n<=0 .or. m<n .or. Ldfjac<n .or. Ftol<zero .or.              &
         & Xtol<zero .or. Gtol<zero .or. Maxfev<=0 .or. Factor<=zero )  &
         & goto 200
      if ( Mode==2 ) then
         do j = 1 , n
            if ( Diag(j)<=zero ) goto 200
         enddo
      endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(m,n,x,Fvec,Wa3,iflag)
      Nfev = 1
      if ( iflag<0 ) goto 200
      fnorm = enorm(m,Fvec)
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
 100  if ( Nprint>0 ) then
         iflag = 0
         if ( mod(iter-1,Nprint)==0 ) call fcn(m,n,x,Fvec,Wa3,iflag)
         if ( iflag<0 ) goto 200
      endif
!
!        compute the qr factorization of the jacobian matrix
!        calculated one row at a time, while simultaneously
!        forming (q transpose)*fvec and storing the first
!        n components in qtf.
!
      do j = 1 , n
         Qtf(j) = zero
         do i = 1 , n
            Fjac(i,j) = zero
         enddo
      enddo
      iflag = 2
      do i = 1 , m
         call fcn(m,n,x,Fvec,Wa3,iflag)
         if ( iflag<0 ) goto 200
         temp = Fvec(i)
         call rwupdt(n,Fjac,Ldfjac,Wa3,Qtf,temp,Wa1,Wa2)
         iflag = iflag + 1
      enddo
      Njev = Njev + 1
!
!        if the jacobian is rank deficient, call qrfac to
!        reorder its columns and update the components of qtf.
!
      sing = .false.
      do j = 1 , n
         if ( Fjac(j,j)==zero ) sing = .true.
         Ipvt(j) = j
         Wa2(j) = enorm(j,Fjac(1,j))
      enddo
      if ( sing ) then
         call qrfac(n,n,Fjac,Ldfjac,.true.,Ipvt,n,Wa1,Wa2,Wa3)
         do j = 1 , n
            if ( Fjac(j,j)/=zero ) then
               sum = zero
               do i = j , n
                  sum = sum + Fjac(i,j)*Qtf(i)
               enddo
               temp = -sum/Fjac(j,j)
               do i = j , n
                  Qtf(i) = Qtf(i) + Fjac(i,j)*temp
               enddo
            endif
            Fjac(j,j) = Wa1(j)
         enddo
      endif
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      if ( iter==1 ) then
         if ( Mode/=2 ) then
            do j = 1 , n
               Diag(j) = Wa2(j)
               if ( Wa2(j)==zero ) Diag(j) = one
            enddo
         endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do j = 1 , n
            Wa3(j) = Diag(j)*x(j)
         enddo
         xnorm = enorm(n,Wa3)
         delta = Factor*xnorm
         if ( delta==zero ) delta = Factor
      endif
!
!        compute the norm of the scaled gradient.
!
      gnorm = zero
      if ( fnorm/=zero ) then
         do j = 1 , n
            l = Ipvt(j)
            if ( Wa2(l)/=zero ) then
               sum = zero
               do i = 1 , j
                  sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
               enddo
               gnorm = dmax1(gnorm,dabs(sum/Wa2(l)))
            endif
         enddo
      endif
!
!        test for convergence of the gradient norm.
!
      if ( gnorm<=Gtol ) Info = 4
      if ( Info==0 ) then
!
!        rescale if necessary.
!
         if ( Mode/=2 ) then
            do j = 1 , n
               Diag(j) = dmax1(Diag(j),Wa2(j))
            enddo
         endif
!
!        beginning of the inner loop.
!
!
!           determine the levenberg-marquardt parameter.
!
 150     call lmpar(n,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1,Wa2,Wa3,  &
                  & Wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
         do j = 1 , n
            Wa1(j) = -Wa1(j)
            Wa2(j) = x(j) + Wa1(j)
            Wa3(j) = Diag(j)*Wa1(j)
         enddo
         pnorm = enorm(n,Wa3)
!
!           on the first iteration, adjust the initial step bound.
!
         if ( iter==1 ) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
         iflag = 1
         call fcn(m,n,Wa2,Wa4,Wa3,iflag)
         Nfev = Nfev + 1
         if ( iflag>=0 ) then
            fnorm1 = enorm(m,Wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if ( p1*fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            do j = 1 , n
               Wa3(j) = zero
               l = Ipvt(j)
               temp = Wa1(l)
               do i = 1 , j
                  Wa3(i) = Wa3(i) + Fjac(i,j)*temp
               enddo
            enddo
            temp1 = enorm(n,Wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2+temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if ( prered/=zero ) ratio = actred/prered
!
!           update the step bound.
!
            if ( ratio<=p25 ) then
               if ( actred>=zero ) temp = p5
               if ( actred<zero ) temp = p5*dirder/(dirder+p5*actred)
               if ( p1*fnorm1>=fnorm .or. temp<p1 ) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
            elseif ( par==zero .or. ratio>=p75 ) then
               delta = pnorm/p5
               par = p5*par
            endif
!
!           test for successful iteration.
!
            if ( ratio>=p0001 ) then
!
!           successful iteration. update x, fvec, and their norms.
!
               do j = 1 , n
                  x(j) = Wa2(j)
                  Wa2(j) = Diag(j)*x(j)
               enddo
               do i = 1 , m
                  Fvec(i) = Wa4(i)
               enddo
               xnorm = enorm(n,Wa2)
               fnorm = fnorm1
               iter = iter + 1
            endif
!
!           tests for convergence.
!
            if ( dabs(actred)<=Ftol .and. prered<=Ftol .and.            &
               & p5*ratio<=one ) Info = 1
            if ( delta<=Xtol*xnorm ) Info = 2
            if ( dabs(actred)<=Ftol .and. prered<=Ftol .and.            &
               & p5*ratio<=one .and. Info==2 ) Info = 3
            if ( Info==0 ) then
!
!           tests for termination and stringent tolerances.
!
               if ( Nfev>=Maxfev ) Info = 5
               if ( dabs(actred)<=epsmch .and. prered<=epsmch .and.     &
                  & p5*ratio<=one ) Info = 6
               if ( delta<=epsmch*xnorm ) Info = 7
               if ( gnorm<=epsmch ) Info = 8
               if ( Info==0 ) then
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
!
!        end of the outer loop.
!
                  if ( ratio>=p0001 ) goto 100
                  goto 150
               endif
            endif
         endif
      endif
!
!     termination, either normal or user imposed.
!
 200  if ( iflag<0 ) Info = iflag
      iflag = 0
      if ( Nprint>0 ) call fcn(m,n,x,Fvec,Wa3,iflag)
!
!     last card of subroutine lmstr.
!
      end

      subroutine lmstr1(fcn,m,n,x,Fvec,Fjac,Ldfjac,Tol,Info,Ipvt,Wa,Lwa)
      implicit none

      integer m , n , Ldfjac , Info , Lwa
      integer Ipvt(n)
      double precision Tol
      double precision x(n) , Fvec(m) , Fjac(Ldfjac,n) , Wa(Lwa)
      external fcn
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
      integer maxfev , mode , nfev , njev , nprint
      double precision factor , ftol , gtol , xtol , zero
      data factor , zero/1.0d2 , 0.0d0/
      Info = 0
!
!     check the input parameters for errors.
!
      if ( n>0 .and. m>=n .and. Ldfjac>=n .and. Tol>=zero .and.         &
         & Lwa>=5*n+m ) then
!
!     call lmstr.
!
         maxfev = 100*(n+1)
         ftol = Tol
         xtol = Tol
         gtol = zero
         mode = 1
         nprint = 0
         call lmstr(fcn,m,n,x,Fvec,Fjac,Ldfjac,ftol,xtol,gtol,maxfev,   &
                  & Wa(1),mode,factor,nprint,Info,nfev,njev,Ipvt,Wa(n+1)&
                  & ,Wa(2*n+1),Wa(3*n+1),Wa(4*n+1),Wa(5*n+1))
         if ( Info==8 ) Info = 4
      endif
!
!     last card of subroutine lmstr1.
!
      end

      subroutine qform(m,n,q,Ldq,Wa)
      implicit none

      integer m , n , Ldq
      double precision q(Ldq,m) , Wa(m)
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
      integer i , j , jm1 , k , l , minmn , np1
      double precision one , sum , temp , zero
      data one , zero/1.0d0 , 0.0d0/
!
!     zero out upper triangle of q in the first min(m,n) columns.
!
      minmn = min0(m,n)
      if ( minmn>=2 ) then
         do j = 2 , minmn
            jm1 = j - 1
            do i = 1 , jm1
               q(i,j) = zero
            enddo
         enddo
      endif
!
!     initialize remaining columns to those of the identity matrix.
!
      np1 = n + 1
      if ( m>=np1 ) then
         do j = np1 , m
            do i = 1 , m
               q(i,j) = zero
            enddo
            q(j,j) = one
         enddo
      endif
!
!     accumulate q from its factored form.
!
      do l = 1 , minmn
         k = minmn - l + 1
         do i = k , m
            Wa(i) = q(i,k)
            q(i,k) = zero
         enddo
         q(k,k) = one
         if ( Wa(k)/=zero ) then
            do j = k , m
               sum = zero
               do i = k , m
                  sum = sum + q(i,j)*Wa(i)
               enddo
               temp = sum/Wa(k)
               do i = k , m
                  q(i,j) = q(i,j) - temp*Wa(i)
               enddo
            enddo
         endif
      enddo
!
!     last card of subroutine qform.
!
      end

      subroutine qrfac(m,n,a,Lda,Pivot,Ipvt,Lipvt,Rdiag,Acnorm,Wa)
      implicit none

      integer m , n , Lda , Lipvt
      integer Ipvt(Lipvt)
      logical Pivot
      double precision a(Lda,n) , Rdiag(n) , Acnorm(n) , Wa(n)
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
      integer i , j , jp1 , k , kmax , minmn
      double precision ajnorm , epsmch , one , p05 , sum , temp , zero
      double precision dpmpar , enorm
      data one , p05 , zero/1.0d0 , 5.0d-2 , 0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
!     compute the initial column norms and initialize several arrays.
!
      do j = 1 , n
         Acnorm(j) = enorm(m,a(1,j))
         Rdiag(j) = Acnorm(j)
         Wa(j) = Rdiag(j)
         if ( Pivot ) Ipvt(j) = j
      enddo
!
!     reduce a to r with householder transformations.
!
      minmn = min0(m,n)
      do j = 1 , minmn
         if ( Pivot ) then
!
!        bring the column of largest norm into the pivot position.
!
            kmax = j
            do k = j , n
               if ( Rdiag(k)>Rdiag(kmax) ) kmax = k
            enddo
            if ( kmax/=j ) then
               do i = 1 , m
                  temp = a(i,j)
                  a(i,j) = a(i,kmax)
                  a(i,kmax) = temp
               enddo
               Rdiag(kmax) = Rdiag(j)
               Wa(kmax) = Wa(j)
               k = Ipvt(j)
               Ipvt(j) = Ipvt(kmax)
               Ipvt(kmax) = k
            endif
         endif
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = enorm(m-j+1,a(j,j))
         if ( ajnorm/=zero ) then
            if ( a(j,j)<zero ) ajnorm = -ajnorm
            do i = j , m
               a(i,j) = a(i,j)/ajnorm
            enddo
            a(j,j) = a(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
            jp1 = j + 1
            if ( n>=jp1 ) then
               do k = jp1 , n
                  sum = zero
                  do i = j , m
                     sum = sum + a(i,j)*a(i,k)
                  enddo
                  temp = sum/a(j,j)
                  do i = j , m
                     a(i,k) = a(i,k) - temp*a(i,j)
                  enddo
                  if ( .not.(.not.Pivot .or. Rdiag(k)==zero) ) then
                     temp = a(j,k)/Rdiag(k)
                     Rdiag(k) = Rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
                     if ( p05*(Rdiag(k)/Wa(k))**2<=epsmch ) then
                        Rdiag(k) = enorm(m-j,a(jp1,k))
                        Wa(k) = Rdiag(k)
                     endif
                  endif
               enddo
            endif
         endif
         Rdiag(j) = -ajnorm
      enddo
!
!     last card of subroutine qrfac.
!
      end

      subroutine qrsolv(n,r,Ldr,Ipvt,Diag,Qtb,x,Sdiag,Wa)
      implicit none

      integer n , Ldr
      integer Ipvt(n)
      double precision r(Ldr,n) , Diag(n) , Qtb(n) , x(n) , Sdiag(n) ,  &
                     & Wa(n)
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
      integer i , j , jp1 , k , kp1 , l , nsing
      double precision cos , cotan , p5 , p25 , qtbpj , sin , sum ,     &
                     & tan , temp , zero
      data p5 , p25 , zero/5.0d-1 , 2.5d-1 , 0.0d0/
!
!     copy r and (q transpose)*b to preserve input and initialize s.
!     in particular, save the diagonal elements of r in x.
!
      do j = 1 , n
         do i = j , n
            r(i,j) = r(j,i)
         enddo
         x(j) = r(j,j)
         Wa(j) = Qtb(j)
      enddo
!
!     eliminate the diagonal matrix d using a givens rotation.
!
      do j = 1 , n
!
!        prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
!
         l = Ipvt(j)
         if ( Diag(l)/=zero ) then
            do k = j , n
               Sdiag(k) = zero
            enddo
            Sdiag(j) = Diag(l)
!
!        the transformations to eliminate the row of d
!        modify only a single element of (q transpose)*b
!        beyond the first n, which is initially zero.
!
            qtbpj = zero
            do k = j , n
!
!           determine a givens rotation which eliminates the
!           appropriate element in the current row of d.
!
               if ( Sdiag(k)/=zero ) then
                  if ( dabs(r(k,k))>=dabs(Sdiag(k)) ) then
                     tan = Sdiag(k)/r(k,k)
                     cos = p5/dsqrt(p25+p25*tan**2)
                     sin = cos*tan
                  else
                     cotan = r(k,k)/Sdiag(k)
                     sin = p5/dsqrt(p25+p25*cotan**2)
                     cos = sin*cotan
                  endif
!
!           compute the modified diagonal element of r and
!           the modified element of ((q transpose)*b,0).
!
                  r(k,k) = cos*r(k,k) + sin*Sdiag(k)
                  temp = cos*Wa(k) + sin*qtbpj
                  qtbpj = -sin*Wa(k) + cos*qtbpj
                  Wa(k) = temp
!
!           accumulate the tranformation in the row of s.
!
                  kp1 = k + 1
                  if ( n>=kp1 ) then
                     do i = kp1 , n
                        temp = cos*r(i,k) + sin*Sdiag(i)
                        Sdiag(i) = -sin*r(i,k) + cos*Sdiag(i)
                        r(i,k) = temp
                     enddo
                  endif
               endif
            enddo
         endif
!
!        store the diagonal element of s and restore
!        the corresponding diagonal element of r.
!
         Sdiag(j) = r(j,j)
         r(j,j) = x(j)
      enddo
!
!     solve the triangular system for z. if the system is
!     singular, then obtain a least squares solution.
!
      nsing = n
      do j = 1 , n
         if ( Sdiag(j)==zero .and. nsing==n ) nsing = j - 1
         if ( nsing<n ) Wa(j) = zero
      enddo
      if ( nsing>=1 ) then
         do k = 1 , nsing
            j = nsing - k + 1
            sum = zero
            jp1 = j + 1
            if ( nsing>=jp1 ) then
               do i = jp1 , nsing
                  sum = sum + r(i,j)*Wa(i)
               enddo
            endif
            Wa(j) = (Wa(j)-sum)/Sdiag(j)
         enddo
      endif
!
!     permute the components of z back to components of x.
!
      do j = 1 , n
         l = Ipvt(j)
         x(l) = Wa(j)
      enddo
!
!     last card of subroutine qrsolv.
!
      end

      subroutine r1mpyq(m,n,a,Lda,v,w)
      implicit none

      integer m , n , Lda
      double precision a(Lda,n) , v(n) , w(n)
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
      integer i , j , nmj , nm1
      double precision cos , one , sin , temp
      data one/1.0d0/
!
!     apply the first set of givens rotations to a.
!
      nm1 = n - 1
      if ( nm1>=1 ) then
         do nmj = 1 , nm1
            j = n - nmj
            if ( dabs(v(j))>one ) cos = one/v(j)
            if ( dabs(v(j))>one ) sin = dsqrt(one-cos**2)
            if ( dabs(v(j))<=one ) sin = v(j)
            if ( dabs(v(j))<=one ) cos = dsqrt(one-sin**2)
            do i = 1 , m
               temp = cos*a(i,j) - sin*a(i,n)
               a(i,n) = sin*a(i,j) + cos*a(i,n)
               a(i,j) = temp
            enddo
         enddo
!
!     apply the second set of givens rotations to a.
!
         do j = 1 , nm1
            if ( dabs(w(j))>one ) cos = one/w(j)
            if ( dabs(w(j))>one ) sin = dsqrt(one-cos**2)
            if ( dabs(w(j))<=one ) sin = w(j)
            if ( dabs(w(j))<=one ) cos = dsqrt(one-sin**2)
            do i = 1 , m
               temp = cos*a(i,j) + sin*a(i,n)
               a(i,n) = -sin*a(i,j) + cos*a(i,n)
               a(i,j) = temp
            enddo
         enddo
      endif
!
!     last card of subroutine r1mpyq.
!
      end

      subroutine r1updt(m,n,s,Ls,u,v,w,Sing)
      implicit none

      integer m , n , Ls
      logical Sing
      double precision s(Ls) , u(m) , v(n) , w(m)
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
      integer i , j , jj , l , nmj , nm1
      double precision cos , cotan , giant , one , p5 , p25 , sin ,     &
                     & tan , tau , temp , zero
      double precision dpmpar
      data one , p5 , p25 , zero/1.0d0 , 5.0d-1 , 2.5d-1 , 0.0d0/
!
!     giant is the largest magnitude.
!
      giant = dpmpar(3)
!
!     initialize the diagonal element pointer.
!
      jj = (n*(2*m-n+1))/2 - (m-n)
!
!     move the nontrivial part of the last column of s into w.
!
      l = jj
      do i = n , m
         w(i) = s(l)
         l = l + 1
      enddo
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
      nm1 = n - 1
      if ( nm1>=1 ) then
         do nmj = 1 , nm1
            j = n - nmj
            jj = jj - (m-j+1)
            w(j) = zero
            if ( v(j)/=zero ) then
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
               if ( dabs(v(n))>=dabs(v(j)) ) then
                  tan = v(j)/v(n)
                  cos = p5/dsqrt(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               else
                  cotan = v(n)/v(j)
                  sin = p5/dsqrt(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  if ( dabs(cos)*giant>one ) tau = one/cos
               endif
!
!        apply the transformation to v and store the information
!        necessary to recover the givens rotation.
!
               v(n) = sin*v(j) + cos*v(n)
               v(j) = tau
!
!        apply the transformation to s and extend the spike in w.
!
               l = jj
               do i = j , m
                  temp = cos*s(l) - sin*w(i)
                  w(i) = sin*s(l) + cos*w(i)
                  s(l) = temp
                  l = l + 1
               enddo
            endif
         enddo
      endif
!
!     add the spike from the rank 1 update to w.
!
      do i = 1 , m
         w(i) = w(i) + v(n)*u(i)
      enddo
!
!     eliminate the spike.
!
      Sing = .false.
      if ( nm1>=1 ) then
         do j = 1 , nm1
            if ( w(j)/=zero ) then
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
               if ( dabs(s(jj))>=dabs(w(j)) ) then
                  tan = w(j)/s(jj)
                  cos = p5/dsqrt(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               else
                  cotan = s(jj)/w(j)
                  sin = p5/dsqrt(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  if ( dabs(cos)*giant>one ) tau = one/cos
               endif
!
!        apply the transformation to s and reduce the spike in w.
!
               l = jj
               do i = j , m
                  temp = cos*s(l) + sin*w(i)
                  w(i) = -sin*s(l) + cos*w(i)
                  s(l) = temp
                  l = l + 1
               enddo
!
!        store the information necessary to recover the
!        givens rotation.
!
               w(j) = tau
            endif
!
!        test for zero diagonal elements in the output s.
!
            if ( s(jj)==zero ) Sing = .true.
            jj = jj + (m-j+1)
         enddo
      endif
!
!     move w back into the last column of the output s.
!
      l = jj
      do i = n , m
         s(l) = w(i)
         l = l + 1
      enddo
      if ( s(jj)==zero ) Sing = .true.
!
!     last card of subroutine r1updt.
!
      end

      subroutine rwupdt(n,r,Ldr,w,b,Alpha,Cos,Sin)
      implicit none

      integer n , Ldr
      double precision Alpha
      double precision r(Ldr,n) , w(n) , b(n) , Cos(n) , Sin(n)
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
      integer i , j , jm1
      double precision cotan , one , p5 , p25 , rowj , tan , temp , zero
      data one , p5 , p25 , zero/1.0d0 , 5.0d-1 , 2.5d-1 , 0.0d0/
!
      do j = 1 , n
         rowj = w(j)
         jm1 = j - 1
!
!        apply the previous transformations to
!        r(i,j), i=1,2,...,j-1, and to w(j).
!
         if ( jm1>=1 ) then
            do i = 1 , jm1
               temp = Cos(i)*r(i,j) + Sin(i)*rowj
               rowj = -Sin(i)*r(i,j) + Cos(i)*rowj
               r(i,j) = temp
            enddo
         endif
!
!        determine a givens rotation which eliminates w(j).
!
         Cos(j) = one
         Sin(j) = zero
         if ( rowj/=zero ) then
            if ( dabs(r(j,j))>=dabs(rowj) ) then
               tan = rowj/r(j,j)
               Cos(j) = p5/dsqrt(p25+p25*tan**2)
               Sin(j) = Cos(j)*tan
            else
               cotan = r(j,j)/rowj
               Sin(j) = p5/dsqrt(p25+p25*cotan**2)
               Cos(j) = Sin(j)*cotan
            endif
!
!        apply the current transformation to r(j,j), b(j), and alpha.
!
            r(j,j) = Cos(j)*r(j,j) + Sin(j)*rowj
            temp = Cos(j)*b(j) + Sin(j)*Alpha
            Alpha = -Sin(j)*b(j) + Cos(j)*Alpha
            b(j) = temp
         endif
      enddo
!
!     last card of subroutine rwupdt.
!
      end
