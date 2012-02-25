      double precision function dpmpar(i)
      integer i
c     **********
c
c     Function dpmpar
c
c     This function provides double precision machine parameters
c     when the appropriate set of data statements is activated (by
c     removing the c from column 1) and all other data statements are
c     rendered inactive. Most of the parameter values were obtained
c     from the corresponding Bell Laboratories Port Library function.
c
c     The function statement is
c
c       double precision function dpmpar(i)
c
c     where
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. If the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     Argonne National Laboratory. MINPACK Project. November 1996.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
c
c     **********

      double precision dmach(3)
      dmach(1) = epsilon(1d0)
      dmach(2) = tiny(1d0)
      dmach(3) = huge(1d0)
c
      dpmpar = dmach(i)
      return
c
c     Last card of function dpmpar.
c
      end
