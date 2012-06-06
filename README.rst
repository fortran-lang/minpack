Minpack
=======

Information
-----------

This repository contains the original double precision Minpack from netlib.org,
together with CMake makefiles and examples.

About Minpack
-------------

Minpack includes software for solving nonlinear equations and
nonlinear least squares problems.  Five algorithmic paths each include
a core subroutine and an easy-to-use driver.  The algorithms proceed
either from an analytic specification of the Jacobian matrix or
directly from the problem functions.  The paths include facilities for
systems of equations with a banded Jacobian matrix, for least squares
problems with a large amount of data, and for checking the consistency
of the Jacobian matrix with the functions.

Jorge Moré, Burt Garbow, and Ken Hillstrom at Argonne National Laboratory.

Documentation
-------------

Minpack contains 4 subroutines for solution of systems of nonlinear equations:

* ``hybrd``, ``hybrd1``: Jacobian matrix is calculated by a forward difference
  approximation
* ``hybrj``, ``hybrj1``: Jacobian matrix is provided by the user

and 6 subroutines for nonlinear least squares problems:

* ``lmdif``, ``lmdif1``: Jacobian matrix is calculated by a forward difference
  approximation
* ``lmder``, ``lmder1``: Jacobian matrix is provided by the user
* ``lmstr``, ``lmstr1``: Jacobian matrix is provided by the user, one row per
  call (uses less memory)

The routines without ``1`` in the name expose all parameters to the user (`core
subroutines`), routines with ``1`` only expose the essential parameters and set
default values for the rest (`easy-to-use driver`). Finally:

* ``chkder``: checks the consistency of the Jacobian matrix with the functions

More general documentation is given in
the 1980 Argonne technical report written by the authors of Minpack,
`Chapters 1-3 <http://www.mcs.anl.gov/~more/ANL8074a.pdf>`_.
The `Chapter 4 <http://www.mcs.anl.gov/~more/ANL8074b.pdf>`_ (also available in
the file ``ex/file06``) contains detailed documentation for all these routines
together with an example of usage.  Ready to use examples of usage are in the
``examples`` directory.

Other files in the ``ex`` directory are original examples of usage of various
routines (single and double precision), but are not compiled by default.
