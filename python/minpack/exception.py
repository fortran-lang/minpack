# -*- coding: utf-8 -*-
"""
Possible exceptions for the `minpack` package mapping the info codes produced
by the library to exceptions.
"""

from typing import Optional, Type


class MinpackError(Exception):
    """
    Exception raised when Minpack returns an error.
    """

    pass


class MinpackInputError(MinpackError):
    """
    Exception raised when Minpack input is invalid.
    """

    pass


class MinpackMaxIterations(MinpackError):
    """
    The maximum number of calls to the objective function is reached.
    """

    pass


class MinpackFunctionTolerance(MinpackError):
    """
    `ftol` is too small. No further reduction in the sum of squares is possible.
    """

    pass


class MinpackSolutionTolerance(MinpackError):
    """
    `xtol` is too small. No further improvement in the approximate
    solution x is possible.
    """

    pass


class MinpackJacobianTolerance(MinpackError):
    """
    `gtol` is too small. `fvec` is orthogonal to the columns of
    the Jacobian to machine precision.
    """

    pass


class MinpackSlowProgress(MinpackError):
    """
    Iteration is not making good progress, as measured by the improvement
    from the last ten iterations.
    """

    pass


class MinpackSlowProgressJacobian(MinpackError):
    """
    Iteration is not making good progress, as measured by the improvement
    from the last five jacobian evaluations.
    """


def info_hy(info: int) -> Optional[Type[MinpackError]]:
    """
    Get possible errors for `hybrd` and `hybrj` drivers.
    """

    return {
        0: MinpackInputError,
        2: MinpackMaxIterations,
        3: MinpackFunctionTolerance,
        4: MinpackSlowProgressJacobian,
        5: MinpackSlowProgress,
    }.get(info)


def info_lm(info: int) -> Optional[Type[MinpackError]]:
    """
    Get possible errors for `lmdif`, `lmder`, and `lmstr` drivers.
    """

    return {
        0: MinpackInputError,
        5: MinpackMaxIterations,
        6: MinpackFunctionTolerance,
        7: MinpackSolutionTolerance,
        8: MinpackJacobianTolerance,
    }.get(info)
