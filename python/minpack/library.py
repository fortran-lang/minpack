# -*- coding: utf-8 -*-
"""
Low-level Python bindings to the Minpack library.

This module forwards the CFFI generated bindings to the Minpack library and provides
a Pythonic interface to the C API.
"""

import numpy as np
import functools
import math
from typing import Optional

from .typing import (
    CallableHybrd,
    CallableHybrj,
    CallableLmder,
    CallableLmdif,
    CallableLmstr,
)
from ._libminpack import ffi, lib
from .exception import info_hy, info_lm


class UserData:
    """
    Carry Python callable object through callback and propagate exceptions without
    disturbing foreign runtime.
    """

    def __init__(self, fcn, **kwargs):
        self.fcn = fcn
        self.exception = None
        self.kwargs = kwargs


@ffi.def_extern()
def func(n, x, fvec, iflag, data) -> None:
    """
    Entry point for callback from minpack_hybrd and minpack_hybrd1 library functions.
    Restores type information for NDArray objects and calls the user-provided callback.
    """
    if iflag[0] <= 0:
        return
    handle: UserData = ffi.from_handle(data)
    try:
        handle.fcn(
            np.frombuffer(ffi.buffer(x, n * real.itemsize), dtype=real),
            np.frombuffer(ffi.buffer(fvec, n * real.itemsize), dtype=real),
            **handle.kwargs,
        )
    except BaseException as e:
        iflag[0] = -1
        handle.exception = e


@ffi.def_extern()
def fcn_hybrj(n, x, fvec, fjac, ldfjac, iflag, data) -> None:
    """
    Entry point for callback from minpack_hybrj and minpack_hybrj1 library functions.
    Restores type information for NDArray objects and calls the user-provided callback.
    """
    if iflag[0] <= 0:
        return
    handle: UserData = ffi.from_handle(data)
    try:
        fjac = np.frombuffer(ffi.buffer(fjac, ldfjac * n * real.itemsize), dtype=real)
        handle.fcn(
            np.frombuffer(ffi.buffer(x, n * real.itemsize), dtype=real),
            np.frombuffer(ffi.buffer(fvec, n * real.itemsize), dtype=real),
            np.reshape(fjac, (n, ldfjac)),
            iflag[0] == 2,
            **handle.kwargs,
        )
    except BaseException as e:
        iflag[0] = -1
        handle.exception = e


@ffi.def_extern()
def fcn_lmder(m, n, x, fvec, fjac, ldfjac, iflag, data) -> None:
    """
    Entry point for callback from minpack_lmder and minpack_lmder1 library functions.
    Restores type information for NDArray objects and calls the user-provided callback.
    """
    if iflag[0] <= 0:
        return
    handle: UserData = ffi.from_handle(data)
    try:
        fjac = np.frombuffer(ffi.buffer(fjac, ldfjac * n * real.itemsize), dtype=real)
        handle.fcn(
            np.frombuffer(ffi.buffer(x, n * real.itemsize), dtype=real),
            np.frombuffer(ffi.buffer(fvec, m * real.itemsize), dtype=real),
            np.reshape(fjac, (n, ldfjac)),
            iflag[0] == 2,
            **handle.kwargs,
        )
    except BaseException as e:
        iflag[0] = -1
        handle.exception = e


@ffi.def_extern()
def func2(m, n, x, fvec, iflag, data) -> None:
    """
    Entry point for callback from minpack_lmdif and minpack_lmdif1 library functions.
    Restores type information for NDArray objects and calls the user-provided callback.
    """
    if iflag[0] <= 0:
        return
    handle: UserData = ffi.from_handle(data)
    try:
        handle.fcn(
            np.frombuffer(ffi.buffer(x, n * real.itemsize), dtype=real),
            np.frombuffer(ffi.buffer(fvec, m * real.itemsize), dtype=real),
            **handle.kwargs,
        )
    except BaseException as e:
        iflag[0] = -1
        handle.exception = e


@ffi.def_extern()
def fcn_lmstr(m, n, x, fvec, fjrow, iflag, data) -> None:
    """
    Entry point for callback from minpack_lmstr and minpack_lmstr1 library functions.
    Restores type information for NDArray objects and calls the user-provided callback.
    """
    if iflag[0] <= 0:
        return
    handle: UserData = ffi.from_handle(data)
    try:
        handle.fcn(
            np.frombuffer(ffi.buffer(x, n * real.itemsize), dtype=real),
            np.frombuffer(ffi.buffer(fvec, m * real.itemsize), dtype=real),
            np.frombuffer(ffi.buffer(fjrow, n * real.itemsize), dtype=real),
            iflag[0] - 2 if iflag[0] > 1 else None,
            **handle.kwargs,
        )
    except BaseException as e:
        iflag[0] = -1
        handle.exception = e


def extern_python(dec):
    """
    Meta-decorator to attach a CFFI extern "Python" callback to a decorator
    handling the Python side of the callback.
    """

    def layer(*args, **kwargs):
        def wrapper(func):
            return dec(func, *args, **kwargs)

        return wrapper

    return layer


@extern_python
def cffi_callback(func, callback):
    """
    Attach Python callback to a library function with extern "Python" callback.

    This decorator wraps the user-provided Python callback in a `UserData` object
    to carry it through the foreign runtime. It also propagates exceptions from
    the user-provided Python callback back through the foreign runtime and re-raises.
    """

    @functools.wraps(func)
    def entry_point(fcn, *args, **kwargs):
        data = UserData(fcn, **kwargs)
        handle = ffi.new_handle(data)
        func(callback, *args, handle)
        if data.exception is not None:
            raise data.exception

    return entry_point


real = np.dtype("f8")
minpack_hybrd1 = cffi_callback(lib.func)(lib.minpack_hybrd1)
minpack_hybrd = cffi_callback(lib.func)(lib.minpack_hybrd)
minpack_hybrj1 = cffi_callback(lib.fcn_hybrj)(lib.minpack_hybrj1)
minpack_hybrj = cffi_callback(lib.fcn_hybrj)(lib.minpack_hybrj)
minpack_lmder1 = cffi_callback(lib.fcn_lmder)(lib.minpack_lmder1)
minpack_lmder = cffi_callback(lib.fcn_lmder)(lib.minpack_lmder)
minpack_lmdif1 = cffi_callback(lib.func2)(lib.minpack_lmdif1)
minpack_lmdif = cffi_callback(lib.func2)(lib.minpack_lmdif)
minpack_lmstr1 = cffi_callback(lib.fcn_lmstr)(lib.minpack_lmstr1)
minpack_lmstr = cffi_callback(lib.fcn_lmstr)(lib.minpack_lmstr)
minpack_chkder = lib.minpack_chkder


def hybrd1(
    fcn: CallableHybrd,
    x: np.ndarray,
    fvec: np.ndarray,
    tol: float = math.sqrt(np.finfo(real).eps),
    **kwargs,
) -> int:
    """
    Find a zero of a system of n nonlinear functions in n variables
    by a modification of the Powell hybrid method.
    This is done by using the more general nonlinear equation solver `hybrd`.
    The user must provide a subroutine which calculates the functions.
    The Jacobian is then calculated by a forward-difference approximation.

    Parameters
    ----------
    func : callable ``f(x, fvec)``
        A function that takes at least one (possibly vector) argument,
        and returns a value of the same length.
    x : ndarray
        The starting estimate for the roots of ``func(x) = 0``.
    fvec: ndarray
        Function evaluated at the output
    tol : float
        The calculation will terminate if the relative error between two
        consecutive iterates is at most `tol`.

    Returns
    -------
    info : int
        Set to 1 if a solution was found.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSlowProgressJacobian
        When the Jacobian is not changing.
    MinpackSlowProgress
        When the function is not changing.

    Examples
    --------
    >>> import numpy as np
    >>> from minpack import hybrd
    >>>
    >>> def fcn(x, fvec) -> None:
    ...     fvec[0] = x[0] * np.cos(x[1]) - 4
    ...     fvec[1] = x[1] * x[0] - x[1] - 5
    ...
    >>> x = np.array(2 * [1.0])
    >>> fvec = np.zeros(2, dtype=np.float64)
    >>> hybrd1(fcn, x, fvec)
    1
    >>> x
    array([6.50409711, 0.90841421])
    >>> np.isclose(fvec, [0.0, 0.0])  # fvec should be almost 0.0.
    array([ True,  True])
    """
    n = x.size
    lwa = n * (3 * n + 13) // 2
    wa = np.zeros(lwa, dtype=real)
    info = ffi.new("int *")

    minpack_hybrd1(
        fcn,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        tol,
        info,
        ffi.cast("double*", wa.ctypes.data),
        lwa,
        **kwargs,
    )
    ex = info_hy(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def hybrd(
    fcn: CallableHybrd,
    x: np.ndarray,
    fvec: np.ndarray,
    xtol: float = math.sqrt(np.finfo(real).eps),
    *,
    maxfev: Optional[int] = None,
    ml: Optional[int] = None,
    mu: Optional[int] = None,
    epsfcn: float = 0.0,
    diag: Optional[np.ndarray] = None,
    mode: int = 2,
    factor: float = 100.0,
    nprint: int = 0,
    fjac: Optional[np.ndarray] = None,
    r: Optional[np.ndarray] = None,
    qtf: Optional[np.ndarray] = None,
    **kwargs,
) -> int:
    """
    Find a zero of a system of n nonlinear functions in n variables
    by a modification of the Powell hybrid method.
    The user must provide a subroutine which calculates the functions.
    The Jacobian is then calculated by a forward-difference approximation.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSlowProgressJacobian
        When the Jacobian is not changing.
    MinpackSlowProgress
        When the function is not changing.
    """
    n = x.size
    info = ffi.new("int *")
    nfev = ffi.new("int *")
    if maxfev is None:
        maxfev = 200 * (n + 1)
    if ml is None:
        ml = n - 1
    if mu is None:
        mu = n - 1
    if diag is None:
        diag = np.ones(n, dtype=real)
    if fjac is None:
        fjac = np.zeros((n, n), dtype=real)
    if r is None:
        r = np.zeros(n * (n + 1) // 2, dtype=real)
    if qtf is None:
        qtf = np.zeros(n, dtype=real)
    wa1 = np.zeros(n, dtype=real)
    wa2 = np.zeros(n, dtype=real)
    wa3 = np.zeros(n, dtype=real)
    wa4 = np.zeros(n, dtype=real)

    minpack_hybrd(
        fcn,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        xtol,
        maxfev,
        ml,
        mu,
        epsfcn,
        ffi.cast("double*", diag.ctypes.data),
        mode,
        factor,
        nprint,
        info,
        nfev,
        ffi.cast("double*", fjac.ctypes.data),
        n,
        ffi.cast("double*", r.ctypes.data),
        r.size,
        ffi.cast("double*", qtf.ctypes.data),
        ffi.cast("double*", wa1.ctypes.data),
        ffi.cast("double*", wa2.ctypes.data),
        ffi.cast("double*", wa3.ctypes.data),
        ffi.cast("double*", wa4.ctypes.data),
        **kwargs,
    )
    ex = info_hy(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def hybrj1(
    fcn: CallableHybrj,
    x: np.ndarray,
    fvec: np.ndarray,
    fjac: np.ndarray,
    tol: float = math.sqrt(np.finfo(real).eps),
    **kwargs,
) -> int:
    """
    Find a zero of a system of n nonlinear functions in n variables
    by a modification of the Powell hybrid method.
    This is done by using the more general nonlinear equation solver `hybrj`.
    The user must provide a subroutine which calculates the functions
    and the Jacobian.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSlowProgressJacobian
        When the Jacobian is not changing.
    MinpackSlowProgress
        When the function is not changing.
    """
    n = x.size
    lwa = (n * (n + 13)) // 2
    wa = np.zeros(lwa, dtype=real)
    info = ffi.new("int *")

    minpack_hybrj1(
        fcn,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        n,
        tol,
        info,
        ffi.cast("double*", wa.ctypes.data),
        lwa,
        **kwargs,
    )
    ex = info_hy(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def hybrj(
    fcn: CallableHybrj,
    x: np.ndarray,
    fvec: np.ndarray,
    fjac: np.ndarray,
    xtol: float = math.sqrt(np.finfo(real).eps),
    *,
    maxfev: Optional[int] = None,
    diag: Optional[np.ndarray] = None,
    mode: int = 2,
    factor: float = 100.0,
    nprint: int = 0,
    r: Optional[np.ndarray] = None,
    qtf: Optional[np.ndarray] = None,
    **kwargs,
) -> int:
    """
    Find a zero of a system of n nonlinear functions in n variables
    by a modification of the Powell hybrid method.
    The user must provide a subroutine which calculates the functions
    and the Jacobian.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSlowProgressJacobian
        When the Jacobian is not changing.
    MinpackSlowProgress
        When the function is not changing.
    """
    n = x.size
    info = ffi.new("int *")
    nfev = ffi.new("int *")
    njev = ffi.new("int *")
    if maxfev is None:
        maxfev = 200 * (n + 1)
    if diag is None:
        diag = np.ones(n, dtype=real)
    if fjac is None:
        fjac = np.zeros((n, n), dtype=real)
    if r is None:
        r = np.zeros(n * (n + 1) // 2, dtype=real)
    if qtf is None:
        qtf = np.zeros(n, dtype=real)
    wa1 = np.zeros(n, dtype=real)
    wa2 = np.zeros(n, dtype=real)
    wa3 = np.zeros(n, dtype=real)
    wa4 = np.zeros(n, dtype=real)

    minpack_hybrj(
        fcn,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        n,
        xtol,
        maxfev,
        ffi.cast("double*", diag.ctypes.data),
        mode,
        factor,
        nprint,
        info,
        nfev,
        njev,
        ffi.cast("double*", r.ctypes.data),
        r.size,
        ffi.cast("double*", qtf.ctypes.data),
        ffi.cast("double*", wa1.ctypes.data),
        ffi.cast("double*", wa2.ctypes.data),
        ffi.cast("double*", wa3.ctypes.data),
        ffi.cast("double*", wa4.ctypes.data),
        **kwargs,
    )
    ex = info_hy(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def lmder1(
    fcn: CallableLmder,
    x: np.ndarray,
    fvec: np.ndarray,
    fjac: np.ndarray,
    tol: float = math.sqrt(np.finfo(real).eps),
    **kwargs,
) -> int:
    """
    Minimize the sum of the squares of m nonlinear functions in n variables
    by a modification of the Levenberg-Marquardt algorithm.
    This is done by using the more general least-squares solver `lmder`.
    The user must provide a subroutine which calculates the functions and the Jacobian.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSolutionTolerance
        When no further improvement in the approximate solution is possible.
    MinpackJacobianTolerance
        The solution is orthogonal to the jacobian.
    """
    n = x.size
    m = fvec.size
    lwa = 5 * n + m
    wa = np.zeros(lwa, dtype=real)
    ipvt = np.zeros(n, dtype=np.int32)
    info = ffi.new("int *")

    minpack_lmder1(
        fcn,
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        m,
        tol,
        info,
        ffi.cast("int*", ipvt.ctypes.data),
        ffi.cast("double*", wa.ctypes.data),
        lwa,
        **kwargs,
    )
    ex = info_lm(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def lmder(
    fcn: CallableLmder,
    x: np.ndarray,
    fvec: np.ndarray,
    fjac: np.ndarray,
    ftol: float = math.sqrt(np.finfo(real).eps),
    xtol: float = math.sqrt(np.finfo(real).eps),
    *,
    gtol: float = 0.0,
    maxfev: Optional[int] = None,
    diag: Optional[np.ndarray] = None,
    mode: int = 1,
    factor=100.0,
    nprint=0,
    ipvt: Optional[np.ndarray] = None,
    qtf: Optional[np.ndarray] = None,
    **kwargs,
) -> int:
    """
    Minimize the sum of the squares of m nonlinear functions in n variables
    by a modification of the Levenberg-Marquardt algorithm.
    The user must provide a subroutine which calculates the functions and the Jacobian.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSolutionTolerance
        When no further improvement in the approximate solution is possible.
    MinpackJacobianTolerance
        The solution is orthogonal to the jacobian.
    """
    n = x.size
    m = fvec.size
    info = ffi.new("int *")
    nfev = ffi.new("int *")
    njev = ffi.new("int *")
    if diag is None:
        diag = np.ones(n, dtype=real)
    if maxfev is None:
        maxfev = 100 * (n + 1)
    if ipvt is None:
        ipvt = np.zeros(n, dtype=np.int32)
    if qtf is None:
        qtf = np.zeros(n, dtype=real)
    wa1 = np.zeros(n, dtype=real)
    wa2 = np.zeros(n, dtype=real)
    wa3 = np.zeros(n, dtype=real)
    wa4 = np.zeros(m, dtype=real)

    minpack_lmder(
        fcn,
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        m,
        ftol,
        xtol,
        gtol,
        maxfev,
        ffi.cast("double*", diag.ctypes.data),
        mode,
        factor,
        nprint,
        info,
        nfev,
        njev,
        ffi.cast("int*", ipvt.ctypes.data),
        ffi.cast("double*", qtf.ctypes.data),
        ffi.cast("double*", wa1.ctypes.data),
        ffi.cast("double*", wa2.ctypes.data),
        ffi.cast("double*", wa3.ctypes.data),
        ffi.cast("double*", wa4.ctypes.data),
        **kwargs,
    )
    ex = info_lm(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def lmdif1(
    fcn: CallableLmdif,
    x: np.ndarray,
    fvec: np.ndarray,
    tol: float = math.sqrt(np.finfo(real).eps),
    **kwargs,
) -> int:
    """
    Minimize the sum of the squares of m nonlinear functions in n variables
    by a modification of the Levenberg-Marquardt algorithm.
    This is done by using the more general least-squares solver `lmdif`.
    The user must provide a subroutine which calculates the functions.
    The jacobian is then calculated by a forward-difference approximation.

    Parameters
    ----------
    func : callable ``f(x, fvec)``
        Should take at least one (possibly length n vector) argument and
        compute m floating point numbers in fvec. It must not return NaNs or
        fitting might fail. m must be greater than or equal to n.
    x : ndarray
        The starting estimate for the minimization.
    fvec : ndarray
        The function evaluated at the output.
    tol : float, optional
        Relative error desired in the sum of squares and the approximate solution.

    Returns
    -------
    info : int
        An integer flag. If it is equal to 1, 2, 3 or 4, the solution was
        found. Otherwise, the solution was not found.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSolutionTolerance
        When no further improvement in the approximate solution is possible.
    MinpackJacobianTolerance
        The solution is orthogonal to the jacobian.

    Example
    -------
    >>> import numpy as np
    >>> from minpack import lmdif1
    >>>
    >>> def func(x, fvec):
    ...     fvec[:] = 2*(x-3)**2+1
    ...
    >>> x = np.array(0.0)
    >>> fvec = np.zeros(1, dtype=np.float64)
    >>> lmdif1(func, x, fvec)
    1
    >>> x
    array(2.99999999)
    """
    n = x.size
    m = fvec.size
    lwa = m * n + 5 * n + m
    wa = np.zeros(lwa, dtype=real)
    ipvt = np.zeros(n, dtype=np.int32)
    info = ffi.new("int *")

    minpack_lmdif1(
        fcn,
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        tol,
        info,
        ffi.cast("int*", ipvt.ctypes.data),
        ffi.cast("double*", wa.ctypes.data),
        lwa,
        **kwargs,
    )
    ex = info_lm(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def lmdif(
    fcn: CallableLmdif,
    x: np.ndarray,
    fvec: np.ndarray,
    ftol: float = math.sqrt(np.finfo(real).eps),
    xtol: float = math.sqrt(np.finfo(real).eps),
    *,
    gtol: float = 0.0,
    maxfev: Optional[int] = None,
    epsfcn: float = 0.0,
    diag: Optional[np.ndarray] = None,
    mode=1,
    factor: float = 100.0,
    nprint: int = 0,
    fjac: Optional[np.ndarray] = None,
    ipvt: Optional[np.ndarray] = None,
    qtf: Optional[np.ndarray] = None,
    **kwargs,
) -> int:
    """
    Minimize the sum of the squares of m nonlinear functions in n variables
    by a modification of the Levenberg-Marquardt algorithm.
    The user must provide a subroutine which calculates the functions.
    The jacobian is then calculated by a forward-difference approximation.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSolutionTolerance
        When no further improvement in the approximate solution is possible.
    MinpackJacobianTolerance
        The solution is orthogonal to the jacobian.
    """
    n = x.size
    m = fvec.size
    info = ffi.new("int *")
    nfev = ffi.new("int *")
    if maxfev is None:
        maxfev = 200 * (n + 1)
    if diag is None:
        diag = np.ones(n, dtype=real)
    if fjac is None:
        fjac = np.zeros((n, m), dtype=real)
    if ipvt is None:
        ipvt = np.zeros(n, dtype=np.int32)
    if qtf is None:
        qtf = np.zeros(n, dtype=real)
    wa1 = np.zeros(n, dtype=real)
    wa2 = np.zeros(n, dtype=real)
    wa3 = np.zeros(n, dtype=real)
    wa4 = np.zeros(m, dtype=real)

    minpack_lmdif(
        fcn,
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ftol,
        xtol,
        gtol,
        maxfev,
        epsfcn,
        ffi.cast("double*", diag.ctypes.data),
        mode,
        factor,
        nprint,
        info,
        nfev,
        ffi.cast("double*", fjac.ctypes.data),
        m,
        ffi.cast("int*", ipvt.ctypes.data),
        ffi.cast("double*", qtf.ctypes.data),
        ffi.cast("double*", wa1.ctypes.data),
        ffi.cast("double*", wa2.ctypes.data),
        ffi.cast("double*", wa3.ctypes.data),
        ffi.cast("double*", wa4.ctypes.data),
        **kwargs,
    )
    ex = info_lm(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def lmstr1(
    fcn: CallableLmstr,
    x: np.ndarray,
    fvec: np.ndarray,
    fjac: np.ndarray,
    tol: float = math.sqrt(np.finfo(real).eps),
    **kwargs,
) -> int:
    """
    Minimize the sum of the squares of m nonlinear functions in n variables by
    a modification of the Levenberg-Marquardt algorithm which uses minimal storage.
    This is done by using the more general least-squares solver `lmstr`.
    The user must provide a subroutine which calculates the functions and
    the rows of the Jacobian.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSolutionTolerance
        When no further improvement in the approximate solution is possible.
    MinpackJacobianTolerance
        The solution is orthogonal to the jacobian.
    """
    n = x.size
    m = fvec.size
    lwa = m * n + 5 * n + m
    wa = np.zeros(lwa, dtype=real)
    ipvt = np.zeros(n, dtype=np.int32)
    info = ffi.new("int *")

    minpack_lmstr1(
        fcn,
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        n,
        tol,
        info,
        ffi.cast("int*", ipvt.ctypes.data),
        ffi.cast("double*", wa.ctypes.data),
        lwa,
        **kwargs,
    )
    ex = info_lm(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def lmstr(
    fcn: CallableLmstr,
    x: np.ndarray,
    fvec: np.ndarray,
    fjac: np.ndarray,
    ftol: float = math.sqrt(np.finfo(real).eps),
    xtol: float = math.sqrt(np.finfo(real).eps),
    *,
    gtol: float = 0.0,
    maxfev: Optional[int] = None,
    diag: Optional[np.ndarray] = None,
    mode: int = 1,
    factor: float = 100.0,
    nprint=0,
    ipvt: Optional[np.ndarray] = None,
    qtf: Optional[np.ndarray] = None,
    **kwargs,
) -> int:
    """
    Minimize the sum of the squares of m nonlinear functions in n variables by
    a modification of the Levenberg-Marquardt algorithm which uses minimal storage.
    The user must provide a subroutine which calculates the functions and
    the rows of the Jacobian.

    Raises
    ------
    MinpackInputError
        In case of invalid input parameters.
    MinpackMaxIterations
        When the maximum number of iterations is exceeded.
    MinpackFunctionTolerance
        When the function tolerance cannot be satisfied.
    MinpackSolutionTolerance
        When no further improvement in the approximate solution is possible.
    MinpackJacobianTolerance
        The solution is orthogonal to the jacobian.
    """
    n = x.size
    m = fvec.size
    info = ffi.new("int *")
    nfev = ffi.new("int *")
    njev = ffi.new("int *")
    if maxfev is None:
        maxfev = 100 * (n + 1)
    if diag is None:
        diag = np.ones(n, dtype=real)
    if ipvt is None:
        ipvt = np.zeros(n, dtype=np.int32)
    if qtf is None:
        qtf = np.zeros(n, dtype=real)
    wa1 = np.zeros(n, dtype=real)
    wa2 = np.zeros(n, dtype=real)
    wa3 = np.zeros(n, dtype=real)
    wa4 = np.zeros(m, dtype=real)

    minpack_lmstr(
        fcn,
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        n,
        ftol,
        xtol,
        gtol,
        maxfev,
        ffi.cast("double*", diag.ctypes.data),
        mode,
        factor,
        nprint,
        info,
        nfev,
        njev,
        ffi.cast("int*", ipvt.ctypes.data),
        ffi.cast("double*", qtf.ctypes.data),
        ffi.cast("double*", wa1.ctypes.data),
        ffi.cast("double*", wa2.ctypes.data),
        ffi.cast("double*", wa3.ctypes.data),
        ffi.cast("double*", wa4.ctypes.data),
        **kwargs,
    )
    ex = info_lm(info[0])
    if ex is not None:
        raise ex(ex.__doc__)
    return info[0]


def chkder(x, fvec, fjac, xp, fvecp, check, error):
    """
    This subroutine checks the gradients of m nonlinear functions
    in n variables, evaluated at a point x, for consistency with
    the functions themselves.

    The subroutine does not perform reliably if cancellation or
    rounding errors cause a severe loss of significance in the
    evaluation of a function. Therefore, none of the components
    of x should be unusually small (in particular, zero) or any
    other value which may cause loss of significance.
    """
    if not fvec.size == fjac.shape[-1] == fvecp.size == error.size:
        raise ValueError("fvec, fjac, fvecp, error must have the same size")
    if not x.size == fjac.shape[0] == xp.size:
        raise ValueError("x, fjac, xp must have the same size")

    m = fvec.size
    n = x.size
    ldfjac = fjac.shape[-1]

    minpack_chkder(
        m,
        n,
        ffi.cast("double*", x.ctypes.data),
        ffi.cast("double*", fvec.ctypes.data),
        ffi.cast("double*", fjac.ctypes.data),
        ldfjac,
        ffi.cast("double*", xp.ctypes.data),
        ffi.cast("double*", fvecp.ctypes.data),
        2 if check else 1,
        ffi.cast("double*", error.ctypes.data),
    )
