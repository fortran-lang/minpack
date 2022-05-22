# -*- coding: utf-8 -*-

import pytest
import minpack.library
import numpy as np
from math import sqrt


@pytest.mark.parametrize("driver", [minpack.library.hybrd1, minpack.library.hybrd])
def test_hybrd(driver):
    def fcn(x, fvec) -> None:
        for k in range(x.size):
            tmp = (3.0 - 2.0 * x[k]) * x[k]
            tmp1 = x[k - 1] if k > 0 else 0.0
            tmp2 = x[k + 1] if k < len(x) - 1 else 0.0
            fvec[k] = tmp - tmp1 - 2.0 * tmp2 + 1.0

    x = np.array(9 * [-1.0])
    fvec = np.zeros(9, dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    assert driver(fcn, x, fvec, tol) == 1

    assert pytest.approx(x, abs=10 * tol) == [
        -0.5706545,
        -0.6816283,
        -0.7017325,
        -0.7042129,
        -0.7013690,
        -0.6918656,
        -0.6657920,
        -0.5960342,
        -0.4164121,
    ]


@pytest.mark.parametrize("driver", [minpack.library.hybrd1, minpack.library.hybrd])
def test_hybrd_exception(driver):
    class DummyException(Exception):
        ...

    def fcn(x, fvec) -> None:
        raise DummyException()

    x = np.array(9 * [-1.0])
    fvec = np.zeros(9, dtype=np.float64)

    with pytest.raises(DummyException):
        driver(fcn, x, fvec)


@pytest.mark.parametrize("driver", [minpack.library.hybrj1, minpack.library.hybrj])
def test_hybrj(driver):
    def fcn(x, fvec, fjac, jacobian: bool) -> None:

        if jacobian:
            for k in range(x.size):
                for j in range(x.size):
                    fjac[k, j] = 0.0
                fjac[k, k] = 3.0 - 4.0 * x[k]
                if k > 0:
                    fjac[k, k - 1] = -1.0
                if k < x.size - 1:
                    fjac[k, k + 1] = -2.0
        else:
            for k in range(x.size):
                tmp = (3.0 - 2.0 * x[k]) * x[k]
                tmp1 = x[k - 1] if k > 0 else 0.0
                tmp2 = x[k + 1] if k < len(x) - 1 else 0.0
                fvec[k] = tmp - tmp1 - 2.0 * tmp2 + 1.0

    x = np.array(9 * [-1.0])
    fvec = np.zeros(9, dtype=np.float64)
    fjac = np.zeros((9, 9), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    assert driver(fcn, x, fvec, fjac, tol) == 1

    assert pytest.approx(x, abs=10 * tol) == [
        -0.5706545,
        -0.6816283,
        -0.7017325,
        -0.7042129,
        -0.7013690,
        -0.6918656,
        -0.6657920,
        -0.5960342,
        -0.4164121,
    ]


@pytest.mark.parametrize("driver", [minpack.library.hybrj1, minpack.library.hybrj])
def test_hybrj_exception(driver):
    class DummyException(Exception):
        ...

    def fcn(x, fvec, fjac, jacobian) -> None:
        raise DummyException()

    x = np.array(9 * [-1.0])
    fvec = np.zeros(9, dtype=np.float64)
    fjac = np.zeros((9, 9), dtype=np.float64)

    with pytest.raises(DummyException):
        driver(fcn, x, fvec, fjac)


@pytest.mark.parametrize("driver", [minpack.library.lmder1, minpack.library.lmder])
def test_lmder(driver):
    y = np.array(
        [
            1.4e-1,
            1.8e-1,
            2.2e-1,
            2.5e-1,
            2.9e-1,
            3.2e-1,
            3.5e-1,
            3.9e-1,
            3.7e-1,
            5.8e-1,
            7.3e-1,
            9.6e-1,
            1.34e0,
            2.1e0,
            4.39e0,
        ]
    )

    def fcn(x, fvec, fjac, jacobian: bool, y) -> None:
        if jacobian:
            for i in range(fvec.size):
                tmp1, tmp2 = i + 1, 16 - i - 1
                tmp3 = tmp2 if i >= 8 else tmp1
                tmp4 = (x[1] * tmp2 + x[2] * tmp3) ** 2
                fjac[0, i] = -1.0
                fjac[1, i] = tmp1 * tmp2 / tmp4
                fjac[2, i] = tmp1 * tmp3 / tmp4
        else:
            for i in range(fvec.size):
                tmp1, tmp2 = i + 1, 16 - i - 1
                tmp3 = tmp2 if i >= 8 else tmp1
                fvec[i] = y[i] - (x[0] + tmp1 / (x[1] * tmp2 + x[2] * tmp3))

    x = np.array([1.0, 1.0, 1.0])
    fvec = np.zeros(15, dtype=np.float64)
    fjac = np.zeros((3, 15), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    xp = np.zeros(3, dtype=np.float64)
    fvecp = np.zeros(15, dtype=np.float64)
    err = np.zeros(15, dtype=np.float64)
    minpack.library.chkder(x, fvec, fjac, xp, fvecp, False, err)
    fcn(x, fvec, fjac, False, y=y)
    fcn(x, fvec, fjac, True, y=y)
    fcn(xp, fvecp, fjac, False, y=y)
    minpack.library.chkder(x, fvec, fjac, xp, fvecp, True, err)

    assert pytest.approx(err) == 15 * [1.0]

    assert driver(fcn, x, fvec, fjac, tol, y=y) == 1

    assert pytest.approx(x, abs=100 * tol) == [0.8241058e-1, 0.1133037e1, 0.2343695e1]


@pytest.mark.parametrize("driver", [minpack.library.lmder1, minpack.library.lmder])
def test_lmder_exception(driver):
    class DummyException(Exception):
        ...

    def fcn(x, fvec, fjac, jacobian: bool) -> None:
        raise DummyException()

    x = np.array([1.0, 1.0, 1.0])
    fvec = np.zeros(15, dtype=np.float64)
    fjac = np.zeros((3, 15), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    with pytest.raises(DummyException):
        driver(fcn, x, fvec, fjac, tol)


@pytest.mark.parametrize("driver", [minpack.library.lmdif1, minpack.library.lmdif])
def test_lmdif(driver):
    y = np.array(
        [
            1.4e-1,
            1.8e-1,
            2.2e-1,
            2.5e-1,
            2.9e-1,
            3.2e-1,
            3.5e-1,
            3.9e-1,
            3.7e-1,
            5.8e-1,
            7.3e-1,
            9.6e-1,
            1.34e0,
            2.1e0,
            4.39e0,
        ]
    )

    def fcn(x, fvec, y) -> None:
        for i in range(fvec.size):
            tmp1, tmp2 = i + 1, 16 - i - 1
            tmp3 = tmp2 if i >= 8 else tmp1
            fvec[i] = y[i] - (x[0] + tmp1 / (x[1] * tmp2 + x[2] * tmp3))

    x = np.array([1.0, 1.0, 1.0])
    fvec = np.zeros(15, dtype=np.float64)
    fjac = np.zeros((3, 15), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    assert driver(fcn, x, fvec, tol, y=y) == 1

    assert pytest.approx(x, abs=100 * tol) == [0.8241058e-1, 0.1133037e1, 0.2343695e1]


@pytest.mark.parametrize("driver", [minpack.library.lmdif1, minpack.library.lmdif])
def test_lmdif_exception(driver):
    class DummyException(Exception):
        ...

    def fcn(x, fvec) -> None:
        raise DummyException()

    x = np.array([1.0, 1.0, 1.0])
    fvec = np.zeros(15, dtype=np.float64)
    fjac = np.zeros((3, 15), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    with pytest.raises(DummyException):
        driver(fcn, x, fvec, tol)


@pytest.mark.parametrize("driver", [minpack.library.lmstr1, minpack.library.lmstr])
def test_lmstr(driver):
    def fcn(x, fvec, fjrow, row) -> None:
        if row is None:
            fvec[0] = 10.0 * (x[1] - x[0] ** 2)
            fvec[1] = 1.0 - x[0]
        else:
            fjrow[0] = -20.0 * x[0] if row == 0 else -1.0
            fjrow[1] = 10.0 if row == 0 else 0.0

    x = np.array([-1.2, 1.0])
    fvec = np.zeros(2, dtype=np.float64)
    fjac = np.zeros((2, 2), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    assert driver(fcn, x, fvec, fjac, tol) == 4

    assert pytest.approx(x, abs=100 * tol) == 2 * [1.0]


@pytest.mark.parametrize("driver", [minpack.library.lmstr1, minpack.library.lmstr])
def test_lmstr_exception(driver):
    class DummyException(Exception):
        ...

    def fcn(x, fvec, fjrow, row) -> None:
        raise DummyException()

    x = np.array([-1.2, 1.0])
    fvec = np.zeros(2, dtype=np.float64)
    fjac = np.zeros((2, 2), dtype=np.float64)
    tol = sqrt(np.finfo(np.float64).eps)

    with pytest.raises(DummyException):
        driver(fcn, x, fvec, fjac, tol)
