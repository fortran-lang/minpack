# -*- coding: utf-8 -*-
"""
Callback signatures for Minpack drivers.
"""

import numpy as np
from typing import Optional
try:
    from typing import Protocol
except ImportError:
    from typing_extensions import Protocol


class CallableHybrd(Protocol):
    def __call__(self, x: np.ndarray, fvec: np.ndarray, **kwargs) -> None:
        ...


class CallableHybrj(Protocol):
    def __call__(
        self,
        x: np.ndarray,
        fvec: np.ndarray,
        fjac: np.ndarray,
        jacobian: bool,
        **kwargs
    ) -> None:
        ...


class CallableLmder(Protocol):
    def __call__(
        self,
        x: np.ndarray,
        fvec: np.ndarray,
        fjac: np.ndarray,
        jacobian: bool,
        **kwargs
    ) -> None:
        ...


class CallableLmdif(Protocol):
    def __call__(self, x: np.ndarray, fvec: np.ndarray, **kwargs) -> None:
        ...


class CallableLmstr(Protocol):
    def __call__(
        self,
        x: np.ndarray,
        fvec: np.ndarray,
        fjrow: np.ndarray,
        row: Optional[int],
        **kwargs
    ) -> None:
        ...
