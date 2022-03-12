# -*- coding: utf-8 -*-
"""
Callback signatures for Minpack drivers.
"""

import numpy as np
from typing import Optional, Callable

CallableHybrd = Callable[[np.ndarray, np.ndarray], None]
CallableHybrj = Callable[[np.ndarray, np.ndarray, np.ndarray, bool], None]
CallableLmder = Callable[[np.ndarray, np.ndarray, np.ndarray, bool], None]
CallableLmdif = Callable[[np.ndarray, np.ndarray], None]
CallableLmstr = Callable[[np.ndarray, np.ndarray, np.ndarray, Optional[int]], None]
