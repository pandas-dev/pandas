"""
Isolate pandas's exposure to NumPy
"""

import numpy as np

Array = np.ndarray

bool = np.bool_

_dtypes = {
    'int': [8, 16, 32, 64],
    'uint': [8, 16, 32, 64],
    'float': [16, 32, 64]
}

_lift_types = []

for _k, _v in _dtypes.items():
    for _i in _v:
        _lift_types.append(_k + str(_i))

for _t in _lift_types:
    globals()[_t] = getattr(np, _t)

_lift_function = ['empty', 'arange', 'array', 'putmask', 'where']

for _f in _lift_function:
    globals()[_f] = getattr(np, _f)

_lift_random = ['randn', 'rand']

for _f in _lift_random:
    globals()[_f] = getattr(np.random, _f)

NA = np.nan

