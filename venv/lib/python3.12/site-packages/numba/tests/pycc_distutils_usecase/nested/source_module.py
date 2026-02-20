import numpy as np

from numba.pycc import CC


cc = CC('pycc_compiled_module')

_const = 42


# This ones references a global variable at compile time
@cc.export('get_const', 'i8()')
def get_const():
    return _const


# This one needs NRT and an environment
@cc.export('ones', 'f8[:](i4)')
def ones(n):
    return np.ones(n)
