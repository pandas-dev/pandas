import cmath

import numpy as np

from numba import float32
from numba.types import unicode_type, i8
from numba.pycc import CC, exportmany, export
from numba.tests.support import has_blas
from numba import typed


#
# New API
#

cc = CC('pycc_test_simple')
cc.use_nrt = False

# Note the first signature omits the return type
@cc.export('multf', (float32, float32))
@cc.export('multi', 'i4(i4, i4)')
def mult(a, b):
    return a * b

# Test imported C globals such as Py_None, PyExc_ZeroDivisionError
@cc.export('get_none', 'none()')
def get_none():
    return None

@cc.export('div', 'f8(f8, f8)')
def div(x, y):
    return x / y

_two = 2

# This one can't be compiled by the legacy API as it doesn't execute
# the script in a proper module.
@cc.export('square', 'i8(i8)')
def square(u):
    return u ** _two

# These ones need helperlib
cc_helperlib = CC('pycc_test_helperlib')
cc_helperlib.use_nrt = False

@cc_helperlib.export('power', 'i8(i8, i8)')
def power(u, v):
    return u ** v

@cc_helperlib.export('sqrt', 'c16(c16)')
def sqrt(u):
    return cmath.sqrt(u)

@cc_helperlib.export('size', 'i8(f8[:])')
def size(arr):
    return arr.size

# Exercise linking to Numpy math functions
@cc_helperlib.export('np_sqrt', 'f8(f8)')
def np_sqrt(u):
    return np.sqrt(u)

@cc_helperlib.export('spacing', 'f8(f8)')
def np_spacing(u):
    return np.spacing(u)


# This one clashes with libc random() unless pycc is careful with naming.
@cc_helperlib.export('random', 'f8(i4)')
def random_impl(seed):
    if seed != -1:
        np.random.seed(seed)
    return np.random.random()

# These ones need NRT
cc_nrt = CC('pycc_test_nrt')

@cc_nrt.export('zero_scalar', 'f8(i4)')
def zero_scalar(n):
    arr = np.zeros(n)
    return arr[-1]

if has_blas:
    # This one also needs BLAS
    @cc_nrt.export('vector_dot', 'f8(i4)')
    def vector_dot(n):
        a = np.linspace(1, n, n)
        return np.dot(a, a)

# This one needs an environment
@cc_nrt.export('zeros', 'f8[:](i4)')
def zeros(n):
    return np.zeros(n)

# requires list dtor, #issue3535
@cc_nrt.export('np_argsort', 'intp[:](float64[:])')
def np_argsort(arr):
    return np.argsort(arr)

#
# Legacy API
#

exportmany(['multf f4(f4,f4)', 'multi i4(i4,i4)'])(mult)
# Needs to link to helperlib to due with complex arguments
# export('multc c16(c16,c16)')(mult)
export('mult f8(f8, f8)')(mult)


@cc_nrt.export('dict_usecase', 'intp[:](intp[:])')
def dict_usecase(arr):
    d = typed.Dict()
    for i in range(arr.size):
        d[i] = arr[i]
    out = np.zeros_like(arr)
    for k, v in d.items():
        out[k] = k * v
    return out

# checks for issue #6386
@cc_nrt.export('internal_str_dict', i8(unicode_type))
def internal_str_dict(x):
    d = typed.Dict.empty(unicode_type,i8)
    if(x not in d):
        d[x] = len(d)
    return len(d)

@cc_nrt.export('hash_str', i8(unicode_type))
def internal_str_dict(x):
    return hash(x)

@cc_nrt.export('hash_literal_str_A', i8())
def internal_str_dict():
    return hash("A")
