"""
Separate module with function samples for serialization tests,
to avoid issues with __main__.
"""

import math

from numba import jit
from numba.core import types


@jit((types.int32, types.int32))
def add_with_sig(a, b):
    return a + b

@jit
def add_without_sig(a, b):
    return a + b

@jit(nopython=True)
def add_nopython(a, b):
    return a + b

@jit(nopython=True)
def add_nopython_fail(a, b):
    object()
    return a + b

def closure(a):
    @jit(nopython=True)
    def inner(b, c):
        return a + b + c
    return inner

K = 3.0

from math import sqrt

def closure_with_globals(x, **jit_args):
    @jit(**jit_args)
    def inner(y):
        # Exercise a builtin function and a module-level constant
        k = max(K, K + 1)
        # Exercise two functions from another module, one accessed with
        # dotted notation, one imported explicitly.
        return math.hypot(x, y) + sqrt(k)
    return inner

@jit(nopython=True)
def other_function(x, y):
    return math.hypot(x, y)

@jit(forceobj=True)
def get_global_objmode(x):
    return K * x

import numpy as np
import numpy.random as nprand

@jit(nopython=True)
def get_renamed_module(x):
    nprand.seed(42)
    return np.cos(x), nprand.random()


def closure_calling_other_function(x):
    @jit(nopython=True)
    def inner(y, z):
        return other_function(x, y) + z
    return inner

def closure_calling_other_closure(x):
    @jit(nopython=True)
    def other_inner(y):
        return math.hypot(x, y)

    @jit(nopython=True)
    def inner(y):
        return other_inner(y) + x
    return inner


# A dynamic function calling a builtin function
def _get_dyn_func(**jit_args):
    code = """
        def dyn_func(x):
            res = 0
            for i in range(x):
                res += x
            return res
        """
    ns = {}
    exec(code.strip(), ns)
    return jit(**jit_args)(ns['dyn_func'])

dyn_func = _get_dyn_func(nopython=True)
dyn_func_objmode = _get_dyn_func(forceobj=True)
