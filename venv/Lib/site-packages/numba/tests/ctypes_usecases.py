from ctypes import *
import sys

import numpy as np


is_windows = sys.platform.startswith('win32')

# We can't rely on libc availability on Windows anymore, so we use our
# own compiled wrappers (see https://bugs.python.org/issue23606).

from numba import _helperlib
libnumba = CDLL(_helperlib.__file__)
del _helperlib

# A typed C function (cdecl under Windows)

c_sin = libnumba._numba_test_sin
c_sin.argtypes = [c_double]
c_sin.restype = c_double

def use_c_sin(x):
    return c_sin(x)

c_cos = libnumba._numba_test_cos
c_cos.argtypes = [c_double]
c_cos.restype = c_double

def use_two_funcs(x):
    return c_sin(x) - c_cos(x)

# Typed C functions accepting an array pointer
# (either as a "void *" or as a typed pointer)

c_vsquare = libnumba._numba_test_vsquare
c_vsquare.argtypes = [c_int, c_void_p, c_void_p]

c_vcube = libnumba._numba_test_vsquare
c_vcube.argtypes = [c_int, POINTER(c_double), POINTER(c_double)]

def use_c_vsquare(x):
    out = np.empty_like(x)
    c_vsquare(x.size, x.ctypes, out.ctypes)
    return out

def use_c_vcube(x):
    out = np.empty_like(x)
    c_vcube(x.size, x.ctypes, out.ctypes)
    return out

# An untyped C function

c_untyped = libnumba._numba_test_exp

def use_c_untyped(x):
    return c_untyped(x)

# A C function wrapped in a CFUNCTYPE

ctype_wrapping = CFUNCTYPE(c_double, c_double)(use_c_sin)

def use_ctype_wrapping(x):
    return ctype_wrapping(x)

# A Python API function

savethread = pythonapi.PyEval_SaveThread
savethread.argtypes = []
savethread.restype = c_void_p

restorethread = pythonapi.PyEval_RestoreThread
restorethread.argtypes = [c_void_p]
restorethread.restype = None

if is_windows:
    # A function with the stdcall calling convention
    c_sleep = windll.kernel32.Sleep
    c_sleep.argtypes = [c_uint]
    c_sleep.restype = None

    def use_c_sleep(x):
        c_sleep(x)


def use_c_pointer(x):
    """
    Running in Python will cause a segfault.
    """
    threadstate = savethread()
    x += 1
    restorethread(threadstate)
    return x


def use_func_pointer(fa, fb, x):
    if x > 0:
        return fa(x)
    else:
        return fb(x)


mydct = {'what': 1232121}

def call_me_maybe(arr):
    return mydct[arr[0].decode('ascii')]

# Create a callback into the python interpreter
py_call_back = CFUNCTYPE(c_int, py_object)(call_me_maybe)


def take_array_ptr(ptr):
    return ptr

c_take_array_ptr = CFUNCTYPE(c_void_p, c_void_p)(take_array_ptr)
