'''
Implements the cuda module as called from within an executing kernel
(@cuda.jit-decorated function).
'''

from contextlib import contextmanager
import sys
import threading
import traceback
from numba.core import types
import numpy as np

from numba.np import numpy_support

from .vector_types import vector_types


class Dim3(object):
    '''
    Used to implement thread/block indices/dimensions
    '''
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return '(%s, %s, %s)' % (self.x, self.y, self.z)

    def __repr__(self):
        return 'Dim3(%s, %s, %s)' % (self.x, self.y, self.z)

    def __iter__(self):
        yield self.x
        yield self.y
        yield self.z


class GridGroup:
    '''
    Used to implement the grid group.
    '''

    def sync(self):
        # Synchronization of the grid group is equivalent to synchronization of
        # the thread block, because we only support cooperative grids with one
        # block.
        threading.current_thread().syncthreads()


class FakeCUDACg:
    '''
    CUDA Cooperative Groups
    '''
    def this_grid(self):
        return GridGroup()


class FakeCUDALocal(object):
    '''
    CUDA Local arrays
    '''
    def array(self, shape, dtype):
        if isinstance(dtype, types.Type):
            dtype = numpy_support.as_dtype(dtype)
        return np.empty(shape, dtype)


class FakeCUDAConst(object):
    '''
    CUDA Const arrays
    '''
    def array_like(self, ary):
        return ary


class FakeCUDAShared(object):
    '''
    CUDA Shared arrays.

    Limitations: assumes that only one call to cuda.shared.array is on a line,
    and that that line is only executed once per thread. i.e.::

        a = cuda.shared.array(...); b = cuda.shared.array(...)

    will erroneously alias a and b, and::

        for i in range(10):
            sharedarrs[i] = cuda.shared.array(...)

    will alias all arrays created at that point (though it is not certain that
    this would be supported by Numba anyway).
    '''

    def __init__(self, dynshared_size):
        self._allocations = {}
        self._dynshared_size = dynshared_size
        self._dynshared = np.zeros(dynshared_size, dtype=np.byte)

    def array(self, shape, dtype):
        if isinstance(dtype, types.Type):
            dtype = numpy_support.as_dtype(dtype)
        # Dynamic shared memory is requested with size 0 - this all shares the
        # same underlying memory
        if shape == 0:
            # Count must be the maximum number of whole elements that fit in the
            # buffer (Numpy complains if the buffer is not a multiple of the
            # element size)
            count = self._dynshared_size // dtype.itemsize
            return np.frombuffer(self._dynshared.data, dtype=dtype, count=count)

        # Otherwise, identify allocations by source file and line number
        # We pass the reference frame explicitly to work around
        # http://bugs.python.org/issue25108
        stack = traceback.extract_stack(sys._getframe())
        caller = stack[-2][0:2]
        res = self._allocations.get(caller)
        if res is None:
            res = np.empty(shape, dtype)
            self._allocations[caller] = res
        return res


addlock = threading.Lock()
sublock = threading.Lock()
andlock = threading.Lock()
orlock = threading.Lock()
xorlock = threading.Lock()
maxlock = threading.Lock()
minlock = threading.Lock()
compare_and_swaplock = threading.Lock()
caslock = threading.Lock()
inclock = threading.Lock()
declock = threading.Lock()
exchlock = threading.Lock()


class FakeCUDAAtomic(object):
    def add(self, array, index, val):
        with addlock:
            old = array[index]
            array[index] += val
        return old

    def sub(self, array, index, val):
        with sublock:
            old = array[index]
            array[index] -= val
        return old

    def and_(self, array, index, val):
        with andlock:
            old = array[index]
            array[index] &= val
        return old

    def or_(self, array, index, val):
        with orlock:
            old = array[index]
            array[index] |= val
        return old

    def xor(self, array, index, val):
        with xorlock:
            old = array[index]
            array[index] ^= val
        return old

    def inc(self, array, index, val):
        with inclock:
            old = array[index]
            if old >= val:
                array[index] = 0
            else:
                array[index] += 1
        return old

    def dec(self, array, index, val):
        with declock:
            old = array[index]
            if (old == 0) or (old > val):
                array[index] = val
            else:
                array[index] -= 1
        return old

    def exch(self, array, index, val):
        with exchlock:
            old = array[index]
            array[index] = val
        return old

    def max(self, array, index, val):
        with maxlock:
            old = array[index]
            array[index] = max(old, val)
        return old

    def min(self, array, index, val):
        with minlock:
            old = array[index]
            array[index] = min(old, val)
        return old

    def nanmax(self, array, index, val):
        with maxlock:
            old = array[index]
            array[index] = np.nanmax([array[index], val])
        return old

    def nanmin(self, array, index, val):
        with minlock:
            old = array[index]
            array[index] = np.nanmin([array[index], val])
        return old

    def compare_and_swap(self, array, old, val):
        with compare_and_swaplock:
            index = (0,) * array.ndim
            loaded = array[index]
            if loaded == old:
                array[index] = val
            return loaded

    def cas(self, array, index, old, val):
        with caslock:
            loaded = array[index]
            if loaded == old:
                array[index] = val
            return loaded


class FakeCUDAFp16(object):
    def hadd(self, a, b):
        return a + b

    def hsub(self, a, b):
        return a - b

    def hmul(self, a, b):
        return a * b

    def hdiv(self, a, b):
        return a / b

    def hfma(self, a, b, c):
        return a * b + c

    def hneg(self, a):
        return -a

    def habs(self, a):
        return abs(a)

    def hsin(self, x):
        return np.sin(x, dtype=np.float16)

    def hcos(self, x):
        return np.cos(x, dtype=np.float16)

    def hlog(self, x):
        return np.log(x, dtype=np.float16)

    def hlog2(self, x):
        return np.log2(x, dtype=np.float16)

    def hlog10(self, x):
        return np.log10(x, dtype=np.float16)

    def hexp(self, x):
        return np.exp(x, dtype=np.float16)

    def hexp2(self, x):
        return np.exp2(x, dtype=np.float16)

    def hexp10(self, x):
        return np.float16(10 ** x)

    def hsqrt(self, x):
        return np.sqrt(x, dtype=np.float16)

    def hrsqrt(self, x):
        return np.float16(x ** -0.5)

    def hceil(self, x):
        return np.ceil(x, dtype=np.float16)

    def hfloor(self, x):
        return np.ceil(x, dtype=np.float16)

    def hrcp(self, x):
        return np.reciprocal(x, dtype=np.float16)

    def htrunc(self, x):
        return np.trunc(x, dtype=np.float16)

    def hrint(self, x):
        return np.rint(x, dtype=np.float16)

    def heq(self, a, b):
        return a == b

    def hne(self, a, b):
        return a != b

    def hge(self, a, b):
        return a >= b

    def hgt(self, a, b):
        return a > b

    def hle(self, a, b):
        return a <= b

    def hlt(self, a, b):
        return a < b

    def hmax(self, a, b):
        return max(a, b)

    def hmin(self, a, b):
        return min(a, b)


class FakeCUDAModule(object):
    '''
    An instance of this class will be injected into the __globals__ for an
    executing function in order to implement calls to cuda.*. This will fail to
    work correctly if the user code does::

        from numba import cuda as something_else

    In other words, the CUDA module must be called cuda.
    '''

    def __init__(self, grid_dim, block_dim, dynshared_size):
        self.gridDim = Dim3(*grid_dim)
        self.blockDim = Dim3(*block_dim)
        self._cg = FakeCUDACg()
        self._local = FakeCUDALocal()
        self._shared = FakeCUDAShared(dynshared_size)
        self._const = FakeCUDAConst()
        self._atomic = FakeCUDAAtomic()
        self._fp16 = FakeCUDAFp16()
        # Insert the vector types into the kernel context
        # Note that we need to do this in addition to exposing them as module
        # variables in `simulator.__init__.py`, because the test cases need
        # to access the actual cuda module as well as the fake cuda module
        # for vector types.
        for name, svty in vector_types.items():
            setattr(self, name, svty)
            for alias in svty.aliases:
                setattr(self, alias, svty)

    @property
    def cg(self):
        return self._cg

    @property
    def local(self):
        return self._local

    @property
    def shared(self):
        return self._shared

    @property
    def const(self):
        return self._const

    @property
    def atomic(self):
        return self._atomic

    @property
    def fp16(self):
        return self._fp16

    @property
    def threadIdx(self):
        return threading.current_thread().threadIdx

    @property
    def blockIdx(self):
        return threading.current_thread().blockIdx

    @property
    def warpsize(self):
        return 32

    @property
    def laneid(self):
        return threading.current_thread().thread_id % 32

    def syncthreads(self):
        threading.current_thread().syncthreads()

    def threadfence(self):
        # No-op
        pass

    def threadfence_block(self):
        # No-op
        pass

    def threadfence_system(self):
        # No-op
        pass

    def syncthreads_count(self, val):
        return threading.current_thread().syncthreads_count(val)

    def syncthreads_and(self, val):
        return threading.current_thread().syncthreads_and(val)

    def syncthreads_or(self, val):
        return threading.current_thread().syncthreads_or(val)

    def popc(self, val):
        return bin(val).count("1")

    def fma(self, a, b, c):
        return a * b + c

    def cbrt(self, a):
        return a ** (1 / 3)

    def brev(self, val):
        return int('{:032b}'.format(val)[::-1], 2)

    def clz(self, val):
        s = '{:032b}'.format(val)
        return len(s) - len(s.lstrip('0'))

    def ffs(self, val):
        # The algorithm is:
        # 1. Count the number of trailing zeros.
        # 2. Add 1, because the LSB is numbered 1 rather than 0, and so on.
        # 3. If we've counted 32 zeros (resulting in 33), there were no bits
        #    set so we need to return zero.
        s = '{:032b}'.format(val)
        r = (len(s) - len(s.rstrip('0')) + 1) % 33
        return r

    def selp(self, a, b, c):
        return b if a else c

    def grid(self, n):
        bdim = self.blockDim
        bid = self.blockIdx
        tid = self.threadIdx
        x = bid.x * bdim.x + tid.x
        if n == 1:
            return x
        y = bid.y * bdim.y + tid.y
        if n == 2:
            return (x, y)
        z = bid.z * bdim.z + tid.z
        if n == 3:
            return (x, y, z)

        raise RuntimeError("Global ID has 1-3 dimensions. %d requested" % n)

    def gridsize(self, n):
        bdim = self.blockDim
        gdim = self.gridDim
        x = bdim.x * gdim.x
        if n == 1:
            return x
        y = bdim.y * gdim.y
        if n == 2:
            return (x, y)
        z = bdim.z * gdim.z
        if n == 3:
            return (x, y, z)

        raise RuntimeError("Global grid has 1-3 dimensions. %d requested" % n)


@contextmanager
def swapped_cuda_module(fn, fake_cuda_module):
    from numba import cuda

    fn_globs = fn.__globals__
    # get all globals that is the "cuda" module
    orig = dict((k, v) for k, v in fn_globs.items() if v is cuda)
    # build replacement dict
    repl = dict((k, fake_cuda_module) for k, v in orig.items())
    # replace
    fn_globs.update(repl)
    try:
        yield
    finally:
        # revert
        fn_globs.update(orig)
