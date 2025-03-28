# -*- coding: utf-8 -*-

from numba.np.ufunc.decorators import Vectorize, GUVectorize, vectorize, guvectorize
from numba.np.ufunc._internal import PyUFunc_None, PyUFunc_Zero, PyUFunc_One
from numba.np.ufunc import _internal, array_exprs
from numba.np.ufunc.parallel import (threading_layer, get_num_threads,
                                     set_num_threads, get_thread_id,
                                     set_parallel_chunksize,
                                     get_parallel_chunksize)


if hasattr(_internal, 'PyUFunc_ReorderableNone'):
    PyUFunc_ReorderableNone = _internal.PyUFunc_ReorderableNone
del _internal, array_exprs


def _init():

    def init_cuda_vectorize():
        from numba.cuda.vectorizers import CUDAVectorize
        return CUDAVectorize

    def init_cuda_guvectorize():
        from numba.cuda.vectorizers import CUDAGUFuncVectorize
        return CUDAGUFuncVectorize

    Vectorize.target_registry.ondemand['cuda'] = init_cuda_vectorize
    GUVectorize.target_registry.ondemand['cuda'] = init_cuda_guvectorize


_init()
del _init
