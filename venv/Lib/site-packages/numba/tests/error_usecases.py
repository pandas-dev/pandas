import numba as nb


@nb.jit(nopython=True, parallel=True)
def foo():
    pass
