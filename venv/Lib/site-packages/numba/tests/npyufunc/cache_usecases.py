import numba as nb


#
# UFunc
#

def direct_ufunc_cache_usecase(**kwargs):
    @nb.vectorize(["intp(intp)", "float64(float64)"], cache=True, **kwargs)
    def ufunc(inp):
        return inp * 2

    return ufunc


def indirect_ufunc_cache_usecase(**kwargs):
    @nb.njit(cache=True)
    def indirect_ufunc_core(inp):
        return inp * 3

    @nb.vectorize(["intp(intp)", "float64(float64)", "complex64(complex64)"],
                  **kwargs)
    def ufunc(inp):
        return indirect_ufunc_core(inp)

    return ufunc


#
# DUFunc
#

def direct_dufunc_cache_usecase(**kwargs):
    @nb.vectorize(cache=True, **kwargs)
    def ufunc(inp):
        return inp * 2

    return ufunc


def indirect_dufunc_cache_usecase(**kwargs):
    @nb.njit(cache=True)
    def indirect_ufunc_core(inp):
        return inp * 3

    @nb.vectorize(**kwargs)
    def ufunc(inp):
        return indirect_ufunc_core(inp)

    return ufunc


#
# GUFunc
#

def direct_gufunc_cache_usecase(**kwargs):
    @nb.guvectorize(["(intp, intp[:])", "(float64, float64[:])"],
                    "()->()", cache=True, **kwargs)
    def gufunc(inp, out):
        out[0] = inp * 2

    return gufunc


def indirect_gufunc_cache_usecase(**kwargs):
    @nb.njit(cache=True)
    def core(x):
        return x * 3

    @nb.guvectorize(["(intp, intp[:])", "(float64, float64[:])",
                     "(complex64, complex64[:])"], "()->()", **kwargs)
    def gufunc(inp, out):
        out[0] = core(inp)

    return gufunc
