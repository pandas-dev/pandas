import inspect

from numba.np.ufunc import _internal
from numba.np.ufunc.parallel import ParallelUFuncBuilder, ParallelGUFuncBuilder

from numba.core.registry import DelayedRegistry
from numba.np.ufunc import dufunc
from numba.np.ufunc import gufunc


class _BaseVectorize(object):

    @classmethod
    def get_identity(cls, kwargs):
        return kwargs.pop('identity', None)

    @classmethod
    def get_cache(cls, kwargs):
        return kwargs.pop('cache', False)

    @classmethod
    def get_writable_args(cls, kwargs):
        return kwargs.pop('writable_args', ())

    @classmethod
    def get_target_implementation(cls, kwargs):
        target = kwargs.pop('target', 'cpu')
        try:
            return cls.target_registry[target]
        except KeyError:
            raise ValueError("Unsupported target: %s" % target)


class Vectorize(_BaseVectorize):
    target_registry = DelayedRegistry({'cpu': dufunc.DUFunc,
                                       'parallel': ParallelUFuncBuilder,})

    def __new__(cls, func, **kws):
        identity = cls.get_identity(kws)
        cache = cls.get_cache(kws)
        imp = cls.get_target_implementation(kws)
        return imp(func, identity=identity, cache=cache, targetoptions=kws)


class GUVectorize(_BaseVectorize):
    target_registry = DelayedRegistry({'cpu': gufunc.GUFunc,
                                       'parallel': ParallelGUFuncBuilder,})

    def __new__(cls, func, signature, **kws):
        identity = cls.get_identity(kws)
        cache = cls.get_cache(kws)
        imp = cls.get_target_implementation(kws)
        writable_args = cls.get_writable_args(kws)
        if imp is gufunc.GUFunc:
            is_dyn = kws.pop('is_dynamic', False)
            return imp(func, signature, identity=identity, cache=cache,
                       is_dynamic=is_dyn, targetoptions=kws,
                       writable_args=writable_args)
        else:
            return imp(func, signature, identity=identity, cache=cache,
                       targetoptions=kws, writable_args=writable_args)


def vectorize(ftylist_or_function=(), **kws):
    """vectorize(ftylist_or_function=(), target='cpu', identity=None, **kws)

    A decorator that creates a NumPy ufunc object using Numba compiled
    code.  When no arguments or only keyword arguments are given,
    vectorize will return a Numba dynamic ufunc (DUFunc) object, where
    compilation/specialization may occur at call-time.

    Args
    -----
    ftylist_or_function: function or iterable

        When the first argument is a function, signatures are dealt
        with at call-time.

        When the first argument is an iterable of type signatures,
        which are either function type object or a string describing
        the function type, signatures are finalized at decoration
        time.

    Keyword Args
    ------------

    target: str
            A string for code generation target.  Default to "cpu".

    identity: int, str, or None
        The identity (or unit) value for the element-wise function
        being implemented.  Allowed values are None (the default), 0, 1,
        and "reorderable".

    cache: bool
        Turns on caching.


    Returns
    --------

    A NumPy universal function

    Examples
    -------
        @vectorize(['float32(float32, float32)',
                    'float64(float64, float64)'], identity=0)
        def sum(a, b):
            return a + b

        @vectorize
        def sum(a, b):
            return a + b

        @vectorize(identity=1)
        def mul(a, b):
            return a * b

    """
    if isinstance(ftylist_or_function, str):
        # Common user mistake
        ftylist = [ftylist_or_function]
    elif inspect.isfunction(ftylist_or_function):
        return dufunc.DUFunc(ftylist_or_function, **kws)
    elif ftylist_or_function is not None:
        ftylist = ftylist_or_function

    def wrap(func):
        vec = Vectorize(func, **kws)
        for sig in ftylist:
            vec.add(sig)
        if len(ftylist) > 0:
            vec.disable_compile()
        return vec.build_ufunc()

    return wrap


def guvectorize(*args, **kwargs):
    """guvectorize(ftylist, signature, target='cpu', identity=None, **kws)

    A decorator to create NumPy generalized-ufunc object from Numba compiled
    code.

    Args
    -----
    ftylist: iterable
        An iterable of type signatures, which are either
        function type object or a string describing the
        function type.

    signature: str
        A NumPy generalized-ufunc signature.
        e.g. "(m, n), (n, p)->(m, p)"

    identity: int, str, or None
        The identity (or unit) value for the element-wise function
        being implemented.  Allowed values are None (the default), 0, 1,
        and "reorderable".

    cache: bool
        Turns on caching.

    writable_args: tuple
        a tuple of indices of input variables that are writable.

    target: str
            A string for code generation target.  Defaults to "cpu".

    Returns
    --------

    A NumPy generalized universal-function

    Example
    -------
        @guvectorize(['void(int32[:,:], int32[:,:], int32[:,:])',
                      'void(float32[:,:], float32[:,:], float32[:,:])'],
                      '(x, y),(x, y)->(x, y)')
        def add_2d_array(a, b, c):
            for i in range(c.shape[0]):
                for j in range(c.shape[1]):
                    c[i, j] = a[i, j] + b[i, j]

    """
    if len(args) == 1:
        ftylist = []
        signature = args[0]
        kwargs.setdefault('is_dynamic', True)
    elif len(args) == 2:
        ftylist = args[0]
        signature = args[1]
    else:
        raise TypeError('guvectorize() takes one or two positional arguments')

    if isinstance(ftylist, str):
        # Common user mistake
        ftylist = [ftylist]

    def wrap(func):
        guvec = GUVectorize(func, signature, **kwargs)
        for fty in ftylist:
            guvec.add(fty)
        if len(ftylist) > 0:
            guvec.disable_compile()
        return guvec.build_ufunc()

    return wrap
