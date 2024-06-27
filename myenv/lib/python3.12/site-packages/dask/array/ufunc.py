from __future__ import annotations

from functools import partial
from operator import getitem

import numpy as np

from dask import core
from dask.array.core import Array, apply_infer_dtype, asarray, blockwise, elemwise
from dask.base import is_dask_collection, normalize_token
from dask.highlevelgraph import HighLevelGraph
from dask.utils import derived_from, funcname


def wrap_elemwise(numpy_ufunc, source=np):
    """Wrap up numpy function into dask.array"""

    def wrapped(*args, **kwargs):
        dsk = [arg for arg in args if hasattr(arg, "_elemwise")]
        if len(dsk) > 0:
            return dsk[0]._elemwise(numpy_ufunc, *args, **kwargs)
        else:
            return numpy_ufunc(*args, **kwargs)

    # functools.wraps cannot wrap ufunc in Python 2.x
    wrapped.__name__ = numpy_ufunc.__name__
    return derived_from(source)(wrapped)


class da_frompyfunc:
    """A serializable `frompyfunc` object"""

    def __init__(self, func, nin, nout):
        self._ufunc = np.frompyfunc(func, nin, nout)
        self._func = func
        self.nin = nin
        self.nout = nout
        self._name = funcname(func)
        self.__name__ = "frompyfunc-%s" % self._name

    def __repr__(self):
        return "da.frompyfunc<%s, %d, %d>" % (self._name, self.nin, self.nout)

    def __dask_tokenize__(self):
        return (normalize_token(self._func), self.nin, self.nout)

    def __reduce__(self):
        return (da_frompyfunc, (self._func, self.nin, self.nout))

    def __call__(self, *args, **kwargs):
        return self._ufunc(*args, **kwargs)

    def __getattr__(self, a):
        if not a.startswith("_"):
            return getattr(self._ufunc, a)
        raise AttributeError(f"{type(self).__name__!r} object has no attribute {a!r}")

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        o.update(dir(self._ufunc))
        return list(o)


@derived_from(np)
def frompyfunc(func, nin, nout):
    if nout > 1:
        raise NotImplementedError("frompyfunc with more than one output")
    return ufunc(da_frompyfunc(func, nin, nout))


class ufunc:
    _forward_attrs = {
        "nin",
        "nargs",
        "nout",
        "ntypes",
        "identity",
        "signature",
        "types",
    }

    def __init__(self, ufunc):
        if not isinstance(ufunc, (np.ufunc, da_frompyfunc)):
            raise TypeError(
                "must be an instance of `ufunc` or "
                "`da_frompyfunc`, got `%s" % type(ufunc).__name__
            )
        self._ufunc = ufunc
        self.__name__ = ufunc.__name__
        if isinstance(ufunc, np.ufunc):
            derived_from(np)(self)

    def __dask_tokenize__(self):
        return self.__name__, normalize_token(self._ufunc)

    def __getattr__(self, key):
        if key in self._forward_attrs:
            return getattr(self._ufunc, key)
        raise AttributeError(f"{type(self).__name__!r} object has no attribute {key!r}")

    def __dir__(self):
        return list(self._forward_attrs.union(dir(type(self)), self.__dict__))

    def __repr__(self):
        return repr(self._ufunc)

    def __call__(self, *args, **kwargs):
        dsks = [arg for arg in args if hasattr(arg, "_elemwise")]
        if len(dsks) > 0:
            for dsk in dsks:
                result = dsk._elemwise(self._ufunc, *args, **kwargs)
                if type(result) != type(NotImplemented):
                    return result
            raise TypeError(
                "Parameters of such types are not supported by " + self.__name__
            )
        else:
            return self._ufunc(*args, **kwargs)

    @derived_from(np.ufunc)
    def outer(self, A, B, **kwargs):
        if self.nin != 2:
            raise ValueError("outer product only supported for binary functions")
        if "out" in kwargs:
            raise ValueError("`out` kwarg not supported")

        A_is_dask = is_dask_collection(A)
        B_is_dask = is_dask_collection(B)
        if not A_is_dask and not B_is_dask:
            return self._ufunc.outer(A, B, **kwargs)
        elif (
            A_is_dask
            and not isinstance(A, Array)
            or B_is_dask
            and not isinstance(B, Array)
        ):
            raise NotImplementedError(
                "Dask objects besides `dask.array.Array` "
                "are not supported at this time."
            )

        A = asarray(A)
        B = asarray(B)
        ndim = A.ndim + B.ndim
        out_inds = tuple(range(ndim))
        A_inds = out_inds[: A.ndim]
        B_inds = out_inds[A.ndim :]

        dtype = apply_infer_dtype(
            self._ufunc.outer, [A, B], kwargs, "ufunc.outer", suggest_dtype=False
        )

        if "dtype" in kwargs:
            func = partial(self._ufunc.outer, dtype=kwargs.pop("dtype"))
        else:
            func = self._ufunc.outer

        return blockwise(
            func,
            out_inds,
            A,
            A_inds,
            B,
            B_inds,
            dtype=dtype,
            token=self.__name__ + ".outer",
            **kwargs,
        )


# ufuncs, copied from this page:
# https://docs.scipy.org/doc/numpy/reference/ufuncs.html

# math operations
add = ufunc(np.add)
subtract = ufunc(np.subtract)
multiply = ufunc(np.multiply)
divide = ufunc(np.divide)
logaddexp = ufunc(np.logaddexp)
logaddexp2 = ufunc(np.logaddexp2)
true_divide = ufunc(np.true_divide)
floor_divide = ufunc(np.floor_divide)
negative = ufunc(np.negative)
positive = ufunc(np.positive)
power = ufunc(np.power)
float_power = ufunc(np.float_power)
remainder = ufunc(np.remainder)
mod = ufunc(np.mod)
# fmod: see below
conj = conjugate = ufunc(np.conjugate)
exp = ufunc(np.exp)
exp2 = ufunc(np.exp2)
log = ufunc(np.log)
log2 = ufunc(np.log2)
log10 = ufunc(np.log10)
log1p = ufunc(np.log1p)
expm1 = ufunc(np.expm1)
sqrt = ufunc(np.sqrt)
square = ufunc(np.square)
cbrt = ufunc(np.cbrt)
reciprocal = ufunc(np.reciprocal)

# trigonometric functions
sin = ufunc(np.sin)
cos = ufunc(np.cos)
tan = ufunc(np.tan)
arcsin = ufunc(np.arcsin)
arccos = ufunc(np.arccos)
arctan = ufunc(np.arctan)
arctan2 = ufunc(np.arctan2)
hypot = ufunc(np.hypot)
sinh = ufunc(np.sinh)
cosh = ufunc(np.cosh)
tanh = ufunc(np.tanh)
arcsinh = ufunc(np.arcsinh)
arccosh = ufunc(np.arccosh)
arctanh = ufunc(np.arctanh)
deg2rad = ufunc(np.deg2rad)
rad2deg = ufunc(np.rad2deg)

# comparison functions
greater = ufunc(np.greater)
greater_equal = ufunc(np.greater_equal)
less = ufunc(np.less)
less_equal = ufunc(np.less_equal)
not_equal = ufunc(np.not_equal)
equal = ufunc(np.equal)
isneginf = partial(equal, -np.inf)
isposinf = partial(equal, np.inf)
logical_and = ufunc(np.logical_and)
logical_or = ufunc(np.logical_or)
logical_xor = ufunc(np.logical_xor)
logical_not = ufunc(np.logical_not)
maximum = ufunc(np.maximum)
minimum = ufunc(np.minimum)
fmax = ufunc(np.fmax)
fmin = ufunc(np.fmin)

# bitwise functions
bitwise_and = ufunc(np.bitwise_and)
bitwise_or = ufunc(np.bitwise_or)
bitwise_xor = ufunc(np.bitwise_xor)
bitwise_not = ufunc(np.bitwise_not)
invert = bitwise_not
left_shift = ufunc(np.left_shift)
right_shift = ufunc(np.right_shift)

# floating functions
isfinite = ufunc(np.isfinite)
isinf = ufunc(np.isinf)
isnan = ufunc(np.isnan)
signbit = ufunc(np.signbit)
copysign = ufunc(np.copysign)
nextafter = ufunc(np.nextafter)
spacing = ufunc(np.spacing)
# modf: see below
ldexp = ufunc(np.ldexp)
# frexp: see below
fmod = ufunc(np.fmod)
floor = ufunc(np.floor)
ceil = ufunc(np.ceil)
trunc = ufunc(np.trunc)

# more math routines, from this page:
# https://docs.scipy.org/doc/numpy/reference/routines.math.html
degrees = ufunc(np.degrees)
radians = ufunc(np.radians)
rint = ufunc(np.rint)
fabs = ufunc(np.fabs)
sign = ufunc(np.sign)
absolute = ufunc(np.absolute)
abs = absolute

# non-ufunc elementwise functions
clip = wrap_elemwise(np.clip)
isreal = wrap_elemwise(np.isreal)
iscomplex = wrap_elemwise(np.iscomplex)
real = wrap_elemwise(np.real)
imag = wrap_elemwise(np.imag)
fix = wrap_elemwise(np.fix)
i0 = wrap_elemwise(np.i0)
sinc = wrap_elemwise(np.sinc)
nan_to_num = wrap_elemwise(np.nan_to_num)


@derived_from(np)
def angle(x, deg=0):
    deg = bool(deg)
    if hasattr(x, "_elemwise"):
        return x._elemwise(np.angle, x, deg)
    return np.angle(x, deg=deg)


@derived_from(np)
def frexp(x):
    # Not actually object dtype, just need to specify something
    tmp = elemwise(np.frexp, x, dtype=object)
    left = "mantissa-" + tmp.name
    right = "exponent-" + tmp.name
    ldsk = {
        (left,) + key[1:]: (getitem, key, 0)
        for key in core.flatten(tmp.__dask_keys__())
    }
    rdsk = {
        (right,) + key[1:]: (getitem, key, 1)
        for key in core.flatten(tmp.__dask_keys__())
    }

    a = np.empty_like(getattr(x, "_meta", x), shape=(1,) * x.ndim, dtype=x.dtype)
    l, r = np.frexp(a)

    graph = HighLevelGraph.from_collections(left, ldsk, dependencies=[tmp])
    L = Array(graph, left, chunks=tmp.chunks, meta=l)
    graph = HighLevelGraph.from_collections(right, rdsk, dependencies=[tmp])
    R = Array(graph, right, chunks=tmp.chunks, meta=r)
    return L, R


@derived_from(np)
def modf(x):
    # Not actually object dtype, just need to specify something
    tmp = elemwise(np.modf, x, dtype=object)
    left = "modf1-" + tmp.name
    right = "modf2-" + tmp.name
    ldsk = {
        (left,) + key[1:]: (getitem, key, 0)
        for key in core.flatten(tmp.__dask_keys__())
    }
    rdsk = {
        (right,) + key[1:]: (getitem, key, 1)
        for key in core.flatten(tmp.__dask_keys__())
    }

    a = np.ones_like(getattr(x, "_meta", x), shape=(1,) * x.ndim, dtype=x.dtype)
    l, r = np.modf(a)

    graph = HighLevelGraph.from_collections(left, ldsk, dependencies=[tmp])
    L = Array(graph, left, chunks=tmp.chunks, meta=l)
    graph = HighLevelGraph.from_collections(right, rdsk, dependencies=[tmp])
    R = Array(graph, right, chunks=tmp.chunks, meta=r)
    return L, R


@derived_from(np)
def divmod(x, y):
    res1 = x // y
    res2 = x % y
    return res1, res2
