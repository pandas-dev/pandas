"""xarray specific universal functions."""

import textwrap
from abc import ABC, abstractmethod

import numpy as np

import xarray as xr
from xarray.core.groupby import GroupBy


def _walk_array_namespaces(obj, namespaces):
    if isinstance(obj, xr.DataTree):
        # TODO: DataTree doesn't actually support ufuncs yet
        for node in obj.subtree:
            _walk_array_namespaces(node.dataset, namespaces)
    elif isinstance(obj, xr.Dataset):
        for name in obj.data_vars:
            _walk_array_namespaces(obj[name], namespaces)
    elif isinstance(obj, GroupBy):
        _walk_array_namespaces(next(iter(obj))[1], namespaces)
    elif isinstance(obj, xr.DataArray | xr.Variable):
        _walk_array_namespaces(obj.data, namespaces)
    else:
        namespace = getattr(obj, "__array_namespace__", None)
        if namespace is not None:
            namespaces.add(namespace())

    return namespaces


def get_array_namespace(*args):
    xps = set()
    for arg in args:
        _walk_array_namespaces(arg, xps)

    xps.discard(np)
    if len(xps) > 1:
        names = [module.__name__ for module in xps]
        raise ValueError(f"Mixed array types {names} are not supported.")

    return next(iter(xps)) if xps else np


class _ufunc_wrapper(ABC):
    def __init__(self, name):
        self.__name__ = name
        if hasattr(np, name):
            self._create_doc()

    @abstractmethod
    def __call__(self, *args, **kwargs):
        raise NotImplementedError

    def _create_doc(self):
        doc = getattr(np, self.__name__).__doc__
        doc = _remove_unused_reference_labels(
            _skip_signature(_dedent(doc), self.__name__)
        )
        self.__doc__ = (
            f"xarray specific variant of :py:func:`numpy.{self.__name__}`. "
            "Handles xarray objects by dispatching to the appropriate "
            "function for the underlying array type.\n\n"
            f"Documentation from numpy:\n\n{doc}"
        )


class _unary_ufunc(_ufunc_wrapper):
    """Wrapper for dispatching unary ufuncs."""

    def __call__(self, x, /, **kwargs):
        xp = get_array_namespace(x)
        func = getattr(xp, self.__name__)
        return xr.apply_ufunc(func, x, dask="allowed", **kwargs)


class _binary_ufunc(_ufunc_wrapper):
    """Wrapper for dispatching binary ufuncs."""

    def __call__(self, x, y, /, **kwargs):
        xp = get_array_namespace(x, y)
        func = getattr(xp, self.__name__)
        return xr.apply_ufunc(func, x, y, dask="allowed", **kwargs)


def _skip_signature(doc, name):
    if not isinstance(doc, str):
        return doc

    # numpy creates some functions as aliases and copies the docstring exactly,
    # so check the actual name to handle this case
    np_name = getattr(np, name).__name__
    if doc.startswith(np_name):
        signature_end = doc.find("\n\n")
        doc = doc[signature_end + 2 :]

    return doc


def _remove_unused_reference_labels(doc):
    if not isinstance(doc, str):
        return doc

    max_references = 5
    for num in range(max_references):
        label = f".. [{num}]"
        reference = f"[{num}]_"
        index = f"{num}.    "

        if label not in doc or reference in doc:
            continue

        doc = doc.replace(label, index)

    return doc


def _dedent(doc):
    if not isinstance(doc, str):
        return doc

    return textwrap.dedent(doc)


# These can be auto-generated from the public numpy ufuncs:
# {name for name in dir(np) if isinstance(getattr(np, name), np.ufunc)}

# Generalized ufuncs that use core dimensions or produce multiple output
# arrays are not currently supported, and left commented out below.

# UNARY
abs = _unary_ufunc("abs")
absolute = _unary_ufunc("absolute")
acos = _unary_ufunc("acos")
acosh = _unary_ufunc("acosh")
arccos = _unary_ufunc("arccos")
arccosh = _unary_ufunc("arccosh")
arcsin = _unary_ufunc("arcsin")
arcsinh = _unary_ufunc("arcsinh")
arctan = _unary_ufunc("arctan")
arctanh = _unary_ufunc("arctanh")
asin = _unary_ufunc("asin")
asinh = _unary_ufunc("asinh")
atan = _unary_ufunc("atan")
atanh = _unary_ufunc("atanh")
bitwise_count = _unary_ufunc("bitwise_count")
bitwise_invert = _unary_ufunc("bitwise_invert")
bitwise_not = _unary_ufunc("bitwise_not")
cbrt = _unary_ufunc("cbrt")
ceil = _unary_ufunc("ceil")
conj = _unary_ufunc("conj")
conjugate = _unary_ufunc("conjugate")
cos = _unary_ufunc("cos")
cosh = _unary_ufunc("cosh")
deg2rad = _unary_ufunc("deg2rad")
degrees = _unary_ufunc("degrees")
exp = _unary_ufunc("exp")
exp2 = _unary_ufunc("exp2")
expm1 = _unary_ufunc("expm1")
fabs = _unary_ufunc("fabs")
floor = _unary_ufunc("floor")
# frexp = _unary_ufunc("frexp")
invert = _unary_ufunc("invert")
isfinite = _unary_ufunc("isfinite")
isinf = _unary_ufunc("isinf")
isnan = _unary_ufunc("isnan")
isnat = _unary_ufunc("isnat")
log = _unary_ufunc("log")
log10 = _unary_ufunc("log10")
log1p = _unary_ufunc("log1p")
log2 = _unary_ufunc("log2")
logical_not = _unary_ufunc("logical_not")
# modf = _unary_ufunc("modf")
negative = _unary_ufunc("negative")
positive = _unary_ufunc("positive")
rad2deg = _unary_ufunc("rad2deg")
radians = _unary_ufunc("radians")
reciprocal = _unary_ufunc("reciprocal")
rint = _unary_ufunc("rint")
sign = _unary_ufunc("sign")
signbit = _unary_ufunc("signbit")
sin = _unary_ufunc("sin")
sinh = _unary_ufunc("sinh")
spacing = _unary_ufunc("spacing")
sqrt = _unary_ufunc("sqrt")
square = _unary_ufunc("square")
tan = _unary_ufunc("tan")
tanh = _unary_ufunc("tanh")
trunc = _unary_ufunc("trunc")

# BINARY
add = _binary_ufunc("add")
arctan2 = _binary_ufunc("arctan2")
atan2 = _binary_ufunc("atan2")
bitwise_and = _binary_ufunc("bitwise_and")
bitwise_left_shift = _binary_ufunc("bitwise_left_shift")
bitwise_or = _binary_ufunc("bitwise_or")
bitwise_right_shift = _binary_ufunc("bitwise_right_shift")
bitwise_xor = _binary_ufunc("bitwise_xor")
copysign = _binary_ufunc("copysign")
divide = _binary_ufunc("divide")
# divmod = _binary_ufunc("divmod")
equal = _binary_ufunc("equal")
float_power = _binary_ufunc("float_power")
floor_divide = _binary_ufunc("floor_divide")
fmax = _binary_ufunc("fmax")
fmin = _binary_ufunc("fmin")
fmod = _binary_ufunc("fmod")
gcd = _binary_ufunc("gcd")
greater = _binary_ufunc("greater")
greater_equal = _binary_ufunc("greater_equal")
heaviside = _binary_ufunc("heaviside")
hypot = _binary_ufunc("hypot")
lcm = _binary_ufunc("lcm")
ldexp = _binary_ufunc("ldexp")
left_shift = _binary_ufunc("left_shift")
less = _binary_ufunc("less")
less_equal = _binary_ufunc("less_equal")
logaddexp = _binary_ufunc("logaddexp")
logaddexp2 = _binary_ufunc("logaddexp2")
logical_and = _binary_ufunc("logical_and")
logical_or = _binary_ufunc("logical_or")
logical_xor = _binary_ufunc("logical_xor")
# matmul = _binary_ufunc("matmul")
maximum = _binary_ufunc("maximum")
minimum = _binary_ufunc("minimum")
mod = _binary_ufunc("mod")
multiply = _binary_ufunc("multiply")
nextafter = _binary_ufunc("nextafter")
not_equal = _binary_ufunc("not_equal")
pow = _binary_ufunc("pow")
power = _binary_ufunc("power")
remainder = _binary_ufunc("remainder")
right_shift = _binary_ufunc("right_shift")
subtract = _binary_ufunc("subtract")
true_divide = _binary_ufunc("true_divide")
# vecdot = _binary_ufunc("vecdot")

# elementwise non-ufunc
angle = _unary_ufunc("angle")
isreal = _unary_ufunc("isreal")
iscomplex = _unary_ufunc("iscomplex")


__all__ = [
    "abs",
    "absolute",
    "acos",
    "acosh",
    "add",
    "angle",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctan2",
    "arctanh",
    "asin",
    "asinh",
    "atan",
    "atan2",
    "atanh",
    "bitwise_and",
    "bitwise_count",
    "bitwise_invert",
    "bitwise_left_shift",
    "bitwise_not",
    "bitwise_or",
    "bitwise_right_shift",
    "bitwise_xor",
    "cbrt",
    "ceil",
    "conj",
    "conjugate",
    "copysign",
    "cos",
    "cosh",
    "deg2rad",
    "degrees",
    "divide",
    "equal",
    "exp",
    "exp2",
    "expm1",
    "fabs",
    "float_power",
    "floor",
    "floor_divide",
    "fmax",
    "fmin",
    "fmod",
    "gcd",
    "greater",
    "greater_equal",
    "heaviside",
    "hypot",
    "invert",
    "iscomplex",
    "isfinite",
    "isinf",
    "isnan",
    "isnat",
    "isreal",
    "lcm",
    "ldexp",
    "left_shift",
    "less",
    "less_equal",
    "log",
    "log1p",
    "log2",
    "log10",
    "logaddexp",
    "logaddexp2",
    "logical_and",
    "logical_not",
    "logical_or",
    "logical_xor",
    "maximum",
    "minimum",
    "mod",
    "multiply",
    "negative",
    "nextafter",
    "not_equal",
    "positive",
    "pow",
    "power",
    "rad2deg",
    "radians",
    "reciprocal",
    "remainder",
    "right_shift",
    "rint",
    "sign",
    "signbit",
    "sin",
    "sinh",
    "spacing",
    "sqrt",
    "square",
    "subtract",
    "tan",
    "tanh",
    "true_divide",
    "trunc",
]
