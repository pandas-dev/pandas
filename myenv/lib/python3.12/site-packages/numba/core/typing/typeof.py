from collections import namedtuple
from functools import singledispatch
import ctypes
import enum

import numpy as np
from numpy.random.bit_generator import BitGenerator

from numba.core import types, utils, errors
from numba.np import numpy_support


# terminal color markup
_termcolor = errors.termcolor()


class Purpose(enum.Enum):
    # Value being typed is used as an argument
    argument = 1
    # Value being typed is used as a constant
    constant = 2


_TypeofContext = namedtuple("_TypeofContext", ("purpose",))


def typeof(val, purpose=Purpose.argument):
    """
    Get the Numba type of a Python value for the given purpose.
    """
    # Note the behaviour for Purpose.argument must match _typeof.c.
    c = _TypeofContext(purpose)
    ty = typeof_impl(val, c)
    if ty is None:
        msg = _termcolor.errmsg(
            f"Cannot determine Numba type of {type(val)}")
        raise ValueError(msg)
    return ty


@singledispatch
def typeof_impl(val, c):
    """
    Generic typeof() implementation.
    """
    tp = _typeof_buffer(val, c)
    if tp is not None:
        return tp

    tp = getattr(val, "_numba_type_", None)
    if tp is not None:
        return tp

    # cffi is handled here as it does not expose a public base class
    # for exported functions or CompiledFFI instances.
    from numba.core.typing import cffi_utils
    if cffi_utils.SUPPORTED:
        if cffi_utils.is_cffi_func(val):
            return cffi_utils.make_function_type(val)
        if cffi_utils.is_ffi_instance(val):
            return types.ffi

    return None


def _typeof_buffer(val, c):
    from numba.core.typing import bufproto
    try:
        m = memoryview(val)
    except TypeError:
        return
    # Object has the buffer protocol
    try:
        dtype = bufproto.decode_pep3118_format(m.format, m.itemsize)
    except ValueError:
        return
    type_class = bufproto.get_type_class(type(val))
    layout = bufproto.infer_layout(m)
    return type_class(dtype, m.ndim, layout=layout,
                      readonly=m.readonly)


@typeof_impl.register(ctypes._CFuncPtr)
def _typeof_ctypes_function(val, c):
    from .ctypes_utils import is_ctypes_funcptr, make_function_type
    if is_ctypes_funcptr(val):
        return make_function_type(val)


@typeof_impl.register(type)
def _typeof_type(val, c):
    """
    Type various specific Python types.
    """
    if issubclass(val, BaseException):
        return types.ExceptionClass(val)
    if issubclass(val, tuple) and hasattr(val, "_asdict"):
        return types.NamedTupleClass(val)

    if issubclass(val, np.generic):
        return types.NumberClass(numpy_support.from_dtype(val))

    if issubclass(val, types.Type):
        return types.TypeRef(val)

    from numba.typed import Dict
    if issubclass(val, Dict):
        return types.TypeRef(types.DictType)

    from numba.typed import List
    if issubclass(val, List):
        return types.TypeRef(types.ListType)


@typeof_impl.register(bool)
def _typeof_bool(val, c):
    return types.boolean


@typeof_impl.register(float)
def _typeof_float(val, c):
    return types.float64


@typeof_impl.register(complex)
def _typeof_complex(val, c):
    return types.complex128


@typeof_impl.register(int)
def _typeof_int(val, c):
    # As in _typeof.c
    nbits = utils.bit_length(val)
    if nbits < 32:
        typ = types.intp
    elif nbits < 64:
        typ = types.int64
    elif nbits == 64 and val >= 0:
        typ = types.uint64
    else:
        raise ValueError("Int value is too large: %s" % val)
    return typ


@typeof_impl.register(np.generic)
def _typeof_numpy_scalar(val, c):
    try:
        return numpy_support.map_arrayscalar_type(val)
    except NotImplementedError:
        pass


@typeof_impl.register(str)
def _typeof_str(val, c):
    return types.string


@typeof_impl.register(type((lambda a: a).__code__))
def _typeof_code(val, c):
    return types.code_type


@typeof_impl.register(type(None))
def _typeof_none(val, c):
    return types.none


@typeof_impl.register(type(Ellipsis))
def _typeof_ellipsis(val, c):
    return types.ellipsis


@typeof_impl.register(tuple)
def _typeof_tuple(val, c):
    tys = [typeof_impl(v, c) for v in val]
    if any(ty is None for ty in tys):
        return
    return types.BaseTuple.from_types(tys, type(val))


@typeof_impl.register(list)
def _typeof_list(val, c):
    if len(val) == 0:
        raise ValueError("Cannot type empty list")
    ty = typeof_impl(val[0], c)
    if ty is None:
        raise ValueError(
            f"Cannot type list element type {type(val[0])}")
    return types.List(ty, reflected=True)


@typeof_impl.register(set)
def _typeof_set(val, c):
    if len(val) == 0:
        raise ValueError("Cannot type empty set")
    item = next(iter(val))
    ty = typeof_impl(item, c)
    if ty is None:
        raise ValueError(
            f"Cannot type set element type {type(item)}")
    return types.Set(ty, reflected=True)


@typeof_impl.register(slice)
def _typeof_slice(val, c):
    return types.slice2_type if val.step in (None, 1) else types.slice3_type


@typeof_impl.register(enum.Enum)
@typeof_impl.register(enum.IntEnum)
def _typeof_enum(val, c):
    clsty = typeof_impl(type(val), c)
    return clsty.member_type


@typeof_impl.register(enum.EnumMeta)
def _typeof_enum_class(val, c):
    cls = val
    members = list(cls.__members__.values())
    if len(members) == 0:
        raise ValueError("Cannot type enum with no members")
    dtypes = {typeof_impl(mem.value, c) for mem in members}
    if len(dtypes) > 1:
        raise ValueError("Cannot type heterogeneous enum: "
                         "got value types %s"
                         % ", ".join(sorted(str(ty) for ty in dtypes)))
    if issubclass(val, enum.IntEnum):
        typecls = types.IntEnumClass
    else:
        typecls = types.EnumClass
    return typecls(cls, dtypes.pop())


@typeof_impl.register(np.dtype)
def _typeof_dtype(val, c):
    tp = numpy_support.from_dtype(val)
    return types.DType(tp)


@typeof_impl.register(np.ndarray)
def _typeof_ndarray(val, c):
    if isinstance(val, np.ma.MaskedArray):
        msg = "Unsupported array type: numpy.ma.MaskedArray."
        raise errors.NumbaTypeError(msg)
    try:
        dtype = numpy_support.from_dtype(val.dtype)
    except errors.NumbaNotImplementedError:
        raise errors.NumbaValueError(f"Unsupported array dtype: {val.dtype}")
    layout = numpy_support.map_layout(val)
    readonly = not val.flags.writeable
    return types.Array(dtype, val.ndim, layout, readonly=readonly)


@typeof_impl.register(types.NumberClass)
def _typeof_number_class(val, c):
    return val


@typeof_impl.register(types.Literal)
def _typeof_literal(val, c):
    return val


@typeof_impl.register(types.TypeRef)
def _typeof_typeref(val, c):
    return val


@typeof_impl.register(types.Type)
def _typeof_nb_type(val, c):
    if isinstance(val, types.BaseFunction):
        return val
    elif isinstance(val, (types.Number, types.Boolean)):
        return types.NumberClass(val)
    else:
        return types.TypeRef(val)


@typeof_impl.register(BitGenerator)
def typeof_numpy_random_bitgen(val, c):
    return types.NumPyRandomBitGeneratorType(val)


@typeof_impl.register(np.random.Generator)
def typeof_random_generator(val, c):
    return types.NumPyRandomGeneratorType(val)


@typeof_impl.register(np.polynomial.polynomial.Polynomial)
def typeof_numpy_polynomial(val, c):
    coef = typeof(val.coef)
    domain = typeof(val.domain)
    window = typeof(val.window)
    return types.PolynomialType(coef, domain, window)
