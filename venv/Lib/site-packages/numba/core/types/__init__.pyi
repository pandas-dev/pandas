# This file is provided by @jorenham with minor modifications (e.g. typos).
# See original content at: https://github.com/numba/numba/pull/9945#pullrequestreview-2668923222.
#
# This file has been tested under:
#   - mypy for the use-case in issue #9900
#   - mypy numba/core/types/__init__.pyi
# Testing with mypy.stubtest does not work due to other mypy errors in the code
# base.
from .abstract import *
from .common import Opaque
from .containers import *
from .function_type import *
from .functions import *
from .iterators import *
from .misc import *
from .npytypes import *
from .scalars import *

__all__ = [
    "b1",
    "bool",  # numpy>=2
    "bool_",
    "boolean",
    "byte",
    "c8",
    "c16",
    "char",
    "complex64",
    "complex128",
    "deferred_type",
    "double",
    "f4",
    "f8",
    "ffi",
    "ffi_forced_object",
    "float32",
    "float64",
    "float_",  # numpy<2
    "i1",
    "i2",
    "i4",
    "i8",
    "int8",
    "int16",
    "int32",
    "int64",
    "int_",
    "intc",
    "intp",
    "long_",
    "longlong",
    "none",
    "optional",
    "short",
    "size_t",
    "ssize_t",
    "u1",
    "u2",
    "u4",
    "u8",
    "uchar",
    "uint",
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    "uintc",
    "uintp",
    "ulong",
    "ulonglong",
    "ushort",
    "void",
]

# TODO: Final

pyobject: PyObject = ...
ffi_forced_object: Opaque = ...
ffi: Opaque = ...
none: NoneType = ...
ellipsis: EllipsisType = ...
Any: Phantom = ...
undefined: Undefined = ...
py2_string_type: Opaque = ...
unicode_type: UnicodeType = ...
string: UnicodeType = ...
unknown: Dummy = ...
npy_rng: NumPyRandomGeneratorType = ...
npy_bitgen: NumPyRandomBitGeneratorType = ...

_undef_var: UndefVar = ...

code_type: Opaque = ...
pyfunc_type: Opaque = ...

voidptr: RawPointer = ...

optional = Optional
deferred_type = DeferredType
slice2_type: SliceType = ...
slice3_type: SliceType = ...
void: NoneType = ...

boolean: Boolean = ...
bool_: Boolean = boolean
bool: Boolean = boolean  # numpy>=2

int8: Integer = ...
int16: Integer = ...
int32: Integer = ...
int64: Integer = ...
intp: Integer = ...
intc: Integer = ...
ssize_t: Integer = ...
char: Integer = ...
short: Integer = ...
int_: Integer = ...
long_: Integer = ...
longlong: Integer = ...

byte: Integer = ...
uint8: Integer = ...
uint16: Integer = ...
uint32: Integer = ...
uint64: Integer = ...
uintp: Integer = ...
uintc: Integer = ...
size_t: Integer = ...
uchar: Integer = ...
ushort: Integer = ...
uint: Integer = ...
ulong: Integer = ...
ulonglong: Integer = ...

float16: Float = ...
float32: Float = ...
float64: Float = ...
float_: Float = float32  # numpy<2
double: Float = float64

# TODO: make generic in the wrapped `Float` type
complex64: Complex = ...
complex128: Complex = ...

range_iter32_type: RangeIteratorType = ...
range_iter64_type: RangeIteratorType = ...
unsigned_range_iter64_type: RangeIteratorType = ...
range_state32_type: RangeType = ...
range_state64_type: RangeType = ...
unsigned_range_state64_type: RangeType = ...

signed_domain: frozenset[Integer] = ...
unsigned_domain: frozenset[Integer] = ...
integer_domain: frozenset[Integer] = ...
real_domain: frozenset[Float] = ...
complex_domain: frozenset[Complex] = ...
number_domain: frozenset[Integer | Float | Complex] = ...

c_bool: Boolean = boolean
c_int8: Integer = int8
c_int16: Integer = int16
c_int32: Integer = int32
c_int64: Integer = int64
c_intp: Integer = intp
c_uint8: Integer = uint8
c_uint16: Integer = uint16
c_uint32: Integer = uint32
c_uint64: Integer = uint64
c_uintp: Integer = uintp
c_float16: Float = float16
c_float32: Float = float32
c_float64: Float = float64

np_bool_: Boolean = boolean
np_int8: Integer = int8
np_int16: Integer = int16
np_int32: Integer = int32
np_int64: Integer = int64
np_intp: Integer = intp
np_uint8: Integer = uint8
np_uint16: Integer = uint16
np_uint32: Integer = uint32
np_uint64: Integer = uint64
np_uintp: Integer = uintp
np_float16: Float = float16
np_float32: Float = float32
np_float64: Float = float64
np_float_: Float = float32
np_double: Float = float64
np_complex64: Complex = complex64
np_complex128: Complex = complex128

py_bool: Boolean = boolean
py_int: Integer = intp
py_float: Float = float64
py_complex: Complex = complex128

py_signed_domain: frozenset[Integer] = signed_domain
py_integer_domain: frozenset[Integer] = integer_domain
py_real_domain: frozenset[Float] = real_domain
py_complex_domain: frozenset[Complex] = complex_domain
py_number_domain: frozenset[Integer | Float | Complex] = number_domain

np_signed_domain: frozenset[Integer] = signed_domain
np_unsigned_domain: frozenset[Integer] = unsigned_domain
np_integer_domain: frozenset[Integer] = integer_domain
np_real_domain: frozenset[Float] = real_domain
np_complex_domain: frozenset[Complex] = complex_domain
np_number_domain: frozenset[Integer | Float | Complex] = number_domain

b1: Boolean = bool_
i1: Integer = int8
i2: Integer = int16
i4: Integer = int32
i8: Integer = int64
u1: Integer = uint8
u2: Integer = uint16
u4: Integer = uint32
u8: Integer = uint64
f2: Float = float16
f4: Float = float32
f8: Float = float64
c8: Complex = complex64
c16: Complex = complex128
