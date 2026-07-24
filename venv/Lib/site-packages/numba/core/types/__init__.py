import struct

import numpy as np
from numba.core import utils
import ctypes

from .abstract import *
from .containers import *
from .functions import *
from .iterators import *
from .misc import *
from .npytypes import *
from .scalars import *
from .function_type import *

numpy_version = tuple(map(int, np.__version__.split('.')[:2]))

# Short names

pyobject = PyObject('pyobject')
ffi_forced_object = Opaque('ffi_forced_object')
ffi = Opaque('ffi')
none = NoneType('none')
ellipsis = EllipsisType('...')
Any = Phantom('any')
undefined = Undefined('undefined')
py2_string_type = Opaque('str')
unicode_type = UnicodeType('unicode_type')
string = unicode_type
unknown = Dummy('unknown')
npy_rng = NumPyRandomGeneratorType('rng')
npy_bitgen = NumPyRandomBitGeneratorType('bitgen')

# _undef_var is used to represent undefined variables in the type system.
_undef_var = UndefVar('_undef_var')

code_type = Opaque('code')
pyfunc_type = Opaque('pyfunc')

# No operation is defined on voidptr
# Can only pass it around
voidptr = RawPointer('void*')

# optional types
optional = Optional
deferred_type = DeferredType
slice2_type = SliceType('slice<a:b>', 2)
slice3_type = SliceType('slice<a:b:c>', 3)
void = none

# Need to ignore mypy errors because mypy cannot unify types for both
# the type systems even if they're logically mutually exclusive.
# mypy: ignore-errors


boolean = bool_ = Boolean('bool')
if numpy_version >= (2, 0):
    bool = bool_

byte = uint8 = Integer('uint8')
uint16 = Integer('uint16')
uint32 = Integer('uint32')
uint64 = Integer('uint64')

int8 = Integer('int8')
int16 = Integer('int16')
int32 = Integer('int32')
int64 = Integer('int64')
intp = int32 if utils.MACHINE_BITS == 32 else int64
uintp = uint32 if utils.MACHINE_BITS == 32 else uint64
intc = int32 if struct.calcsize('i') == 4 else int64
uintc = uint32 if struct.calcsize('I') == 4 else uint64
ssize_t = int32 if struct.calcsize('n') == 4 else int64
size_t = uint32 if struct.calcsize('N') == 4 else uint64

float32 = Float('float32')
float64 = Float('float64')
float16 = Float('float16')

complex64 = Complex('complex64', float32)
complex128 = Complex('complex128', float64)

range_iter32_type = RangeIteratorType(int32)
range_iter64_type = RangeIteratorType(int64)
unsigned_range_iter64_type = RangeIteratorType(uint64)
range_state32_type = RangeType(int32)
range_state64_type = RangeType(int64)
unsigned_range_state64_type = RangeType(uint64)

signed_domain = frozenset([int8, int16, int32, int64])
unsigned_domain = frozenset([uint8, uint16, uint32, uint64])
integer_domain = signed_domain | unsigned_domain
real_domain = frozenset([float32, float64])
complex_domain = frozenset([complex64, complex128])
number_domain = real_domain | integer_domain | complex_domain

# Integer Aliases
c_bool = py_bool = np_bool_ = boolean

c_uint8 = np_uint8 = uint8
c_uint16 = np_uint16 = uint16
c_uint32 = np_uint32 = uint32
c_uint64 = np_uint64 = uint64
c_uintp = np_uintp = uintp

c_int8 = np_int8 = int8
c_int16 = np_int16 = int16
c_int32 = np_int32 = int32
c_int64 = np_int64 = int64
c_intp = py_int = np_intp = intp

c_float16 = np_float16 = float16
c_float32 = np_float32 = float32
c_float64 = py_float = np_float64 = float64

np_complex64 = complex64
py_complex = np_complex128 = complex128

# Domain Aliases
py_signed_domain = np_signed_domain = signed_domain
np_unsigned_domain = unsigned_domain
py_integer_domain = np_integer_domain = integer_domain
py_real_domain = np_real_domain = real_domain
py_complex_domain = np_complex_domain = complex_domain
py_number_domain = np_number_domain = number_domain

# Aliases to NumPy type names

b1 = bool_
i1 = int8
i2 = int16
i4 = int32
i8 = int64
u1 = uint8
u2 = uint16
u4 = uint32
u8 = uint64

f2 = float16
f4 = float32
f8 = float64

c8 = complex64
c16 = complex128

np_float_ = float32
np_double = double = float64
if numpy_version < (2, 0):
    float_ = float32

_make_signed = lambda x: globals()["int%d" % (ctypes.sizeof(x) * 8)]
_make_unsigned = lambda x: globals()["uint%d" % (ctypes.sizeof(x) * 8)]

char = np_char = _make_signed(ctypes.c_char)
uchar = np_uchar = byte = _make_unsigned(ctypes.c_ubyte)
short = np_short = _make_signed(ctypes.c_short)
ushort = np_ushort = _make_unsigned(ctypes.c_ushort)
int_ = np_int_ = np_intp
uint = np_uint = np_uintp
intc = np_intc = _make_signed(ctypes.c_int) # C-compat int
uintc = np_uintc = _make_unsigned(ctypes.c_uint) # C-compat uint
long_ = np_long = _make_signed(ctypes.c_long)  # C-compat long
ulong = np_ulong = _make_unsigned(ctypes.c_ulong)  # C-compat ulong
longlong = np_longlong = _make_signed(ctypes.c_longlong)
ulonglong = np_ulonglong = _make_unsigned(ctypes.c_ulonglong)

# This is equivalent to NumPy's `np.dtype('U1').itemsize`,
# which is the size of a single Unicode character in bytes.

# We can't keep this as `ctypes.c_wchar` because its size 
# is platform-dependent (2 bytes on Windows, 4 bytes on Unix).
sizeof_unicode_char = ctypes.sizeof(ctypes.c_byte) * 4

all_str = '''
int8
int16
int32
int64
uint8
uint16
uint32
uint64
intp
uintp
intc
uintc
ssize_t
size_t
boolean
float32
float64
complex64
complex128
bool_
byte
char
uchar
short
ushort
int_
uint
long_
ulong
longlong
ulonglong
float_
double
void
none
b1
i1
i2
i4
i8
u1
u2
u4
u8
f4
f8
c8
c16
optional
ffi_forced_object
ffi
deferred_type
'''

__all__ = all_str.split()
if numpy_version >= (2, 0):
    __all__.remove('float_')
    __all__.append('bool')

from numba.np.types.datetime import NPDatetime, NPTimedelta
