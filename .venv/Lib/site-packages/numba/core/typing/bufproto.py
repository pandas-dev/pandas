"""
Typing support for the buffer protocol (PEP 3118).
"""

import array

from numba.core import types, config
from numba.core.errors import NumbaValueError


_pep3118_int_types = set('bBhHiIlLqQnN')

if config.USE_LEGACY_TYPE_SYSTEM: # Old type system
    _pep3118_scalar_map = {
        'f': types.float32,
        'd': types.float64,
        'Zf': types.complex64,
        'Zd': types.complex128,
        }
else: # New type system
    _pep3118_scalar_map = {
        # TODO: FIXME We need to modify the following Map to use Python Types.
        # However currently here's nothing in Python types that maps
        # to a float32 or a complex64
        # 'f': types.np_float32,
        'd': types.py_float, # 64-bit float
        # 'Zf': types.np_complex64,
        'Zd': types.py_complex, # 128-bit complex
        }

_type_map = {
    bytearray: types.ByteArray,
    array.array: types.PyArray,
    }

_type_map[memoryview] = types.MemoryView
_type_map[bytes] = types.Bytes


def decode_pep3118_format(fmt, itemsize):
    """
    Return the Numba type for an item with format string *fmt* and size
    *itemsize* (in bytes).
    """
    # XXX reuse _dtype_from_pep3118() from np.core._internal?
    if fmt in _pep3118_int_types:
        # Determine int width and signedness
        name = 'int%d' % (itemsize * 8,)
        if fmt.isupper():
            name = 'u' + name
        return types.Integer(name)
    try:
        # For the hard-coded types above, consider "=" the same as "@"
        # (the default).  This is because Numpy sometimes adds "="
        # in front of the PEP 3118 format string.
        return _pep3118_scalar_map[fmt.lstrip('=')]
    except KeyError:
        raise NumbaValueError("unsupported PEP 3118 format %r" % (fmt,))


def get_type_class(typ):
    """
    Get the Numba type class for buffer-compatible Python *typ*.
    """
    try:
        # Look up special case.
        return _type_map[typ]
    except KeyError:
        # Fall back on generic one.
        return types.Buffer


def infer_layout(val):
    """
    Infer layout of the given memoryview *val*.
    """
    return ('C' if val.c_contiguous else
            'F' if val.f_contiguous else
            'A')
