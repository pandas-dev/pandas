# -*- coding: utf-8 -*-
"""
Support for CFFI. Allows checking whether objects are CFFI functions and
obtaining the pointer and numba signature.
"""

from types import BuiltinFunctionType
import ctypes
from functools import partial
import numpy as np

from numba.core import types
from numba.core.errors import TypingError
from numba.core.typing import templates
from numba.np import numpy_support

try:
    import cffi
    ffi = cffi.FFI()
except ImportError:
    ffi = None

SUPPORTED = ffi is not None
_ool_func_types = {}
_ool_func_ptr = {}
_ffi_instances = set()


def is_ffi_instance(obj):
    # Compiled FFI modules have a member, ffi, which is an instance of
    # CompiledFFI, which behaves similarly to an instance of cffi.FFI. In
    # order to simplify handling a CompiledFFI object, we treat them as
    # if they're cffi.FFI instances for typing and lowering purposes.
    try:
        return obj in _ffi_instances or isinstance(obj, cffi.FFI)
    except TypeError: # Unhashable type possible
        return False

def is_cffi_func(obj):
    """Check whether the obj is a CFFI function"""
    try:
        return ffi.typeof(obj).kind == 'function'
    except TypeError:
        try:
            return obj in _ool_func_types
        except:
            return False

def get_pointer(cffi_func):
    """
    Get a pointer to the underlying function for a CFFI function as an
    integer.
    """
    if cffi_func in _ool_func_ptr:
        return _ool_func_ptr[cffi_func]
    return int(ffi.cast("uintptr_t", cffi_func))


_cached_type_map = None

def _type_map():
    """
    Lazily compute type map, as calling ffi.typeof() involves costly
    parsing of C code...
    """
    global _cached_type_map
    if _cached_type_map is None:
        _cached_type_map = {
            ffi.typeof('bool') :                types.boolean,
            ffi.typeof('char') :                types.char,
            ffi.typeof('short') :               types.short,
            ffi.typeof('int') :                 types.intc,
            ffi.typeof('long') :                types.long_,
            ffi.typeof('long long') :           types.longlong,
            ffi.typeof('unsigned char') :       types.uchar,
            ffi.typeof('unsigned short') :      types.ushort,
            ffi.typeof('unsigned int') :        types.uintc,
            ffi.typeof('unsigned long') :       types.ulong,
            ffi.typeof('unsigned long long') :  types.ulonglong,
            ffi.typeof('int8_t') :              types.char,
            ffi.typeof('uint8_t') :             types.uchar,
            ffi.typeof('int16_t') :             types.short,
            ffi.typeof('uint16_t') :            types.ushort,
            ffi.typeof('int32_t') :             types.intc,
            ffi.typeof('uint32_t') :            types.uintc,
            ffi.typeof('int64_t') :             types.longlong,
            ffi.typeof('uint64_t') :            types.ulonglong,
            ffi.typeof('float') :               types.float32,
            ffi.typeof('double') :              types.double,
            ffi.typeof('ssize_t') :             types.intp,
            ffi.typeof('size_t') :              types.uintp,
            ffi.typeof('void') :                types.void,
        }
    return _cached_type_map


def map_type(cffi_type, use_record_dtype=False):
    """
    Map CFFI type to numba type.

    Parameters
    ----------
    cffi_type:
        The CFFI type to be converted.
    use_record_dtype: bool (default: False)
        When True, struct types are mapped to a NumPy Record dtype.

    """
    primed_map_type = partial(map_type, use_record_dtype=use_record_dtype)
    kind = getattr(cffi_type, 'kind', '')
    if kind == 'union':
        raise TypeError("No support for CFFI union")
    elif kind == 'function':
        if cffi_type.ellipsis:
            raise TypeError("vararg function is not supported")
        restype = primed_map_type(cffi_type.result)
        argtypes = [primed_map_type(arg) for arg in cffi_type.args]
        return templates.signature(restype, *argtypes)
    elif kind == 'pointer':
        pointee = cffi_type.item
        if pointee.kind == 'void':
            return types.voidptr
        else:
            return types.CPointer(primed_map_type(pointee))
    elif kind == 'array':
        dtype = primed_map_type(cffi_type.item)
        nelem = cffi_type.length
        return types.NestedArray(dtype=dtype, shape=(nelem,))
    elif kind == 'struct' and use_record_dtype:
        return map_struct_to_record_dtype(cffi_type)
    else:
        result = _type_map().get(cffi_type)
        if result is None:
            raise TypeError(cffi_type)
        return result


def map_struct_to_record_dtype(cffi_type):
    """Convert a cffi type into a NumPy Record dtype
    """
    fields = {
            'names': [],
            'formats': [],
            'offsets': [],
            'itemsize': ffi.sizeof(cffi_type),
    }
    is_aligned = True
    for k, v in cffi_type.fields:
        # guard unsupported values
        if v.bitshift != -1:
            msg = "field {!r} has bitshift, this is not supported"
            raise ValueError(msg.format(k))
        if v.flags != 0:
            msg = "field {!r} has flags, this is not supported"
            raise ValueError(msg.format(k))
        if v.bitsize != -1:
            msg = "field {!r} has bitsize, this is not supported"
            raise ValueError(msg.format(k))
        dtype = numpy_support.as_dtype(
            map_type(v.type, use_record_dtype=True),
        )
        fields['names'].append(k)
        fields['formats'].append(dtype)
        fields['offsets'].append(v.offset)
        # Check alignment
        is_aligned &= (v.offset % dtype.alignment == 0)

    return numpy_support.from_dtype(np.dtype(fields, align=is_aligned))


def make_function_type(cffi_func, use_record_dtype=False):
    """
    Return a Numba type for the given CFFI function pointer.
    """
    cffi_type = _ool_func_types.get(cffi_func) or ffi.typeof(cffi_func)
    if getattr(cffi_type, 'kind', '') == 'struct':
        raise TypeError('No support for CFFI struct values')
    sig = map_type(cffi_type, use_record_dtype=use_record_dtype)
    return types.ExternalFunctionPointer(sig, get_pointer=get_pointer)


registry = templates.Registry()

@registry.register
class FFI_from_buffer(templates.AbstractTemplate):
    key = 'ffi.from_buffer'

    def generic(self, args, kws):
        if kws or len(args) != 1:
            return
        [ary] = args
        if not isinstance(ary, types.Buffer):
            raise TypingError("from_buffer() expected a buffer object, got %s"
                              % (ary,))
        if ary.layout not in ('C', 'F'):
            raise TypingError("from_buffer() unsupported on non-contiguous buffers (got %s)"
                              % (ary,))
        if ary.layout != 'C' and ary.ndim > 1:
            raise TypingError("from_buffer() only supports multidimensional arrays with C layout (got %s)"
                              % (ary,))
        ptr = types.CPointer(ary.dtype)
        return templates.signature(ptr, ary)

@registry.register_attr
class FFIAttribute(templates.AttributeTemplate):
    key = types.ffi

    def resolve_from_buffer(self, ffi):
        return types.BoundFunction(FFI_from_buffer, types.ffi)


def register_module(mod):
    """
    Add typing for all functions in an out-of-line CFFI module to the typemap
    """
    for f in dir(mod.lib):
        f = getattr(mod.lib, f)
        if isinstance(f, BuiltinFunctionType):
            _ool_func_types[f] = mod.ffi.typeof(f)
            addr = mod.ffi.addressof(mod.lib, f.__name__)
            _ool_func_ptr[f] = int(mod.ffi.cast("uintptr_t", addr))
        _ffi_instances.add(mod.ffi)

def register_type(cffi_type, numba_type):
    """
    Add typing for a given CFFI type to the typemap
    """
    tm = _type_map()
    tm[cffi_type] = numba_type
