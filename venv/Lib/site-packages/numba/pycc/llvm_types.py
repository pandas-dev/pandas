import sys
import ctypes
import struct as struct_
import llvmlite.ir
from llvmlite.ir import Constant

_trace_refs_ = hasattr(sys, 'getobjects')
_plat_bits = struct_.calcsize('@P') * 8

_int8 = llvmlite.ir.IntType(8)
_int32 = llvmlite.ir.IntType(32)

_void_star = llvmlite.ir.PointerType(_int8)

_int8_star = _void_star

_sizeof_py_ssize_t = ctypes.sizeof(getattr(ctypes, 'c_size_t'))
_llvm_py_ssize_t = llvmlite.ir.IntType(_sizeof_py_ssize_t * 8)

if _trace_refs_:
    _pyobject_head = llvmlite.ir.LiteralStructType([_void_star, _void_star,
                                                    _llvm_py_ssize_t, _void_star])
    _pyobject_head_init = Constant.literal_struct([
        Constant(_void_star, None),            # _ob_next
        Constant(_void_star, None),            # _ob_prev
        Constant(_llvm_py_ssize_t, 1),         # ob_refcnt
        Constant(_void_star, None),            # ob_type
        ])

else:
    _pyobject_head = llvmlite.ir.LiteralStructType([_llvm_py_ssize_t, _void_star])
    _pyobject_head_init = Constant.literal_struct([
        Constant(_llvm_py_ssize_t, 1),    # ob_refcnt
        Constant(_void_star, None),       # ob_type
        ])

_pyobject_head_p = llvmlite.ir.PointerType(_pyobject_head)
