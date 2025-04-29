from llvmlite.binding import ffi

# FIXME: Remove me once typed pointers are no longer supported.
from llvmlite import opaque_pointers_enabled
from ctypes import c_bool


def create_context():
    return ContextRef(
        ffi.lib.LLVMPY_ContextCreate(opaque_pointers_enabled))


def get_global_context():
    return GlobalContextRef(
        ffi.lib.LLVMPY_GetGlobalContext(opaque_pointers_enabled))


class ContextRef(ffi.ObjectRef):
    def __init__(self, context_ptr):
        super(ContextRef, self).__init__(context_ptr)

    def _dispose(self):
        ffi.lib.LLVMPY_ContextDispose(self)


class GlobalContextRef(ContextRef):
    def _dispose(self):
        pass


# FIXME: Remove argtypes once typed pointers are no longer supported.
ffi.lib.LLVMPY_GetGlobalContext.argtypes = [c_bool]
ffi.lib.LLVMPY_GetGlobalContext.restype = ffi.LLVMContextRef

# FIXME: Remove argtypes once typed pointers are no longer supported.
ffi.lib.LLVMPY_ContextCreate.argtypes = [c_bool]
ffi.lib.LLVMPY_ContextCreate.restype = ffi.LLVMContextRef

ffi.lib.LLVMPY_ContextDispose.argtypes = [ffi.LLVMContextRef]
