from llvmlite.binding import ffi


def create_context():
    return ContextRef(ffi.lib.LLVMPY_ContextCreate())


def get_global_context():
    return GlobalContextRef(ffi.lib.LLVMPY_GetGlobalContext())


class ContextRef(ffi.ObjectRef):
    def __init__(self, context_ptr):
        super(ContextRef, self).__init__(context_ptr)

    def _dispose(self):
        ffi.lib.LLVMPY_ContextDispose(self)


class GlobalContextRef(ContextRef):
    def _dispose(self):
        pass


ffi.lib.LLVMPY_GetGlobalContext.restype = ffi.LLVMContextRef

ffi.lib.LLVMPY_ContextCreate.restype = ffi.LLVMContextRef

ffi.lib.LLVMPY_ContextDispose.argtypes = [ffi.LLVMContextRef]
