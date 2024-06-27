'''
NVVM is not supported in the simulator, but stubs are provided to allow tests
to import correctly.
'''


class NvvmSupportError(ImportError):
    pass


class NVVM(object):
    def __init__(self):
        raise NvvmSupportError('NVVM not supported in the simulator')


CompilationUnit = None
compile_ir = None
set_cuda_kernel = None
get_arch_option = None
LibDevice = None
NvvmError = None


def is_available():
    return False


def get_supported_ccs():
    return ()
