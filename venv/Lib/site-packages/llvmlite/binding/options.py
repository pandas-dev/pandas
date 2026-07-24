from llvmlite.binding import ffi
from llvmlite.binding.common import _encode_string
from ctypes import c_char_p


def set_option(name, option):
    """
    Set the given LLVM "command-line" option.

    For example set_option("test", "-debug-pass=Structure") would display
    all optimization passes when generating code.
    """
    ffi.lib.LLVMPY_SetCommandLine(_encode_string(name),
                                  _encode_string(option))


ffi.lib.LLVMPY_SetCommandLine.argtypes = [c_char_p, c_char_p]
