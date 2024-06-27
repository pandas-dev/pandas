"""
Declarations of the Runtime API functions.
"""

from ctypes import c_int, POINTER

API_PROTOTYPES = {
    # cudaError_t cudaRuntimeGetVersion ( int* runtimeVersion )
    'cudaRuntimeGetVersion': (c_int, POINTER(c_int)),
}
