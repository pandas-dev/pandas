import os
import warnings
from functools import cache
from ctypes import c_int, c_char_p
from llvmlite.binding import ffi

# these are here as they cannot be lazy bound as the module globals make calls
# to the LLVMPY API functions
ffi.lib.LLVMPY_HasSVMLSupport.argtypes = ()
ffi.lib.LLVMPY_HasSVMLSupport.restype = c_int

ffi.lib.LLVMPY_IsStaticLibstdcxxLinkageBuild.argtypes = ()
ffi.lib.LLVMPY_IsStaticLibstdcxxLinkageBuild.restype = c_int

ffi.lib.LLVMPY_IsDynamicLLVMLinkageBuild.argtypes = ()
ffi.lib.LLVMPY_IsDynamicLLVMLinkageBuild.restype = c_int

ffi.lib.LLVMPY_PackageFormat.argtypes = ()
ffi.lib.LLVMPY_PackageFormat.restype = c_char_p

ffi.lib.LLVMPY_LlvmAssertionsState.argtypes = ()
ffi.lib.LLVMPY_LlvmAssertionsState.restype = c_char_p


def _has_svml():
    """
    Returns True if SVML was enabled at FFI support compile time.
    """
    if ffi.lib.LLVMPY_HasSVMLSupport() == 0:
        return False
    else:
        return True


has_svml = _has_svml()


def _build_llvm_linkage_type():
    """
    Returns "static" if the FFI support is statically linked against LLVM,
    returns "dynamic" otherwise.
    """
    if ffi.lib.LLVMPY_IsDynamicLLVMLinkageBuild() == 0:
        return "static"
    else:
        return "dynamic"


build_llvm_linkage_type = _build_llvm_linkage_type()


def _build_libstdcxx_linkage_type():
    """
    Returns "static" if the FFI support is statically linked against libstdc++,
    returns "dynamic" otherwise.
    """
    if ffi.lib.LLVMPY_IsStaticLibstdcxxLinkageBuild() == 1:
        return "static"
    else:
        return "dynamic"


build_libstdcxx_linkage_type = _build_libstdcxx_linkage_type()


def _package_format():
    """
    Returns "wheel", "conda" or "unspecified"
    """
    return ffi.lib.LLVMPY_PackageFormat().decode()


package_format = _package_format()


def _llvm_assertions_state():
    """
    Returns one of "on", "off" or "unknown". Depending on whether it is
    determined that LLVM was build with assertions on, off, or is not known.
    "Is not known" is typically from a dynamic linkage against LLVM in which
    case it's not easily identified whether LLVM was built with assertions.
    """
    return ffi.lib.LLVMPY_LlvmAssertionsState().decode()


llvm_assertions_state = _llvm_assertions_state()


@cache
def get_sysinfo():
    d = dict()
    d["ffi_lib_location"] = ffi.lib._name
    d["package_format"] = package_format
    d["llvm_linkage_type"] = build_llvm_linkage_type
    d["libstdcxx_linkage_type"] = build_libstdcxx_linkage_type
    d["llvm_assertions_state"] = llvm_assertions_state

    # import lief
    HAVE_LIEF = False
    try:
        import lief
        HAVE_LIEF = True
    except ImportError:
        msg = "py-lief package not found, sysinfo is limited as a result"
        warnings.warn(msg)

    d["lief_probe_status"] = HAVE_LIEF
    d["linked_libraries"] = None
    d["canonicalised_linked_libraries"] = None

    def canonicalise_library_type(dso):
        """Canonicalises the representation of the binary::libraries as a
        sequence of strings"""
        # Note lief v16:
        # Mach-O .libraries are DylibCommand instances.
        # Windows PE and Linux ELF .libraries are strings.
        return [getattr(x, "name", x) for x in dso.libraries]

    def canonicalise_library_spelling(libs):
        # This adjusts the library "spelling" so that it just contains the
        # name given to the linker. e.g. `@rpath/somewhere/libfoo.so.1.3`
        # would be canonicalised to "foo".
        fixes = []
        for lib in libs:
            # some libraries, e.g. Mach-O have an @rpath or system path
            # prefix in their name, remove it.
            path_stripped = os.path.split(lib)[-1]
            # Assume all library names contain at least one dot, even if they
            # don't it's fine, the first part is the piece of interest.
            prefix_libname = path_stripped.split(".")[0]
            linker_name = prefix_libname.replace("lib", "").replace("LIB", "")
            # further canonicalize by referring to all libraries in lower case.
            fixes.append(linker_name.lower())
        return fixes

    if HAVE_LIEF:
        dso = lief.parse(d["ffi_lib_location"])
        link_libs = tuple(canonicalise_library_type(dso))
        d["linked_libraries"] = link_libs
        canonicalised_libs = canonicalise_library_spelling(link_libs)
        d["canonicalised_linked_libraries"] = canonicalised_libs

    return d
