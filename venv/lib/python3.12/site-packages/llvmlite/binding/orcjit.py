import ctypes
from ctypes import POINTER, c_bool, c_char_p, c_uint8, c_uint64, c_size_t

from llvmlite.binding import ffi, targets


class _LinkElement(ctypes.Structure):
    _fields_ = [("element_kind", c_uint8),
                ("value", c_char_p),
                ("value_len", c_size_t)]


class _SymbolAddress(ctypes.Structure):
    _fields_ = [("name", c_char_p), ("address", c_uint64)]


class JITLibraryBuilder:
    """
    Create a library for linking by OrcJIT

    OrcJIT operates like a linker: a number of compilation units and
    dependencies are collected together and linked into a single dynamic library
    that can export functions to other libraries or to be consumed directly as
    entry points into JITted code. The native OrcJIT has a lot of memory
    management complications so this API is designed to work well with Python's
    garbage collection.

    The creation of a new library is a bit like a linker command line where
    compilation units, mostly as LLVM IR, and previously constructed libraries
    are linked together, then loaded into memory, and the addresses of exported
    symbols are extracted. Any static initializers are run and the exported
    addresses and a resource tracker is produced. As long as the resource
    tracker is referenced somewhere in Python, the exported addresses will be
    valid. Once the resource tracker is garbage collected, the static
    destructors will run and library will be unloaded from memory.
    """
    def __init__(self):
        self.__entries = []
        self.__exports = set()
        self.__imports = {}

    def add_ir(self, llvmir):
        """
        Adds a compilation unit to the library using LLVM IR as the input
        format.

        This takes a string or an object that can be converted to a string,
        including IRBuilder, that contains LLVM IR.
        """
        self.__entries.append((0, str(llvmir).encode('utf-8')))
        return self

    def add_native_assembly(self, asm):
        """
        Adds a compilation unit to the library using native assembly as the
        input format.

        This takes a string or an object that can be converted to a string that
        contains native assembly, which will be
        parsed by LLVM.
        """
        self.__entries.append((1, str(asm).encode('utf-8')))
        return self

    def add_object_img(self, data):
        """
        Adds a compilation unit to the library using pre-compiled object code.

        This takes the bytes of the contents of an object artifact which will be
        loaded by LLVM.
        """
        self.__entries.append((2, bytes(data)))
        return self

    def add_object_file(self, file_path):
        """
        Adds a compilation unit to the library using pre-compiled object file.

        This takes a string or path-like object that references an object file
        which will be loaded by LLVM.
        """
        with open(file_path, "rb") as f:
            self.__entries.append((2, f.read()))
        return self

    def add_jit_library(self, name):
        """
        Adds an existing JIT library as prerequisite.

        The name of the library must match the one provided in a previous link
        command.
        """
        self.__entries.append((3, str(name).encode('utf-8')))
        return self

    def add_current_process(self):
        """
        Allows the JITted library to access symbols in the current binary.

        That is, it allows exporting the current binary's symbols, including
        loaded libraries, as imports to the JITted
        library.
        """
        self.__entries.append((3, b''))
        return self

    def import_symbol(self, name, address):
        """
        Register the *address* of global symbol *name*.  This will make
        it usable (e.g. callable) from LLVM-compiled functions.
        """
        self.__imports[str(name)] = c_uint64(address)
        return self

    def export_symbol(self, name):
        """
        During linking, extract the address of a symbol that was defined in one
        of the compilation units.

        This allows getting symbols, functions or global variables, out of the
        JIT linked library. The addresses will be
        available when the link method is called.
        """
        self.__exports.add(str(name))
        return self

    def link(self, lljit, library_name):
        """
        Link all the current compilation units into a JITted library and extract
        the address of exported symbols.

        An instance of the OrcJIT instance must be provided and this will be the
        scope that is used to find other JITted libraries that are dependencies
        and also be the place where this library will be defined.

        After linking, the method will return a resource tracker that keeps the
        library alive. This tracker also knows the addresses of any exported
        symbols that were requested.

        The addresses will be valid as long as the resource tracker is
        referenced.

        When the resource tracker is destroyed, the library will be cleaned up,
        however, the name of the library cannot be reused.
        """
        assert not lljit.closed, "Cannot add to closed JIT"
        encoded_library_name = str(library_name).encode('utf-8')
        assert len(encoded_library_name) > 0, "Library cannot be empty"
        elements = (_LinkElement * len(self.__entries))()
        for idx, (kind, value) in enumerate(self.__entries):
            elements[idx].element_kind = c_uint8(kind)
            elements[idx].value = c_char_p(value)
            elements[idx].value_len = c_size_t(len(value))
        exports = (_SymbolAddress * len(self.__exports))()
        for idx, name in enumerate(self.__exports):
            exports[idx].name = name.encode('utf-8')

        imports = (_SymbolAddress * len(self.__imports))()
        for idx, (name, addr) in enumerate(self.__imports.items()):
            imports[idx].name = name.encode('utf-8')
            imports[idx].address = addr

        with ffi.OutputString() as outerr:
            tracker = lljit._capi.LLVMPY_LLJIT_Link(
                lljit._ptr,
                encoded_library_name,
                elements,
                len(self.__entries),
                imports,
                len(self.__imports),
                exports,
                len(self.__exports),
                outerr)
            if not tracker:
                raise RuntimeError(str(outerr))
        return ResourceTracker(tracker,
                               library_name,
                               {name: exports[idx].address
                                for idx, name in enumerate(self.__exports)})


class ResourceTracker(ffi.ObjectRef):
    """
    A resource tracker is created for each loaded JIT library and keeps the
    module alive.

    OrcJIT supports unloading libraries that are no longer used. This resource
    tracker should be stored in any object that reference functions or constants
    for a JITted library. When all references to the resource tracker are
    dropped, this will trigger LLVM to unload the library and destroy any
    functions.

    Failure to keep resource trackers while calling a function or accessing a
    symbol can result in crashes or memory corruption.

    LLVM internally tracks references between different libraries, so only
    "leaf" libraries need to be tracked.
    """
    def __init__(self, ptr, name, addresses):
        self.__addresses = addresses
        self.__name = name
        ffi.ObjectRef.__init__(self, ptr)

    def __getitem__(self, item):
        """
        Get the address of an exported symbol as an integer
        """
        return self.__addresses[item]

    @property
    def name(self):
        return self.__name

    def _dispose(self):
        with ffi.OutputString() as outerr:
            if self._capi.LLVMPY_LLJIT_Dylib_Tracker_Dispose(self, outerr):
                raise RuntimeError(str(outerr))


class LLJIT(ffi.ObjectRef):
    """
    A OrcJIT-based LLVM JIT engine that can compile and run LLVM IR as a
    collection of JITted dynamic libraries

    The C++ OrcJIT API has a lot of memory ownership patterns that do not work
    with Python. This API attempts to provide ones that are safe at the expense
    of some features. Each LLJIT instance is a collection of JIT-compiled
    libraries. In the C++ API, there is a "main" library; this API does not
    provide access to the main library. Use the JITLibraryBuilder to create a
    new named library instead.
    """
    def __init__(self, ptr):
        self._td = None
        ffi.ObjectRef.__init__(self, ptr)

    def lookup(self, dylib, fn):
        """
        Find a function in this dynamic library and construct a new tracking
        object for it

        If the library or function do not exist, an exception will occur.

        Parameters
        ----------
        dylib : str or None
           the name of the library containing the symbol
        fn : str
           the name of the function to get
        """
        assert not self.closed, "Cannot lookup in closed JIT"
        address = ctypes.c_uint64()
        with ffi.OutputString() as outerr:
            tracker = ffi.lib.LLVMPY_LLJITLookup(self,
                                                 dylib.encode("utf-8"),
                                                 fn.encode("utf-8"),
                                                 ctypes.byref(address),
                                                 outerr)
            if not tracker:
                raise RuntimeError(str(outerr))

        return ResourceTracker(tracker, dylib, {fn: address.value})

    @property
    def target_data(self):
        """
        The TargetData for this LLJIT instance.
        """
        if self._td is not None:
            return self._td
        ptr = ffi.lib.LLVMPY_LLJITGetDataLayout(self)
        self._td = targets.TargetData(ptr)
        self._td._owned = True
        return self._td

    def _dispose(self):
        if self._td is not None:
            self._td.detach()
        self._capi.LLVMPY_LLJITDispose(self)


def create_lljit_compiler(target_machine=None, *,
                          use_jit_link=False,
                          suppress_errors=False):
    """
    Create an LLJIT instance
    """
    with ffi.OutputString() as outerr:
        lljit = ffi.lib.LLVMPY_CreateLLJITCompiler(target_machine,
                                                   suppress_errors,
                                                   use_jit_link,
                                                   outerr)
        if not lljit:
            raise RuntimeError(str(outerr))

    return LLJIT(lljit)


ffi.lib.LLVMPY_LLJITLookup.argtypes = [
    ffi.LLVMOrcLLJITRef,
    c_char_p,
    c_char_p,
    POINTER(c_uint64),
    POINTER(c_char_p),
]
ffi.lib.LLVMPY_LLJITLookup.restype = ffi.LLVMOrcDylibTrackerRef

ffi.lib.LLVMPY_LLJITGetDataLayout.argtypes = [
    ffi.LLVMOrcLLJITRef,
]
ffi.lib.LLVMPY_LLJITGetDataLayout.restype = ffi.LLVMTargetDataRef

ffi.lib.LLVMPY_CreateLLJITCompiler.argtypes = [
    ffi.LLVMTargetMachineRef,
    c_bool,
    c_bool,
    POINTER(c_char_p),
]
ffi.lib.LLVMPY_CreateLLJITCompiler.restype = ffi.LLVMOrcLLJITRef

ffi.lib.LLVMPY_LLJITDispose.argtypes = [
    ffi.LLVMOrcLLJITRef,
]


ffi.lib.LLVMPY_LLJIT_Link.argtypes = [
    ffi.LLVMOrcLLJITRef,
    c_char_p,
    POINTER(_LinkElement),
    c_size_t,
    POINTER(_SymbolAddress),
    c_size_t,
    POINTER(_SymbolAddress),
    c_size_t,
    POINTER(c_char_p)
]
ffi.lib.LLVMPY_LLJIT_Link.restype = ffi.LLVMOrcDylibTrackerRef

ffi.lib.LLVMPY_LLJIT_Dylib_Tracker_Dispose.argtypes = [
    ffi.LLVMOrcDylibTrackerRef,
    POINTER(c_char_p)
]
ffi.lib.LLVMPY_LLJIT_Dylib_Tracker_Dispose.restype = c_bool
