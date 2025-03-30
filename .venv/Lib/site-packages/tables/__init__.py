"""PyTables, hierarchical datasets in Python.

:URL: http://www.pytables.org/

PyTables is a package for managing hierarchical datasets and designed
to efficiently cope with extremely large amounts of data.

"""


# start delvewheel patch
def _delvewheel_patch_1_9_0():
    import os
    if os.path.isdir(libs_dir := os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'tables.libs'))):
        os.add_dll_directory(libs_dir)


_delvewheel_patch_1_9_0()
del _delvewheel_patch_1_9_0
# end delvewheel patch

import os
import platform
from ctypes import cdll
from ctypes.util import find_library

# Load the blosc2 library, and if not found in standard locations,
# try this directory (it should be automatically copied in setup.py).
current_dir = os.path.dirname(__file__)
platform_system = platform.system()
blosc2_lib_hardcoded = "libblosc2"
if platform_system == "Linux":
    blosc2_lib_hardcoded += ".so"
elif platform_system == "Darwin":
    blosc2_lib_hardcoded += ".dylib"
else:
    blosc2_lib_hardcoded += ".dll"
blosc2_found = False
blosc2_search_paths = [
    blosc2_lib_hardcoded,
    os.path.join(current_dir, blosc2_lib_hardcoded),
    # delvewheel will put it here
    os.path.join(
        os.path.dirname(current_dir), "tables.libs", blosc2_lib_hardcoded
    ),
]
if find_library("blosc2"):
    blosc2_search_paths.append(find_library("blosc2"))
for blosc2_lib in blosc2_search_paths:
    if blosc2_lib:
        try:
            cdll.LoadLibrary(blosc2_lib)
        except OSError:
            pass
        else:
            blosc2_found = True
            break
if not blosc2_found:
    raise RuntimeError(
        f"Blosc2 library not found. "
        f"I looked for {', '.join(blosc2_search_paths)}"
    )

from ._version import __version__

# Necessary imports to get versions stored on the cython extension
from .utilsextension import get_hdf5_version as _get_hdf5_version

hdf5_version = _get_hdf5_version()
"""The underlying HDF5 library version number.

.. versionadded:: 3.0

"""

from .atom import *
from .file import File, open_file, copy_file
from .leaf import Leaf, ChunkInfo
from .node import Node
from .array import Array
from .group import Group
from .table import Table, Cols, Column
from .tests import print_versions, test
from .carray import CArray
from .earray import EArray
from .flavor import restrict_flavors
from .filters import Filters
from .vlarray import VLArray
from .misc.enum import Enum

# Import the user classes from the proper modules
from .exceptions import *
from .expression import Expr
from .description import *
from .unimplemented import UnImplemented, Unknown
from .utilsextension import (
    blosc_compcode_to_compname_ as blosc_compcode_to_compname,
)
from .utilsextension import (
    blosc2_compcode_to_compname_ as blosc2_compcode_to_compname,
)
from .utilsextension import blosc_get_complib_info_ as blosc_get_complib_info
from .utilsextension import blosc2_get_complib_info_ as blosc2_get_complib_info
from .utilsextension import (
    blosc_compressor_list,
    blosc2_compressor_list,
    is_hdf5_file,
    is_pytables_file,
    which_lib_version,
    set_blosc_max_threads,
    set_blosc2_max_threads,
    silence_hdf5_messages,
)

# List here only the objects we want to be publicly available
__all__ = [
    # Exceptions and warnings:
    "HDF5ExtError",
    "ClosedNodeError",
    "ClosedFileError",
    "FileModeError",
    "NaturalNameWarning",
    "NodeError",
    "NoSuchNodeError",
    "UndoRedoError",
    "UndoRedoWarning",
    "PerformanceWarning",
    "FlavorError",
    "FlavorWarning",
    "FiltersWarning",
    "DataTypeWarning",
    "ChunkError",
    "NotChunkedError",
    "NotChunkAlignedError",
    "NoSuchChunkError",
    # Functions:
    "is_hdf5_file",
    "is_pytables_file",
    "which_lib_version",
    "copy_file",
    "open_file",
    "print_versions",
    "test",
    "split_type",
    "restrict_flavors",
    "set_blosc_max_threads",
    "set_blosc2_max_threads",
    "silence_hdf5_messages",
    # Helper classes:
    "IsDescription",
    "Description",
    "Filters",
    "Cols",
    "Column",
    "ChunkInfo",
    # Types:
    "Enum",
    # Atom types:
    "Atom",
    "StringAtom",
    "BoolAtom",
    "IntAtom",
    "UIntAtom",
    "Int8Atom",
    "UInt8Atom",
    "Int16Atom",
    "UInt16Atom",
    "Int32Atom",
    "UInt32Atom",
    "Int64Atom",
    "UInt64Atom",
    "FloatAtom",
    "Float32Atom",
    "Float64Atom",
    "ComplexAtom",
    "Complex32Atom",
    "Complex64Atom",
    "Complex128Atom",
    "TimeAtom",
    "Time32Atom",
    "Time64Atom",
    "EnumAtom",
    "PseudoAtom",
    "ObjectAtom",
    "VLStringAtom",
    "VLUnicodeAtom",
    # Column types:
    "Col",
    "StringCol",
    "BoolCol",
    "IntCol",
    "UIntCol",
    "Int8Col",
    "UInt8Col",
    "Int16Col",
    "UInt16Col",
    "Int32Col",
    "UInt32Col",
    "Int64Col",
    "UInt64Col",
    "FloatCol",
    "Float32Col",
    "Float64Col",
    "ComplexCol",
    "Complex32Col",
    "Complex64Col",
    "Complex128Col",
    "TimeCol",
    "Time32Col",
    "Time64Col",
    "EnumCol",
    # Node classes:
    "Node",
    "Group",
    "Leaf",
    "Table",
    "Array",
    "CArray",
    "EArray",
    "VLArray",
    "UnImplemented",
    "Unknown",
    # The File class:
    "File",
    # Expr class
    "Expr",
]

if "Float16Atom" in locals():
    # float16 is new in numpy 1.6.0
    __all__.extend(("Float16Atom", "Float16Col"))

if "Float96Atom" in locals():
    __all__.extend(("Float96Atom", "Float96Col"))
    __all__.extend(("Complex192Atom", "Complex192Col"))  # XXX check

if "Float128Atom" in locals():
    __all__.extend(("Float128Atom", "Float128Col"))
    __all__.extend(("Complex256Atom", "Complex256Col"))  # XXX check


def get_pytables_version() -> str:
    warnings.warn(
        "the 'get_pytables_version()' function is deprecated and could be "
        "removed in future versions. Please use 'tables.__version__'",
        DeprecationWarning,
    )
    return __version__


def get_hdf5_version() -> str:
    warnings.warn(
        "the 'get_hdf5_version()' function is deprecated and could be "
        "removed in future versions. Please use 'tables.hdf5_version'",
        DeprecationWarning,
    )
    return hdf5_version