"""PyTables, hierarchical datasets in Python.

:URL: http://www.pytables.org/

PyTables is a package for managing hierarchical datasets and designed
to efficiently cope with extremely large amounts of data.

"""
import os
from ctypes import cdll
from ctypes.util import find_library
import platform


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
blosc2_search_paths = [blosc2_lib_hardcoded,
                       os.path.join(current_dir, blosc2_lib_hardcoded),
                       find_library("blosc2")]
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
    raise RuntimeError("Blosc2 library not found. "
                       f"I looked for \"{', '.join(blosc2_search_paths)}\"")

# Necessary imports to get versions stored on the cython extension
from .utilsextension import get_hdf5_version as _get_hdf5_version

from ._version import __version__

hdf5_version = _get_hdf5_version()
"""The underlying HDF5 library version number.

.. versionadded:: 3.0

"""

from .utilsextension import (
    blosc_compcode_to_compname_ as blosc_compcode_to_compname,
    blosc2_compcode_to_compname_ as blosc2_compcode_to_compname,
    blosc_get_complib_info_ as blosc_get_complib_info,
    blosc2_get_complib_info_ as blosc2_get_complib_info,
)

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

from .misc.enum import Enum
from .atom import *
from .flavor import restrict_flavors
from .description import *
from .filters import Filters

# Import the user classes from the proper modules
from .exceptions import *
from .file import File, open_file, copy_file
from .node import Node
from .group import Group
from .leaf import Leaf
from .table import Table, Cols, Column
from .array import Array
from .carray import CArray
from .earray import EArray
from .vlarray import VLArray
from .unimplemented import UnImplemented, Unknown
from .expression import Expr
from .tests import print_versions, test


# List here only the objects we want to be publicly available
__all__ = [
    # Exceptions and warnings:
    'HDF5ExtError',
    'ClosedNodeError', 'ClosedFileError', 'FileModeError',
    'NaturalNameWarning', 'NodeError', 'NoSuchNodeError',
    'UndoRedoError', 'UndoRedoWarning',
    'PerformanceWarning',
    'FlavorError', 'FlavorWarning',
    'FiltersWarning', 'DataTypeWarning',
    # Functions:
    'is_hdf5_file', 'is_pytables_file', 'which_lib_version',
    'copy_file', 'open_file', 'print_versions', 'test',
    'split_type', 'restrict_flavors',
    'set_blosc_max_threads', 'set_blosc2_max_threads',
    'silence_hdf5_messages',
    # Helper classes:
    'IsDescription', 'Description', 'Filters', 'Cols', 'Column',
    # Types:
    'Enum',
    # Atom types:
    'Atom', 'StringAtom', 'BoolAtom',
    'IntAtom', 'UIntAtom', 'Int8Atom', 'UInt8Atom', 'Int16Atom', 'UInt16Atom',
    'Int32Atom', 'UInt32Atom', 'Int64Atom', 'UInt64Atom',
    'FloatAtom', 'Float32Atom', 'Float64Atom',
    'ComplexAtom', 'Complex32Atom', 'Complex64Atom', 'Complex128Atom',
    'TimeAtom', 'Time32Atom', 'Time64Atom',
    'EnumAtom',
    'PseudoAtom', 'ObjectAtom', 'VLStringAtom', 'VLUnicodeAtom',
    # Column types:
    'Col', 'StringCol', 'BoolCol',
    'IntCol', 'UIntCol', 'Int8Col', 'UInt8Col', 'Int16Col', 'UInt16Col',
    'Int32Col', 'UInt32Col', 'Int64Col', 'UInt64Col',
    'FloatCol', 'Float32Col', 'Float64Col',
    'ComplexCol', 'Complex32Col', 'Complex64Col', 'Complex128Col',
    'TimeCol', 'Time32Col', 'Time64Col',
    'EnumCol',
    # Node classes:
    'Node', 'Group', 'Leaf', 'Table', 'Array', 'CArray', 'EArray', 'VLArray',
    'UnImplemented', 'Unknown',
    # The File class:
    'File',
    # Expr class
    'Expr',
]

if 'Float16Atom' in locals():
    # float16 is new in numpy 1.6.0
    __all__.extend(('Float16Atom', 'Float16Col'))

if 'Float96Atom' in locals():
    __all__.extend(('Float96Atom', 'Float96Col'))
    __all__.extend(('Complex192Atom', 'Complex192Col'))    # XXX check

if 'Float128Atom' in locals():
    __all__.extend(('Float128Atom', 'Float128Col'))
    __all__.extend(('Complex256Atom', 'Complex256Col'))    # XXX check


def get_pytables_version():
    warnings.warn(
        "the 'get_pytables_version()' function is deprecated and could be "
        "removed in future versions. Please use 'tables.__version__'",
        DeprecationWarning)
    return __version__

def get_hdf5_version():
    warnings.warn(
        "the 'get_hdf5_version()' function is deprecated and could be "
        "removed in future versions. Please use 'tables.hdf5_version'",
        DeprecationWarning)
    return hdf5_version