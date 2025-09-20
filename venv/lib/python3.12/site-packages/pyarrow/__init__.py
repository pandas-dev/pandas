# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

# flake8: noqa

"""
PyArrow is the python implementation of Apache Arrow.

Apache Arrow is a cross-language development platform for in-memory data.
It specifies a standardized language-independent columnar memory format for
flat and hierarchical data, organized for efficient analytic operations on
modern hardware. It also provides computational libraries and zero-copy
streaming messaging and interprocess communication.

For more information see the official page at https://arrow.apache.org
"""

import gc as _gc
import importlib as _importlib
import os as _os
import platform as _platform
import sys as _sys
import warnings as _warnings

try:
    from ._generated_version import version as __version__
except ImportError:
    # Package is not installed, parse git tag at runtime
    try:
        import setuptools_scm
        # Code duplicated from setup.py to avoid a dependency on each other

        def parse_git(root, **kwargs):
            """
            Parse function for setuptools_scm that ignores tags for non-C++
            subprojects, e.g. apache-arrow-js-XXX tags.
            """
            from setuptools_scm.git import parse
            kwargs['describe_command'] = \
                "git describe --dirty --tags --long --match 'apache-arrow-[0-9]*.*'"
            return parse(root, **kwargs)
        __version__ = setuptools_scm.get_version('../',
                                                 parse=parse_git)
    except ImportError:
        __version__ = None

import pyarrow.lib as _lib
from pyarrow.lib import (BuildInfo, RuntimeInfo, set_timezone_db_path,
                         MonthDayNano, VersionInfo, cpp_build_info,
                         cpp_version, cpp_version_info, runtime_info,
                         cpu_count, set_cpu_count, enable_signal_handlers,
                         io_thread_count, set_io_thread_count)


def show_versions():
    """
    Print various version information, to help with error reporting.
    """
    def print_entry(label, value):
        print(f"{label: <26}: {value: <8}")

    print("pyarrow version info\n--------------------")
    print_entry("Package kind", cpp_build_info.package_kind
                if len(cpp_build_info.package_kind) > 0
                else "not indicated")
    print_entry("Arrow C++ library version", cpp_build_info.version)
    print_entry("Arrow C++ compiler",
                f"{cpp_build_info.compiler_id} {cpp_build_info.compiler_version}")
    print_entry("Arrow C++ compiler flags", cpp_build_info.compiler_flags)
    print_entry("Arrow C++ git revision", cpp_build_info.git_id)
    print_entry("Arrow C++ git description", cpp_build_info.git_description)
    print_entry("Arrow C++ build type", cpp_build_info.build_type)


def _module_is_available(module):
    try:
        _importlib.import_module(f'pyarrow.{module}')
    except ImportError:
        return False
    else:
        return True


def _filesystem_is_available(fs):
    try:
        import pyarrow.fs
    except ImportError:
        return False

    try:
        getattr(pyarrow.fs, fs)
    except (ImportError, AttributeError):
        return False
    else:
        return True


def show_info():
    """
    Print detailed version and platform information, for error reporting
    """
    show_versions()

    def print_entry(label, value):
        print(f"  {label: <20}: {value: <8}")

    print("\nPlatform:")
    print_entry("OS / Arch", f"{_platform.system()} {_platform.machine()}")
    print_entry("SIMD Level", runtime_info().simd_level)
    print_entry("Detected SIMD Level", runtime_info().detected_simd_level)

    pool = default_memory_pool()
    print("\nMemory:")
    print_entry("Default backend", pool.backend_name)
    print_entry("Bytes allocated", f"{pool.bytes_allocated()} bytes")
    print_entry("Max memory", f"{pool.max_memory()} bytes")
    print_entry("Supported Backends", ', '.join(supported_memory_backends()))

    print("\nOptional modules:")
    modules = ["csv", "cuda", "dataset", "feather", "flight", "fs", "gandiva", "json",
               "orc", "parquet"]
    for module in modules:
        status = "Enabled" if _module_is_available(module) else "-"
        print(f"  {module: <20}: {status: <8}")

    print("\nFilesystems:")
    filesystems = ["AzureFileSystem", "GcsFileSystem",
                   "HadoopFileSystem", "S3FileSystem"]
    for fs in filesystems:
        status = "Enabled" if _filesystem_is_available(fs) else "-"
        print(f"  {fs: <20}: {status: <8}")

    print("\nCompression Codecs:")
    codecs = ["brotli", "bz2", "gzip", "lz4_frame", "lz4", "snappy", "zstd"]
    for codec in codecs:
        status = "Enabled" if Codec.is_available(codec) else "-"
        print(f"  {codec: <20}: {status: <8}")


from pyarrow.lib import (null, bool_,
                         int8, int16, int32, int64,
                         uint8, uint16, uint32, uint64,
                         time32, time64, timestamp, date32, date64, duration,
                         month_day_nano_interval,
                         float16, float32, float64,
                         binary, string, utf8, binary_view, string_view,
                         large_binary, large_string, large_utf8,
                         decimal32, decimal64, decimal128, decimal256,
                         list_, large_list, list_view, large_list_view,
                         map_, struct,
                         union, sparse_union, dense_union,
                         dictionary,
                         run_end_encoded,
                         bool8, fixed_shape_tensor, json_, opaque, uuid,
                         field,
                         type_for_alias,
                         DataType, DictionaryType, StructType,
                         ListType, LargeListType, FixedSizeListType,
                         ListViewType, LargeListViewType,
                         MapType, UnionType, SparseUnionType, DenseUnionType,
                         TimestampType, Time32Type, Time64Type, DurationType,
                         FixedSizeBinaryType,
                         Decimal32Type, Decimal64Type, Decimal128Type, Decimal256Type,
                         BaseExtensionType, ExtensionType,
                         RunEndEncodedType, Bool8Type, FixedShapeTensorType,
                         JsonType, OpaqueType, UuidType,
                         UnknownExtensionType,
                         register_extension_type, unregister_extension_type,
                         DictionaryMemo,
                         KeyValueMetadata,
                         Field,
                         Schema,
                         schema,
                         unify_schemas,
                         Array, Tensor,
                         array, chunked_array, record_batch, nulls, repeat,
                         SparseCOOTensor, SparseCSRMatrix, SparseCSCMatrix,
                         SparseCSFTensor,
                         infer_type, from_numpy_dtype,
                         arange,
                         NullArray,
                         NumericArray, IntegerArray, FloatingPointArray,
                         BooleanArray,
                         Int8Array, UInt8Array,
                         Int16Array, UInt16Array,
                         Int32Array, UInt32Array,
                         Int64Array, UInt64Array,
                         HalfFloatArray, FloatArray, DoubleArray,
                         ListArray, LargeListArray, FixedSizeListArray,
                         ListViewArray, LargeListViewArray,
                         MapArray, UnionArray,
                         BinaryArray, StringArray,
                         LargeBinaryArray, LargeStringArray,
                         BinaryViewArray, StringViewArray,
                         FixedSizeBinaryArray,
                         DictionaryArray,
                         Date32Array, Date64Array, TimestampArray,
                         Time32Array, Time64Array, DurationArray,
                         MonthDayNanoIntervalArray,
                         Decimal32Array, Decimal64Array, Decimal128Array, Decimal256Array,
                         StructArray, ExtensionArray,
                         RunEndEncodedArray, Bool8Array, FixedShapeTensorArray,
                         JsonArray, OpaqueArray, UuidArray,
                         scalar, NA, _NULL as NULL, Scalar,
                         NullScalar, BooleanScalar,
                         Int8Scalar, Int16Scalar, Int32Scalar, Int64Scalar,
                         UInt8Scalar, UInt16Scalar, UInt32Scalar, UInt64Scalar,
                         HalfFloatScalar, FloatScalar, DoubleScalar,
                         Decimal32Scalar, Decimal64Scalar, Decimal128Scalar, Decimal256Scalar,
                         ListScalar, LargeListScalar, FixedSizeListScalar,
                         ListViewScalar, LargeListViewScalar,
                         Date32Scalar, Date64Scalar,
                         Time32Scalar, Time64Scalar,
                         TimestampScalar, DurationScalar,
                         MonthDayNanoIntervalScalar,
                         BinaryScalar, LargeBinaryScalar, BinaryViewScalar,
                         StringScalar, LargeStringScalar, StringViewScalar,
                         FixedSizeBinaryScalar, DictionaryScalar,
                         MapScalar, StructScalar, UnionScalar,
                         RunEndEncodedScalar, Bool8Scalar, ExtensionScalar,
                         FixedShapeTensorScalar, JsonScalar, OpaqueScalar, UuidScalar)

# Buffers, allocation
from pyarrow.lib import (DeviceAllocationType, Device, MemoryManager,
                         default_cpu_memory_manager)

from pyarrow.lib import (Buffer, ResizableBuffer, foreign_buffer, py_buffer,
                         Codec, compress, decompress, allocate_buffer)

from pyarrow.lib import (MemoryPool, LoggingMemoryPool, ProxyMemoryPool,
                         total_allocated_bytes, set_memory_pool,
                         default_memory_pool, system_memory_pool,
                         jemalloc_memory_pool, mimalloc_memory_pool,
                         logging_memory_pool, proxy_memory_pool,
                         log_memory_allocations, jemalloc_set_decay_ms,
                         supported_memory_backends)

# I/O
from pyarrow.lib import (NativeFile, PythonFile,
                         BufferedInputStream, BufferedOutputStream, CacheOptions,
                         CompressedInputStream, CompressedOutputStream,
                         TransformInputStream, transcoding_input_stream,
                         FixedSizeBufferWriter,
                         BufferReader, BufferOutputStream,
                         OSFile, MemoryMappedFile, memory_map,
                         create_memory_map, MockOutputStream,
                         input_stream, output_stream,
                         have_libhdfs)

from pyarrow.lib import (ChunkedArray, RecordBatch, Table, table,
                         concat_arrays, concat_tables, TableGroupBy,
                         RecordBatchReader, concat_batches)

# Exceptions
from pyarrow.lib import (ArrowCancelled,
                         ArrowCapacityError,
                         ArrowException,
                         ArrowKeyError,
                         ArrowIndexError,
                         ArrowInvalid,
                         ArrowIOError,
                         ArrowMemoryError,
                         ArrowNotImplementedError,
                         ArrowTypeError,
                         ArrowSerializationError)

from pyarrow.ipc import serialize_pandas, deserialize_pandas
import pyarrow.ipc as ipc

import pyarrow.types as types


# ----------------------------------------------------------------------
# Deprecations

from pyarrow.util import _deprecate_api, _deprecate_class


# TODO: Deprecate these somehow in the pyarrow namespace
from pyarrow.ipc import (Message, MessageReader, MetadataVersion,
                         RecordBatchFileReader, RecordBatchFileWriter,
                         RecordBatchStreamReader, RecordBatchStreamWriter)

# ----------------------------------------------------------------------
# Returning absolute path to the pyarrow include directory (if bundled, e.g. in
# wheels)


def get_include():
    """
    Return absolute path to directory containing Arrow C++ include
    headers. Similar to numpy.get_include
    """
    return _os.path.join(_os.path.dirname(__file__), 'include')


def _get_pkg_config_executable():
    return _os.environ.get('PKG_CONFIG', 'pkg-config')


def _has_pkg_config(pkgname):
    import subprocess
    try:
        return subprocess.call([_get_pkg_config_executable(),
                                '--exists', pkgname]) == 0
    except FileNotFoundError:
        return False


def _read_pkg_config_variable(pkgname, cli_args):
    import subprocess
    cmd = [_get_pkg_config_executable(), pkgname] + cli_args
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("pkg-config failed: " + err.decode('utf8'))
    return out.rstrip().decode('utf8')


def get_libraries():
    """
    Return list of library names to include in the `libraries` argument for C
    or Cython extensions using pyarrow
    """
    return ['arrow_python', 'arrow']


def create_library_symlinks():
    """
    With Linux and macOS wheels, the bundled shared libraries have an embedded
    ABI version like libarrow.so.17 or libarrow.17.dylib and so linking to them
    with -larrow won't work unless we create symlinks at locations like
    site-packages/pyarrow/libarrow.so. This unfortunate workaround addresses
    prior problems we had with shipping two copies of the shared libraries to
    permit third party projects like turbodbc to build their C++ extensions
    against the pyarrow wheels.

    This function must only be invoked once and only when the shared libraries
    are bundled with the Python package, which should only apply to wheel-based
    installs. It requires write access to the site-packages/pyarrow directory
    and so depending on your system may need to be run with root.
    """
    import glob
    if _sys.platform == 'win32':
        return
    package_cwd = _os.path.dirname(__file__)

    if _sys.platform == 'linux':
        bundled_libs = glob.glob(_os.path.join(package_cwd, '*.so.*'))

        def get_symlink_path(hard_path):
            return hard_path.rsplit('.', 1)[0]
    else:
        bundled_libs = glob.glob(_os.path.join(package_cwd, '*.*.dylib'))

        def get_symlink_path(hard_path):
            return '.'.join((hard_path.rsplit('.', 2)[0], 'dylib'))

    for lib_hard_path in bundled_libs:
        symlink_path = get_symlink_path(lib_hard_path)
        if _os.path.exists(symlink_path):
            continue
        try:
            _os.symlink(lib_hard_path, symlink_path)
        except PermissionError:
            print("Tried creating symlink {}. If you need to link to "
                  "bundled shared libraries, run "
                  "pyarrow.create_library_symlinks() as root")


def get_library_dirs():
    """
    Return lists of directories likely to contain Arrow C++ libraries for
    linking C or Cython extensions using pyarrow
    """
    package_cwd = _os.path.dirname(__file__)
    library_dirs = [package_cwd]

    def append_library_dir(library_dir):
        if library_dir not in library_dirs:
            library_dirs.append(library_dir)

    # Search library paths via pkg-config. This is necessary if the user
    # installed libarrow and the other shared libraries manually and they
    # are not shipped inside the pyarrow package (see also ARROW-2976).
    pkg_config_executable = _os.environ.get('PKG_CONFIG') or 'pkg-config'
    for pkgname in ["arrow", "arrow_python"]:
        if _has_pkg_config(pkgname):
            library_dir = _read_pkg_config_variable(pkgname,
                                                    ["--libs-only-L"])
            # pkg-config output could be empty if Arrow is installed
            # as a system package.
            if library_dir:
                if not library_dir.startswith("-L"):
                    raise ValueError(
                        "pkg-config --libs-only-L returned unexpected "
                        f"value {library_dir!r}")
                append_library_dir(library_dir[2:])

    if _sys.platform == 'win32':
        # TODO(wesm): Is this necessary, or does setuptools within a conda
        # installation add Library\lib to the linker path for MSVC?
        python_base_install = _os.path.dirname(_sys.executable)
        library_dir = _os.path.join(python_base_install, 'Library', 'lib')

        if _os.path.exists(_os.path.join(library_dir, 'arrow.lib')):
            append_library_dir(library_dir)

        # GH-45530: Add pyarrow.libs dir containing delvewheel-mangled
        # msvcp140.dll
        pyarrow_libs_dir = _os.path.abspath(
            _os.path.join(_os.path.dirname(__file__), _os.pardir, "pyarrow.libs")
        )
        if _os.path.exists(pyarrow_libs_dir):
            append_library_dir(pyarrow_libs_dir)

    # ARROW-4074: Allow for ARROW_HOME to be set to some other directory
    if _os.environ.get('ARROW_HOME'):
        append_library_dir(_os.path.join(_os.environ['ARROW_HOME'], 'lib'))
    else:
        # Python wheels bundle the Arrow libraries in the pyarrow directory.
        append_library_dir(_os.path.dirname(_os.path.abspath(__file__)))

    return library_dirs
