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

import os
import subprocess
import sys

import pytest

import pyarrow as pa
from pyarrow.lib import ArrowInvalid


def test_get_include():
    include_dir = pa.get_include()
    assert os.path.exists(os.path.join(include_dir, 'arrow', 'api.h'))


@pytest.mark.skipif('sys.platform != "win32"')
def test_get_library_dirs_win32():
    assert any(os.path.exists(os.path.join(directory, 'arrow.lib'))
               for directory in pa.get_library_dirs())


def test_cpu_count():
    n = pa.cpu_count()
    assert n > 0
    try:
        pa.set_cpu_count(n + 5)
        assert pa.cpu_count() == n + 5
    finally:
        pa.set_cpu_count(n)


def test_io_thread_count():
    n = pa.io_thread_count()
    assert n > 0
    try:
        pa.set_io_thread_count(n + 5)
        assert pa.io_thread_count() == n + 5
    finally:
        pa.set_io_thread_count(n)


@pytest.mark.processes
def test_env_var_io_thread_count():
    # Test that the number of IO threads can be overridden with the
    # ARROW_IO_THREADS environment variable.
    code = """if 1:
        import pyarrow as pa
        print(pa.io_thread_count())
        """

    def run_with_env_var(env_var):
        env = os.environ.copy()
        env['ARROW_IO_THREADS'] = env_var
        res = subprocess.run([sys.executable, "-c", code], env=env,
                             capture_output=True)
        res.check_returncode()
        return res.stdout.decode(), res.stderr.decode()

    out, err = run_with_env_var('17')
    assert out.strip() == '17'
    assert err == ''

    for v in ('-1', 'z'):
        out, err = run_with_env_var(v)
        assert out.strip() == '8'  # default value
        assert ("ARROW_IO_THREADS does not contain a valid number of threads"
                in err.strip())


def test_build_info():
    assert isinstance(pa.build_info.cpp_build_info, pa.CppBuildInfo)
    assert isinstance(pa.cpp_version_info, pa.VersionInfo)
    assert isinstance(pa.cpp_version, str)
    assert isinstance(pa.__version__, str)
    assert pa.build_info.cpp_build_info.version == pa.cpp_version
    assert pa.build_info.cpp_build_info.version_info == pa.cpp_version_info
    assert pa.build_info.cpp_build_info is pa.cpp_build_info

    assert pa.build_info.cpp_build_info.build_type in (
        'debug', 'release', 'minsizerel', 'relwithdebinfo')

    assert isinstance(pa.build_info, pa.BuildInfo)
    assert pa.build_info.build_type in (
        'debug', 'release', 'minsizerel', 'relwithdebinfo')

    # assert pa.version == pa.__version__  # XXX currently false


def test_runtime_info():
    info = pa.runtime_info()
    assert isinstance(info, pa.RuntimeInfo)
    possible_simd_levels = ('none', 'sse4_2', 'avx', 'avx2', 'avx512')
    assert info.simd_level in possible_simd_levels
    assert info.detected_simd_level in possible_simd_levels

    if info.simd_level != 'none':
        env = os.environ.copy()
        env['ARROW_USER_SIMD_LEVEL'] = 'none'
        code = f"""if 1:
            import pyarrow as pa

            info = pa.runtime_info()
            assert info.simd_level == 'none', info.simd_level
            assert info.detected_simd_level == {info.detected_simd_level!r},\
                info.detected_simd_level
            """
        subprocess.check_call([sys.executable, "-c", code], env=env)


@pytest.mark.processes
def test_import_at_shutdown():
    # GH-38626: importing PyArrow at interpreter shutdown would crash
    code = """if 1:
        import atexit

        def import_arrow():
            import pyarrow

        atexit.register(import_arrow)
        """
    subprocess.check_call([sys.executable, "-c", code])


@pytest.mark.skipif(sys.platform == "win32",
                    reason="Path to timezone database is not configurable "
                           "on non-Windows platforms")
def test_set_timezone_db_path_non_windows():
    # set_timezone_db_path raises an error on non-Windows platforms
    with pytest.raises(ArrowInvalid,
                       match="Arrow was set to use OS timezone "
                             "database at compile time"):
        pa.set_timezone_db_path("path")


@pytest.mark.parametrize('klass', [
    pa.Field,
    pa.Schema,
    pa.ChunkedArray,
    pa.RecordBatch,
    pa.Table,
    pa.Buffer,
    pa.Array,
    pa.Tensor,
    pa.DataType,
    pa.ListType,
    pa.LargeListType,
    pa.FixedSizeListType,
    pa.ListViewType,
    pa.LargeListViewType,
    pa.UnionType,
    pa.SparseUnionType,
    pa.DenseUnionType,
    pa.StructType,
    pa.Time32Type,
    pa.Time64Type,
    pa.TimestampType,
    pa.Decimal32Type,
    pa.Decimal64Type,
    pa.Decimal128Type,
    pa.Decimal256Type,
    pa.DictionaryType,
    pa.FixedSizeBinaryType,
    pa.NullArray,
    pa.NumericArray,
    pa.IntegerArray,
    pa.FloatingPointArray,
    pa.BooleanArray,
    pa.Int8Array,
    pa.Int16Array,
    pa.Int32Array,
    pa.Int64Array,
    pa.UInt8Array,
    pa.UInt16Array,
    pa.UInt32Array,
    pa.UInt64Array,
    pa.ListArray,
    pa.LargeListArray,
    pa.MapArray,
    pa.FixedSizeListArray,
    pa.UnionArray,
    pa.BinaryArray,
    pa.StringArray,
    pa.BinaryViewArray,
    pa.StringViewArray,
    pa.FixedSizeBinaryArray,
    pa.DictionaryArray,
    pa.Date32Array,
    pa.Date64Array,
    pa.TimestampArray,
    pa.Time32Array,
    pa.Time64Array,
    pa.DurationArray,
    pa.Decimal128Array,
    pa.Decimal256Array,
    pa.StructArray,
    pa.RunEndEncodedArray,
    pa.Scalar,
    pa.BooleanScalar,
    pa.Int8Scalar,
    pa.Int16Scalar,
    pa.Int32Scalar,
    pa.Int64Scalar,
    pa.UInt8Scalar,
    pa.UInt16Scalar,
    pa.UInt32Scalar,
    pa.UInt64Scalar,
    pa.HalfFloatScalar,
    pa.FloatScalar,
    pa.DoubleScalar,
    pa.Decimal128Scalar,
    pa.Decimal256Scalar,
    pa.Date32Scalar,
    pa.Date64Scalar,
    pa.Time32Scalar,
    pa.Time64Scalar,
    pa.TimestampScalar,
    pa.DurationScalar,
    pa.StringScalar,
    pa.BinaryScalar,
    pa.FixedSizeBinaryScalar,
    pa.BinaryViewScalar,
    pa.StringViewScalar,
    pa.ListScalar,
    pa.LargeListScalar,
    pa.ListViewScalar,
    pa.LargeListViewScalar,
    pa.MapScalar,
    pa.FixedSizeListScalar,
    pa.UnionScalar,
    pa.StructScalar,
    pa.DictionaryScalar,
    pa.RunEndEncodedScalar,
    pa.RecordBatchReader,
    pa.ipc.Message,
    pa.ipc.MessageReader,
    pa.MemoryPool,
    pa.LoggingMemoryPool,
    pa.ProxyMemoryPool,
    pa.Device,
    pa.MemoryManager,
    pa.OpaqueArray,
    pa.OpaqueScalar,
    pa.OpaqueType,
    pa.Bool8Array,
    pa.Bool8Scalar,
    pa.Bool8Type,
    pa.JsonArray,
    pa.JsonScalar,
    pa.JsonType,
])
def test_extension_type_constructor_errors(klass):
    # ARROW-2638: prevent calling extension class constructors directly
    msg = f"Do not call {klass.__name__}'s constructor directly, use .* instead."
    with pytest.raises(TypeError, match=msg):
        klass()
