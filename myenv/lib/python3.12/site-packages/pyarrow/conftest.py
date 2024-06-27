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

import pytest
import pyarrow as pa
from pyarrow import Codec
from pyarrow import fs

import numpy as np

groups = [
    'acero',
    'azure',
    'brotli',
    'bz2',
    'cython',
    'dataset',
    'hypothesis',
    'fastparquet',
    'gandiva',
    'gcs',
    'gdb',
    'gzip',
    'hdfs',
    'large_memory',
    'lz4',
    'memory_leak',
    'nopandas',
    'orc',
    'pandas',
    'parquet',
    'parquet_encryption',
    's3',
    'snappy',
    'substrait',
    'flight',
    'slow',
    'requires_testing_data',
    'zstd',
]

defaults = {
    'acero': False,
    'azure': False,
    'brotli': Codec.is_available('brotli'),
    'bz2': Codec.is_available('bz2'),
    'cython': False,
    'dataset': False,
    'fastparquet': False,
    'flight': False,
    'gandiva': False,
    'gcs': False,
    'gdb': True,
    'gzip': Codec.is_available('gzip'),
    'hdfs': False,
    'hypothesis': False,
    'large_memory': False,
    'lz4': Codec.is_available('lz4'),
    'memory_leak': False,
    'nopandas': False,
    'orc': False,
    'pandas': False,
    'parquet': False,
    'parquet_encryption': False,
    'requires_testing_data': True,
    's3': False,
    'slow': False,
    'snappy': Codec.is_available('snappy'),
    'substrait': False,
    'zstd': Codec.is_available('zstd'),
}

try:
    import cython  # noqa
    defaults['cython'] = True
except ImportError:
    pass

try:
    import fastparquet  # noqa
    defaults['fastparquet'] = True
except ImportError:
    pass

try:
    import pyarrow.gandiva  # noqa
    defaults['gandiva'] = True
except ImportError:
    pass

try:
    import pyarrow.acero  # noqa
    defaults['acero'] = True
except ImportError:
    pass

try:
    import pyarrow.dataset  # noqa
    defaults['dataset'] = True
except ImportError:
    pass

try:
    import pyarrow.orc  # noqa
    defaults['orc'] = True
except ImportError:
    pass

try:
    import pandas  # noqa
    defaults['pandas'] = True
except ImportError:
    defaults['nopandas'] = True

try:
    import pyarrow.parquet  # noqa
    defaults['parquet'] = True
except ImportError:
    pass

try:
    import pyarrow.parquet.encryption  # noqa
    defaults['parquet_encryption'] = True
except ImportError:
    pass

try:
    import pyarrow.flight  # noqa
    defaults['flight'] = True
except ImportError:
    pass

try:
    from pyarrow.fs import AzureFileSystem  # noqa
    defaults['azure'] = True
except ImportError:
    pass

try:
    from pyarrow.fs import GcsFileSystem  # noqa
    defaults['gcs'] = True
except ImportError:
    pass

try:
    from pyarrow.fs import S3FileSystem  # noqa
    defaults['s3'] = True
except ImportError:
    pass

try:
    from pyarrow.fs import HadoopFileSystem  # noqa
    defaults['hdfs'] = True
except ImportError:
    pass

try:
    import pyarrow.substrait  # noqa
    defaults['substrait'] = True
except ImportError:
    pass


# Doctest should ignore files for the modules that are not built
def pytest_ignore_collect(path, config):
    if config.option.doctestmodules:
        # don't try to run doctests on the /tests directory
        if "/pyarrow/tests/" in str(path):
            return True

        doctest_groups = [
            'dataset',
            'orc',
            'parquet',
            'flight',
            'substrait',
        ]

        # handle cuda, flight, etc
        for group in doctest_groups:
            if 'pyarrow/{}'.format(group) in str(path):
                if not defaults[group]:
                    return True

        if 'pyarrow/parquet/encryption' in str(path):
            if not defaults['parquet_encryption']:
                return True

        if 'pyarrow/cuda' in str(path):
            try:
                import pyarrow.cuda  # noqa
                return False
            except ImportError:
                return True

        if 'pyarrow/fs' in str(path):
            try:
                from pyarrow.fs import S3FileSystem  # noqa
                return False
            except ImportError:
                return True

    if getattr(config.option, "doctest_cython", False):
        if "/pyarrow/tests/" in str(path):
            return True
        if "/pyarrow/_parquet_encryption" in str(path):
            return True

    return False


# Save output files from doctest examples into temp dir
@pytest.fixture(autouse=True)
def _docdir(request):

    # Trigger ONLY for the doctests
    doctest_m = request.config.option.doctestmodules
    doctest_c = getattr(request.config.option, "doctest_cython", False)

    if doctest_m or doctest_c:

        # Get the fixture dynamically by its name.
        tmpdir = request.getfixturevalue('tmpdir')

        # Chdir only for the duration of the test.
        with tmpdir.as_cwd():
            yield

    else:
        yield


# Define doctest_namespace for fs module docstring import
@pytest.fixture(autouse=True)
def add_fs(doctest_namespace, request, tmp_path):

    # Trigger ONLY for the doctests
    doctest_m = request.config.option.doctestmodules
    doctest_c = getattr(request.config.option, "doctest_cython", False)

    if doctest_m or doctest_c:
        # fs import
        doctest_namespace["fs"] = fs

        # Creation of an object and file with data
        local = fs.LocalFileSystem()
        path = tmp_path / 'pyarrow-fs-example.dat'
        with local.open_output_stream(str(path)) as stream:
            stream.write(b'data')
        doctest_namespace["local"] = local
        doctest_namespace["local_path"] = str(tmp_path)
        doctest_namespace["path"] = str(path)
    yield


# Define udf fixture for test_udf.py and test_substrait.py
@pytest.fixture(scope="session")
def unary_func_fixture():
    """
    Register a unary scalar function.
    """
    from pyarrow import compute as pc

    def unary_function(ctx, x):
        return pc.call_function("add", [x, 1],
                                memory_pool=ctx.memory_pool)
    func_name = "y=x+1"
    unary_doc = {"summary": "add function",
                 "description": "test add function"}
    pc.register_scalar_function(unary_function,
                                func_name,
                                unary_doc,
                                {"array": pa.int64()},
                                pa.int64())
    return unary_function, func_name


@pytest.fixture(scope="session")
def unary_agg_func_fixture():
    """
    Register a unary aggregate function (mean)
    """
    from pyarrow import compute as pc

    def func(ctx, x):
        return pa.scalar(np.nanmean(x))

    func_name = "mean_udf"
    func_doc = {"summary": "y=avg(x)",
                "description": "find mean of x"}

    pc.register_aggregate_function(func,
                                   func_name,
                                   func_doc,
                                   {
                                       "x": pa.float64(),
                                   },
                                   pa.float64()
                                   )
    return func, func_name


@pytest.fixture(scope="session")
def varargs_agg_func_fixture():
    """
    Register a unary aggregate function
    """
    from pyarrow import compute as pc

    def func(ctx, *args):
        sum = 0.0
        for arg in args:
            sum += np.nanmean(arg)
        return pa.scalar(sum)

    func_name = "sum_mean"
    func_doc = {"summary": "Varargs aggregate",
                "description": "Varargs aggregate"}

    pc.register_aggregate_function(func,
                                   func_name,
                                   func_doc,
                                   {
                                       "x": pa.int64(),
                                       "y": pa.float64()
                                   },
                                   pa.float64()
                                   )
    return func, func_name
