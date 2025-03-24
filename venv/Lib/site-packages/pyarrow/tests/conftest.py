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

import functools
import os
import pathlib
import subprocess
import sys
import time
import urllib.request

import pytest
import hypothesis as h

from ..conftest import groups, defaults

from pyarrow import set_timezone_db_path
from pyarrow.util import find_free_port


# setup hypothesis profiles
h.settings.register_profile('ci', max_examples=1000)
h.settings.register_profile('dev', max_examples=50)
h.settings.register_profile('debug', max_examples=10,
                            verbosity=h.Verbosity.verbose)

# load default hypothesis profile, either set HYPOTHESIS_PROFILE environment
# variable or pass --hypothesis-profile option to pytest, to see the generated
# examples try:
# pytest pyarrow -sv --enable-hypothesis --hypothesis-profile=debug
h.settings.load_profile(os.environ.get('HYPOTHESIS_PROFILE', 'dev'))

# Set this at the beginning before the AWS SDK was loaded to avoid reading in
# user configuration values.
os.environ['AWS_CONFIG_FILE'] = "/dev/null"


if sys.platform == 'win32':
    tzdata_set_path = os.environ.get('PYARROW_TZDATA_PATH', None)
    if tzdata_set_path:
        set_timezone_db_path(tzdata_set_path)


# GH-45295: For ORC, try to populate TZDIR env var from tzdata package resource
# path.
#
# Note this is a different kind of database than what we allow to be set by
# `PYARROW_TZDATA_PATH` and passed to set_timezone_db_path.
if sys.platform == 'win32':
    if os.environ.get('TZDIR', None) is None:
        from importlib import resources
        try:
            os.environ['TZDIR'] = os.path.join(resources.files('tzdata'), 'zoneinfo')
        except ModuleNotFoundError:
            print(
                'Package "tzdata" not found. Not setting TZDIR environment variable.'
            )


def pytest_addoption(parser):
    # Create options to selectively enable test groups
    def bool_env(name, default=None):
        value = os.environ.get(name.upper())
        if not value:  # missing or empty
            return default
        value = value.lower()
        if value in {'1', 'true', 'on', 'yes', 'y'}:
            return True
        elif value in {'0', 'false', 'off', 'no', 'n'}:
            return False
        else:
            raise ValueError('{}={} is not parsable as boolean'
                             .format(name.upper(), value))

    for group in groups:
        default = bool_env('PYARROW_TEST_{}'.format(group), defaults[group])
        parser.addoption('--enable-{}'.format(group),
                         action='store_true', default=default,
                         help=('Enable the {} test group'.format(group)))
        parser.addoption('--disable-{}'.format(group),
                         action='store_true', default=False,
                         help=('Disable the {} test group'.format(group)))


class PyArrowConfig:
    def __init__(self):
        self.is_enabled = {}

    def apply_mark(self, mark):
        group = mark.name
        if group in groups:
            self.requires(group)

    def requires(self, group):
        if not self.is_enabled[group]:
            pytest.skip('{} NOT enabled'.format(group))


def pytest_configure(config):
    # Apply command-line options to initialize PyArrow-specific config object
    config.pyarrow = PyArrowConfig()

    for mark in groups:
        config.addinivalue_line(
            "markers", mark,
        )

        enable_flag = '--enable-{}'.format(mark)
        disable_flag = '--disable-{}'.format(mark)

        is_enabled = (config.getoption(enable_flag) and not
                      config.getoption(disable_flag))
        config.pyarrow.is_enabled[mark] = is_enabled


def pytest_runtest_setup(item):
    # Apply test markers to skip tests selectively
    for mark in item.iter_markers():
        item.config.pyarrow.apply_mark(mark)


@pytest.fixture
def tempdir(tmpdir):
    # convert pytest's LocalPath to pathlib.Path
    return pathlib.Path(tmpdir.strpath)


@pytest.fixture(scope='session')
def base_datadir():
    return pathlib.Path(__file__).parent / 'data'


@pytest.fixture(autouse=True)
def disable_aws_metadata(monkeypatch):
    """Stop the AWS SDK from trying to contact the EC2 metadata server.

    Otherwise, this causes a 5 second delay in tests that exercise the
    S3 filesystem.
    """
    monkeypatch.setenv("AWS_EC2_METADATA_DISABLED", "true")


# TODO(kszucs): move the following fixtures to test_fs.py once the previous
# parquet dataset implementation and hdfs implementation are removed.

@pytest.fixture(scope='session')
def hdfs_connection():
    host = os.environ.get('ARROW_HDFS_TEST_HOST', 'default')
    port = int(os.environ.get('ARROW_HDFS_TEST_PORT', 0))
    user = os.environ.get('ARROW_HDFS_TEST_USER', 'hdfs')
    return host, port, user


@pytest.fixture(scope='session')
def s3_connection():
    host, port = '127.0.0.1', find_free_port()
    access_key, secret_key = 'arrow', 'apachearrow'
    return host, port, access_key, secret_key


def retry(attempts=3, delay=1.0, max_delay=None, backoff=1):
    """
    Retry decorator

    Parameters
    ----------
    attempts : int, default 3
        The number of attempts.
    delay : float, default 1
        Initial delay in seconds.
    max_delay : float, optional
        The max delay between attempts.
    backoff : float, default 1
        The multiplier to delay after each attempt.
    """
    def decorate(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            remaining_attempts = attempts
            curr_delay = delay
            while remaining_attempts > 0:
                try:
                    return func(*args, **kwargs)
                except Exception as err:
                    remaining_attempts -= 1
                    last_exception = err
                    curr_delay *= backoff
                    if max_delay:
                        curr_delay = min(curr_delay, max_delay)
                    time.sleep(curr_delay)
            raise last_exception
        return wrapper
    return decorate


@pytest.fixture(scope='session')
def s3_server(s3_connection, tmpdir_factory):
    @retry(attempts=5, delay=1, backoff=2)
    def minio_server_health_check(address):
        resp = urllib.request.urlopen(f"http://{address}/minio/health/live")
        assert resp.getcode() == 200

    tmpdir = tmpdir_factory.getbasetemp()
    host, port, access_key, secret_key = s3_connection

    address = '{}:{}'.format(host, port)
    env = os.environ.copy()
    env.update({
        'MINIO_ACCESS_KEY': access_key,
        'MINIO_SECRET_KEY': secret_key
    })

    args = ['minio', '--compat', 'server', '--quiet', '--address',
            address, tmpdir]
    proc = None
    try:
        proc = subprocess.Popen(args, env=env)
    except OSError:
        pytest.skip('`minio` command cannot be located')
    else:
        # Wait for the server to startup before yielding
        minio_server_health_check(address)

        yield {
            'connection': s3_connection,
            'process': proc,
            'tempdir': tmpdir
        }
    finally:
        if proc is not None:
            proc.kill()
            proc.wait()


@pytest.fixture(scope='session')
def gcs_server():
    port = find_free_port()
    env = os.environ.copy()
    exe = 'storage-testbench'
    args = [exe, '--port', str(port)]
    proc = None
    try:
        # start server
        proc = subprocess.Popen(args, env=env)
        # Make sure the server is alive.
        if proc.poll() is not None:
            pytest.skip(f"Command {args} did not start server successfully!")
    except OSError as e:
        pytest.skip(f"Command {args} failed to execute: {e}")
    else:
        yield {
            'connection': ('localhost', port),
            'process': proc,
        }
    finally:
        if proc is not None:
            proc.kill()
            proc.wait()


@pytest.fixture(scope='session')
def azure_server(tmpdir_factory):
    port = find_free_port()
    env = os.environ.copy()
    tmpdir = tmpdir_factory.getbasetemp()
    # We only need blob service emulator, not queue or table.
    args = ['azurite-blob', "--location", tmpdir, "--blobPort", str(port)]
    # For old Azurite. We can't install the latest Azurite with old
    # Node.js on old Ubuntu.
    args += ["--skipApiVersionCheck"]
    proc = None
    try:
        proc = subprocess.Popen(args, env=env)
        # Make sure the server is alive.
        if proc.poll() is not None:
            pytest.skip(f"Command {args} did not start server successfully!")
    except (ModuleNotFoundError, OSError) as e:
        pytest.skip(f"Command {args} failed to execute: {e}")
    else:
        yield {
            # Use the standard azurite account_name and account_key.
            # https://learn.microsoft.com/en-us/azure/storage/common/storage-use-emulator#authorize-with-shared-key-credentials
            'connection': ('127.0.0.1', port, 'devstoreaccount1',
                           'Eby8vdM02xNOcqFlqUwJPLlmEtlCDXJ1OUzFT50uSRZ6IFsuFq2'
                           'UVErCz4I6tq/K1SZFPTOtr/KBHBeksoGMGw=='),
            'process': proc,
            'tempdir': tmpdir,
        }
    finally:
        if proc is not None:
            proc.kill()
            proc.wait()


@pytest.fixture(
    params=[
        'builtin_pickle',
        'cloudpickle'
    ],
    scope='session'
)
def pickle_module(request):
    return request.getfixturevalue(request.param)


@pytest.fixture(scope='session')
def builtin_pickle():
    import pickle
    return pickle


@pytest.fixture(scope='session')
def cloudpickle():
    cp = pytest.importorskip('cloudpickle')
    if 'HIGHEST_PROTOCOL' not in cp.__dict__:
        cp.HIGHEST_PROTOCOL = cp.DEFAULT_PROTOCOL
    return cp
