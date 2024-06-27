from __future__ import annotations

import io
import os
import shlex
import subprocess
import sys
import time
from contextlib import contextmanager
from functools import partial

import pytest

s3fs = pytest.importorskip("s3fs")
boto3 = pytest.importorskip("boto3")
moto = pytest.importorskip("moto", minversion="1.3.14")
pytest.importorskip("flask")  # server mode needs flask too
requests = pytest.importorskip("requests")

from fsspec.compression import compr
from fsspec.core import get_fs_token_paths, open_files
from s3fs import S3FileSystem as DaskS3FileSystem
from tlz import concat, valmap

try:
    import fsspec.parquet as fsspec_parquet
except ImportError:
    fsspec_parquet = None

from dask import compute
from dask.bytes.core import read_bytes
from dask.bytes.utils import compress

compute = partial(compute, scheduler="sync")


test_bucket_name = "test"
files = {
    "test/accounts.1.json": (
        b'{"amount": 100, "name": "Alice"}\n'
        b'{"amount": 200, "name": "Bob"}\n'
        b'{"amount": 300, "name": "Charlie"}\n'
        b'{"amount": 400, "name": "Dennis"}\n'
    ),
    "test/accounts.2.json": (
        b'{"amount": 500, "name": "Alice"}\n'
        b'{"amount": 600, "name": "Bob"}\n'
        b'{"amount": 700, "name": "Charlie"}\n'
        b'{"amount": 800, "name": "Dennis"}\n'
    ),
}


@contextmanager
def ensure_safe_environment_variables():
    """
    Get a context manager to safely set environment variables
    All changes will be undone on close, hence environment variables set
    within this contextmanager will neither persist nor change global state.
    """
    saved_environ = dict(os.environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(saved_environ)


@pytest.fixture
def s3so():
    return dict(client_kwargs={"endpoint_url": "http://127.0.0.1:5555/"})


endpoint_uri = "http://127.0.0.1:5555/"


@pytest.fixture(scope="module")
def s3_base():
    with ensure_safe_environment_variables():
        os.environ["AWS_ACCESS_KEY_ID"] = "foobar_key"
        os.environ["AWS_SECRET_ACCESS_KEY"] = "foobar_secret"
        # Ignore any local AWS credentials/config files as they can interfere with moto
        os.environ["AWS_SHARED_CREDENTIALS_FILE"] = ""
        os.environ["AWS_CONFIG_FILE"] = ""

        # pipe to null to avoid logging in terminal
        proc = subprocess.Popen(
            shlex.split("moto_server s3 -p 5555"), stdout=subprocess.DEVNULL
        )

        timeout = 8
        while True:
            try:
                # OK to go once server is accepting connections
                r = requests.get(endpoint_uri)
                if r.ok:
                    break
            except Exception:
                pass
            timeout -= 0.1
            time.sleep(0.1)
            assert timeout > 0, "Timed out waiting for moto server"
        yield

        # shut down external process
        proc.terminate()
        try:
            proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            proc.kill()
            if sys.platform == "win32":
                # belt & braces
                subprocess.call(f"TASKKILL /F /PID {proc.pid} /T")


@pytest.fixture
def s3(s3_base):
    with s3_context() as fs:
        yield fs


@contextmanager
def s3_context(bucket=test_bucket_name, files=files):
    client = boto3.client("s3", endpoint_url=endpoint_uri)
    client.create_bucket(Bucket=bucket, ACL="public-read-write")
    for f, data in files.items():
        client.put_object(Bucket=bucket, Key=f, Body=data)
    fs = s3fs.S3FileSystem(
        anon=True, client_kwargs={"endpoint_url": "http://127.0.0.1:5555/"}
    )
    s3fs.S3FileSystem.clear_instance_cache()
    fs.invalidate_cache()
    try:
        yield fs
    finally:
        fs.rm(bucket, recursive=True)


@pytest.fixture()
def s3_with_yellow_tripdata(s3):
    """
    Fixture with sample yellowtrip CSVs loaded into S3.

    Provides the following CSVs:

    * s3://test/nyc-taxi/2015/yellow_tripdata_2015-01.csv
    * s3://test/nyc-taxi/2014/yellow_tripdata_2015-mm.csv
      for mm from 01 - 12.
    """
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")

    data = {
        "VendorID": {0: 2, 1: 1, 2: 1, 3: 1, 4: 1},
        "tpep_pickup_datetime": {
            0: "2015-01-15 19:05:39",
            1: "2015-01-10 20:33:38",
            2: "2015-01-10 20:33:38",
            3: "2015-01-10 20:33:39",
            4: "2015-01-10 20:33:39",
        },
        "tpep_dropoff_datetime": {
            0: "2015-01-15 19:23:42",
            1: "2015-01-10 20:53:28",
            2: "2015-01-10 20:43:41",
            3: "2015-01-10 20:35:31",
            4: "2015-01-10 20:52:58",
        },
        "passenger_count": {0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
        "trip_distance": {0: 1.59, 1: 3.3, 2: 1.8, 3: 0.5, 4: 3.0},
        "pickup_longitude": {
            0: -73.993896484375,
            1: -74.00164794921875,
            2: -73.96334075927734,
            3: -74.00908660888672,
            4: -73.97117614746094,
        },
        "pickup_latitude": {
            0: 40.7501106262207,
            1: 40.7242431640625,
            2: 40.80278778076172,
            3: 40.71381759643555,
            4: 40.762428283691406,
        },
        "RateCodeID": {0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
        "store_and_fwd_flag": {0: "N", 1: "N", 2: "N", 3: "N", 4: "N"},
        "dropoff_longitude": {
            0: -73.97478485107422,
            1: -73.99441528320312,
            2: -73.95182037353516,
            3: -74.00432586669923,
            4: -74.00418090820312,
        },
        "dropoff_latitude": {
            0: 40.75061798095703,
            1: 40.75910949707031,
            2: 40.82441329956055,
            3: 40.71998596191406,
            4: 40.742652893066406,
        },
        "payment_type": {0: 1, 1: 1, 2: 2, 3: 2, 4: 2},
        "fare_amount": {0: 12.0, 1: 14.5, 2: 9.5, 3: 3.5, 4: 15.0},
        "extra": {0: 1.0, 1: 0.5, 2: 0.5, 3: 0.5, 4: 0.5},
        "mta_tax": {0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5, 4: 0.5},
        "tip_amount": {0: 3.25, 1: 2.0, 2: 0.0, 3: 0.0, 4: 0.0},
        "tolls_amount": {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
        "improvement_surcharge": {0: 0.3, 1: 0.3, 2: 0.3, 3: 0.3, 4: 0.3},
        "total_amount": {0: 17.05, 1: 17.8, 2: 10.8, 3: 4.8, 4: 16.3},
    }
    sample = pd.DataFrame(data)
    df = sample.take(np.arange(5).repeat(10000))
    file = io.BytesIO()
    sfile = io.TextIOWrapper(file)
    df.to_csv(sfile, index=False)

    key = "nyc-taxi/2015/yellow_tripdata_2015-01.csv"
    client = boto3.client("s3", endpoint_url="http://127.0.0.1:5555/")
    client.put_object(Bucket=test_bucket_name, Key=key, Body=file)
    key = "nyc-taxi/2014/yellow_tripdata_2014-{:0>2d}.csv"

    for i in range(1, 13):
        file.seek(0)
        client.put_object(Bucket=test_bucket_name, Key=key.format(i), Body=file)
    yield


def test_get_s3():
    s3 = DaskS3FileSystem(key="key", secret="secret")
    assert s3.key == "key"
    assert s3.secret == "secret"

    s3 = DaskS3FileSystem(username="key", password="secret")
    assert s3.key == "key"
    assert s3.secret == "secret"

    with pytest.raises(KeyError):
        DaskS3FileSystem(key="key", username="key")
    with pytest.raises(KeyError):
        DaskS3FileSystem(secret="key", password="key")


def test_open_files_write(s3, s3so):
    paths = ["s3://" + test_bucket_name + "/more/" + f for f in files]
    fils = open_files(paths, mode="wb", **s3so)
    for fil, data in zip(fils, files.values()):
        with fil as f:
            f.write(data)
    sample, values = read_bytes(
        "s3://" + test_bucket_name + "/more/test/accounts.*", **s3so
    )
    results = compute(*concat(values))
    assert set(list(files.values())) == set(results)


def test_read_bytes(s3, s3so):
    sample, values = read_bytes("s3://" + test_bucket_name + "/test/accounts.*", **s3so)
    assert isinstance(sample, bytes)
    assert sample[:5] == files[sorted(files)[0]][:5]
    assert sample.endswith(b"\n")

    assert isinstance(values, (list, tuple))
    assert isinstance(values[0], (list, tuple))
    assert hasattr(values[0][0], "dask")

    assert sum(map(len, values)) >= len(files)
    results = compute(*concat(values))
    assert set(results) == set(files.values())


def test_read_bytes_sample_delimiter(s3, s3so):
    sample, values = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts.*",
        sample=80,
        delimiter=b"\n",
        **s3so,
    )
    assert sample.endswith(b"\n")
    sample, values = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts.1.json",
        sample=80,
        delimiter=b"\n",
        **s3so,
    )
    assert sample.endswith(b"\n")
    sample, values = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts.1.json",
        sample=2,
        delimiter=b"\n",
        **s3so,
    )
    assert sample.endswith(b"\n")


def test_read_bytes_non_existing_glob(s3, s3so):
    with pytest.raises(IOError):
        read_bytes("s3://" + test_bucket_name + "/non-existing/*", **s3so)


def test_read_bytes_blocksize_none(s3, s3so):
    _, values = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts.*", blocksize=None, **s3so
    )
    assert sum(map(len, values)) == len(files)


def test_read_bytes_blocksize_on_large_data(s3_with_yellow_tripdata, s3so):
    _, L = read_bytes(
        f"s3://{test_bucket_name}/nyc-taxi/2015/yellow_tripdata_2015-01.csv",
        blocksize=None,
        anon=True,
        **s3so,
    )
    assert len(L) == 1

    _, L = read_bytes(
        f"s3://{test_bucket_name}/nyc-taxi/2014/*.csv",
        blocksize=None,
        anon=True,
        **s3so,
    )
    assert len(L) == 12


@pytest.mark.parametrize("blocksize", [5, 15, 45, 1500])
def test_read_bytes_block(s3, blocksize, s3so):
    _, vals = read_bytes(
        "s3://" + test_bucket_name + "/test/account*", blocksize=blocksize, **s3so
    )
    assert list(map(len, vals)) == [
        max((len(v) // blocksize), 1) for v in files.values()
    ]

    results = compute(*concat(vals))
    assert sum(len(r) for r in results) == sum(len(v) for v in files.values())

    ourlines = b"".join(results).split(b"\n")
    testlines = b"".join(files.values()).split(b"\n")
    assert set(ourlines) == set(testlines)


@pytest.mark.parametrize("blocksize", [5, 15, 45, 1500])
def test_read_bytes_delimited(s3, blocksize, s3so):
    _, values = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts*",
        blocksize=blocksize,
        delimiter=b"\n",
        **s3so,
    )
    _, values2 = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts*",
        blocksize=blocksize,
        delimiter=b"foo",
        **s3so,
    )
    assert [a.key for a in concat(values)] != [b.key for b in concat(values2)]

    results = compute(*concat(values))
    res = [r for r in results if r]
    assert all(r.endswith(b"\n") for r in res)
    ourlines = b"".join(res).split(b"\n")
    testlines = b"".join(files[k] for k in sorted(files)).split(b"\n")
    assert ourlines == testlines

    # delimiter not at the end
    d = b"}"
    _, values = read_bytes(
        "s3://" + test_bucket_name + "/test/accounts*",
        blocksize=blocksize,
        delimiter=d,
        **s3so,
    )
    results = compute(*concat(values))
    res = [r for r in results if r]
    # All should end in } except EOF
    assert sum(r.endswith(b"}") for r in res) == len(res) - 2
    ours = b"".join(res)
    test = b"".join(files[v] for v in sorted(files))
    assert ours == test


@pytest.mark.parametrize(
    "fmt,blocksize",
    [(fmt, None) for fmt in compr] + [(fmt, 10) for fmt in compr],
)
def test_compression(s3, fmt, blocksize, s3so):
    if fmt not in compress:
        pytest.skip("compression function not provided")
    s3._cache.clear()
    with s3_context("compress", valmap(compress[fmt], files)):
        if fmt and blocksize:
            with pytest.raises(ValueError):
                read_bytes(
                    "s3://compress/test/accounts.*",
                    compression=fmt,
                    blocksize=blocksize,
                    **s3so,
                )
            return
        sample, values = read_bytes(
            "s3://compress/test/accounts.*",
            compression=fmt,
            blocksize=blocksize,
            **s3so,
        )
        assert sample.startswith(files[sorted(files)[0]][:10])
        assert sample.endswith(b"\n")

        results = compute(*concat(values))
        assert b"".join(results) == b"".join([files[k] for k in sorted(files)])


@pytest.mark.parametrize("mode", ["rt", "rb"])
def test_open_files(s3, mode, s3so):
    myfiles = open_files(
        "s3://" + test_bucket_name + "/test/accounts.*", mode=mode, **s3so
    )
    assert len(myfiles) == len(files)
    for lazy_file, path in zip(myfiles, sorted(files)):
        with lazy_file as f:
            data = f.read()
            sol = files[path]
            assert data == sol if mode == "rb" else sol.decode()


double = lambda x: x * 2


def test_modification_time_read_bytes(s3, s3so):
    with s3_context("compress", files):
        _, a = read_bytes("s3://compress/test/accounts.*", anon=True, **s3so)
        _, b = read_bytes("s3://compress/test/accounts.*", anon=True, **s3so)

        assert [aa._key for aa in concat(a)] == [bb._key for bb in concat(b)]

    with s3_context("compress", valmap(double, files)):
        _, c = read_bytes("s3://compress/test/accounts.*", anon=True, **s3so)

    assert [aa._key for aa in concat(a)] != [cc._key for cc in concat(c)]


@pytest.fixture(
    params=[
        "pyarrow",
        pytest.param(
            "fastparquet", marks=pytest.mark.filterwarnings("ignore::FutureWarning")
        ),
    ]
)
def engine(request):
    import dask.dataframe as dd
    from dask.array.numpy_compat import NUMPY_GE_200

    if NUMPY_GE_200 and request.param == "fastparquet":
        # https://github.com/dask/fastparquet/issues/923
        pytest.skip("fastparquet doesn't work with Numpy 2")

    pytest.importorskip(request.param)

    if dd._dask_expr_enabled() and request.param == "fastparquet":
        pytest.skip("not supported")
    return request.param


@pytest.mark.parametrize("metadata_file", [True, False])
def test_parquet(s3, engine, s3so, metadata_file):
    dd = pytest.importorskip("dask.dataframe")
    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")

    url = "s3://%s/test.parquet" % test_bucket_name
    data = pd.DataFrame(
        {
            "i32": np.arange(1000, dtype=np.int32),
            "i64": np.arange(1000, dtype=np.int64),
            "f": np.arange(1000, dtype=np.float64),
            "bhello": np.random.choice(["hello", "you", "people"], size=1000).astype(
                "O"
            ),
        },
        index=pd.Index(np.arange(1000), name="foo"),
    )
    df = dd.from_pandas(data, chunksize=500)
    df.to_parquet(
        url, engine=engine, storage_options=s3so, write_metadata_file=metadata_file
    )

    files = [f.split("/")[-1] for f in s3.ls(url)]
    if metadata_file:
        assert "_common_metadata" in files
        assert "_metadata" in files
    assert "part.0.parquet" in files

    df2 = dd.read_parquet(
        url, index="foo", calculate_divisions=True, engine=engine, storage_options=s3so
    )
    assert len(df2.divisions) > 1

    dd.utils.assert_eq(data, df2)

    # Check that `open_file_options` arguments are
    # really passed through to fsspec
    if fsspec_parquet:
        # Passing `open_file_options` kwargs will fail
        # if you set an unsupported engine
        with pytest.raises(ValueError):
            dd.read_parquet(
                url,
                engine=engine,
                storage_options=s3so,
                open_file_options={
                    "precache_options": {"method": "parquet", "engine": "foo"},
                },
            ).compute()

        # ...but should work fine if you modify the
        # maximum block-transfer size (max_block)
        dd.read_parquet(
            url,
            engine=engine,
            storage_options=s3so,
            open_file_options={
                "precache_options": {"method": "parquet", "max_block": 8_000},
            },
        ).compute()

    # Check "open_file_func"
    fs = get_fs_token_paths(url, storage_options=s3so)[0]

    def _open(*args, check=True, **kwargs):
        assert check
        return fs.open(*args, **kwargs)

    # Should fail if `check=False`
    with pytest.raises(AssertionError):
        dd.read_parquet(
            url,
            engine=engine,
            storage_options=s3so,
            open_file_options={"open_file_func": _open, "check": False},
        ).compute()

    # Should succeed otherwise
    df3 = dd.read_parquet(
        url,
        engine=engine,
        storage_options=s3so,
        open_file_options={"open_file_func": _open},
    )
    dd.utils.assert_eq(data, df3)

    # Check that `cache_type="all"` result is same
    df4 = dd.read_parquet(
        url,
        engine=engine,
        storage_options=s3so,
        open_file_options={"cache_type": "all"},
    )
    dd.utils.assert_eq(data, df4)


def test_parquet_append(s3, engine, s3so):
    dd = pytest.importorskip("dask.dataframe")
    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")

    url = "s3://%s/test.parquet.append" % test_bucket_name

    data = pd.DataFrame(
        {
            "i32": np.arange(1000, dtype=np.int32),
            "i64": np.arange(1000, dtype=np.int64),
            "f": np.arange(1000, dtype=np.float64),
            "bhello": np.random.choice(["hello", "you", "people"], size=1000).astype(
                "O"
            ),
        },
    )
    df = dd.from_pandas(data, chunksize=500)
    df.to_parquet(
        url,
        engine=engine,
        storage_options=s3so,
        write_index=False,
        write_metadata_file=True,
    )
    df.to_parquet(
        url,
        engine=engine,
        storage_options=s3so,
        write_index=False,
        append=True,
        ignore_divisions=True,
    )

    files = [f.split("/")[-1] for f in s3.ls(url)]
    assert "_common_metadata" in files
    assert "_metadata" in files
    assert "part.0.parquet" in files

    df2 = dd.read_parquet(
        url,
        index=False,
        engine=engine,
        storage_options=s3so,
    )

    dd.utils.assert_eq(
        pd.concat([data, data]),
        df2,
        check_index=False,
    )


def test_parquet_wstoragepars(s3, s3so, engine):
    dd = pytest.importorskip("dask.dataframe")
    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")

    url = "s3://%s/test.parquet" % test_bucket_name

    data = pd.DataFrame({"i32": np.array([0, 5, 2, 5])})
    df = dd.from_pandas(data, chunksize=500)
    df.to_parquet(
        url,
        engine=engine,
        write_index=False,
        storage_options=s3so,
        write_metadata_file=True,
    )

    dd.read_parquet(
        url,
        engine=engine,
        storage_options={"default_fill_cache": False, **s3so},
    )
    assert s3.current().default_fill_cache is False
    dd.read_parquet(
        url,
        engine=engine,
        storage_options={"default_fill_cache": True, **s3so},
    )
    assert s3.current().default_fill_cache is True

    dd.read_parquet(
        url,
        engine=engine,
        storage_options={"default_block_size": 2**20, **s3so},
    )
    assert s3.current().default_block_size == 2**20
    with s3.current().open(url + "/_metadata") as f:
        assert f.blocksize == 2**20
