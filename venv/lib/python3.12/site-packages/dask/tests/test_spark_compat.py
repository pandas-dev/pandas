from __future__ import annotations

import decimal
import signal
import sys
import threading

import pytest

from dask.datasets import timeseries

pytest.importorskip("pandas")
pyspark = pytest.importorskip("pyspark")
pa = pytest.importorskip("pyarrow")
dd = pytest.importorskip("dask.dataframe")

import numpy as np
import pandas as pd

from dask.dataframe.utils import assert_eq

pytestmark = [
    pytest.mark.skipif(
        sys.platform != "linux",
        reason="Unnecessary, and hard to get spark working on non-linux platforms",
    ),
    pytest.mark.skip(
        reason="pyspark doesn't yet have support for pandas 2.0",
    ),
    # we only test with pyarrow strings and pandas 2.0
    pytest.mark.skip_with_pyarrow_strings,  # pyspark doesn't support pandas 2.0
]


@pytest.fixture(
    params=[
        "pyarrow",
    ]
)
def engine(request):
    pytest.importorskip(request.param)
    return request.param


# pyspark auto-converts timezones -- round-tripping timestamps is easier if
# we set everything to UTC.
pdf = timeseries(freq="1h").compute()
pdf.index = pdf.index.tz_localize("UTC")
pdf = pdf.reset_index()


@pytest.fixture(scope="module")
def spark_session():
    # Spark registers a global signal handler that can cause problems elsewhere
    # in the test suite. In particular, the handler fails if the spark session
    # is stopped (a bug in pyspark).
    prev = signal.getsignal(signal.SIGINT)
    # Create a spark session. Note that we set the timezone to UTC to avoid
    # conversion to local time when reading parquet files.
    spark = (
        pyspark.sql.SparkSession.builder.master("local")
        .appName("Dask Testing")
        .config("spark.sql.session.timeZone", "UTC")
        .getOrCreate()
    )
    yield spark

    spark.stop()
    # Make sure we get rid of the signal once we leave stop the session.
    if threading.current_thread() is threading.main_thread():
        signal.signal(signal.SIGINT, prev)


@pytest.mark.parametrize("npartitions", [1, 5, 10])
def test_roundtrip_parquet_spark_to_dask(spark_session, npartitions, tmpdir, engine):
    tmpdir = str(tmpdir)

    sdf = spark_session.createDataFrame(pdf)
    # We are not overwriting any data, but spark complains if the directory
    # already exists (as tmpdir does) and we don't set overwrite
    sdf.repartition(npartitions).write.parquet(tmpdir, mode="overwrite")

    ddf = dd.read_parquet(tmpdir, engine=engine)
    # Papercut: pandas TZ localization doesn't survive roundtrip
    ddf = ddf.assign(timestamp=ddf.timestamp.dt.tz_localize("UTC"))
    assert ddf.npartitions == npartitions

    assert_eq(ddf, pdf, check_index=False)


def test_roundtrip_hive_parquet_spark_to_dask(spark_session, tmpdir, engine):
    tmpdir = str(tmpdir)

    sdf = spark_session.createDataFrame(pdf)
    # not overwriting any data, but spark complains if the directory
    # already exists and we don't set overwrite
    sdf.write.parquet(tmpdir, mode="overwrite", partitionBy="name")

    ddf = dd.read_parquet(tmpdir, engine=engine)
    # Papercut: pandas TZ localization doesn't survive roundtrip
    ddf = ddf.assign(timestamp=ddf.timestamp.dt.tz_localize("UTC"))

    # Partitioning can change the column order. This is mostly okay,
    # but we sort them here to ease comparison
    ddf = ddf.compute().sort_index(axis=1)
    # Dask automatically converts hive-partitioned columns to categories.
    # This is fine, but convert back to strings for comparison.
    ddf = ddf.assign(name=ddf.name.astype("str"))

    assert_eq(ddf, pdf.sort_index(axis=1), check_index=False)


@pytest.mark.parametrize("npartitions", [1, 5, 10])
def test_roundtrip_parquet_dask_to_spark(spark_session, npartitions, tmpdir, engine):
    tmpdir = str(tmpdir)
    ddf = dd.from_pandas(pdf, npartitions=npartitions)

    ddf.to_parquet(tmpdir, engine=engine, write_index=False)

    sdf = spark_session.read.parquet(tmpdir)
    sdf = sdf.toPandas()

    # Papercut: pandas TZ localization doesn't survive roundtrip
    sdf = sdf.assign(timestamp=sdf.timestamp.dt.tz_localize("UTC"))

    assert_eq(sdf, ddf, check_index=False)


def test_roundtrip_parquet_spark_to_dask_extension_dtypes(spark_session, tmpdir):
    tmpdir = str(tmpdir)
    npartitions = 5

    size = 20
    pdf = pd.DataFrame(
        {
            "a": range(size),
            "b": np.random.random(size=size),
            "c": [True, False] * (size // 2),
            "d": ["alice", "bob"] * (size // 2),
        }
    )
    # Note: since we set dtype_backend="numpy_nullable" below, we are expecting *all*
    # of the resulting series to use those dtypes. If there is a mix of nullable
    # and non-nullable dtypes here, then that will result in dtype mismatches
    # in the finale frame.
    pdf = pdf.astype(
        {
            "a": "Int64",
            "b": "Float64",
            "c": "boolean",
            "d": "string",
        }
    )
    # # Ensure all columns are extension dtypes
    assert all([pd.api.types.is_extension_array_dtype(dtype) for dtype in pdf.dtypes])

    sdf = spark_session.createDataFrame(pdf)
    # We are not overwriting any data, but spark complains if the directory
    # already exists (as tmpdir does) and we don't set overwrite
    sdf.repartition(npartitions).write.parquet(tmpdir, mode="overwrite")

    ddf = dd.read_parquet(tmpdir, engine="pyarrow", dtype_backend="numpy_nullable")
    assert all(
        [pd.api.types.is_extension_array_dtype(dtype) for dtype in ddf.dtypes]
    ), ddf.dtypes
    assert_eq(ddf, pdf, check_index=False)


def test_read_decimal_dtype_pyarrow(spark_session, tmpdir):
    tmpdir = str(tmpdir)
    npartitions = 3
    size = 6

    decimal_data = [
        decimal.Decimal("8093.234"),
        decimal.Decimal("8094.234"),
        decimal.Decimal("8095.234"),
        decimal.Decimal("8096.234"),
        decimal.Decimal("8097.234"),
        decimal.Decimal("8098.234"),
    ]
    pdf = pd.DataFrame(
        {
            "a": range(size),
            "b": decimal_data,
        }
    )
    sdf = spark_session.createDataFrame(pdf)
    sdf = sdf.withColumn("b", sdf["b"].cast(pyspark.sql.types.DecimalType(7, 3)))
    # We are not overwriting any data, but spark complains if the directory
    # already exists (as tmpdir does) and we don't set overwrite
    sdf.repartition(npartitions).write.parquet(tmpdir, mode="overwrite")

    ddf = dd.read_parquet(tmpdir, engine="pyarrow", dtype_backend="pyarrow")
    assert ddf.b.dtype.pyarrow_dtype == pa.decimal128(7, 3)
    assert ddf.b.compute().dtype.pyarrow_dtype == pa.decimal128(7, 3)
    expected = pdf.astype(
        {
            "a": "int64[pyarrow]",
            "b": pd.ArrowDtype(pa.decimal128(7, 3)),
        }
    )

    assert_eq(ddf, expected, check_index=False)
