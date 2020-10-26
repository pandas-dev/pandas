""" test parquet compat """
import datetime
from distutils.version import LooseVersion
from io import BytesIO
import os
from warnings import catch_warnings

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm

from pandas.io.parquet import (
    FastParquetImpl,
    PyArrowImpl,
    get_engine,
    read_parquet,
    to_parquet,
)

try:
    import pyarrow  # noqa

    _HAVE_PYARROW = True
except ImportError:
    _HAVE_PYARROW = False

try:
    import fastparquet  # noqa

    _HAVE_FASTPARQUET = True
except ImportError:
    _HAVE_FASTPARQUET = False


pytestmark = pytest.mark.filterwarnings(
    "ignore:RangeIndex.* is deprecated:DeprecationWarning"
)


# setup engines & skips
@pytest.fixture(
    params=[
        pytest.param(
            "fastparquet",
            marks=pytest.mark.skipif(
                not _HAVE_FASTPARQUET, reason="fastparquet is not installed"
            ),
        ),
        pytest.param(
            "pyarrow",
            marks=pytest.mark.skipif(
                not _HAVE_PYARROW, reason="pyarrow is not installed"
            ),
        ),
    ]
)
def engine(request):
    return request.param


@pytest.fixture
def pa():
    if not _HAVE_PYARROW:
        pytest.skip("pyarrow is not installed")
    return "pyarrow"


@pytest.fixture
def fp():
    if not _HAVE_FASTPARQUET:
        pytest.skip("fastparquet is not installed")
    return "fastparquet"


@pytest.fixture
def df_compat():
    return pd.DataFrame({"A": [1, 2, 3], "B": "foo"})


@pytest.fixture
def df_cross_compat():
    df = pd.DataFrame(
        {
            "a": list("abc"),
            "b": list(range(1, 4)),
            # 'c': np.arange(3, 6).astype('u1'),
            "d": np.arange(4.0, 7.0, dtype="float64"),
            "e": [True, False, True],
            "f": pd.date_range("20130101", periods=3),
            # 'g': pd.date_range('20130101', periods=3,
            #                    tz='US/Eastern'),
            # 'h': pd.date_range('20130101', periods=3, freq='ns')
        }
    )
    return df


@pytest.fixture
def df_full():
    return pd.DataFrame(
        {
            "string": list("abc"),
            "string_with_nan": ["a", np.nan, "c"],
            "string_with_none": ["a", None, "c"],
            "bytes": [b"foo", b"bar", b"baz"],
            "unicode": ["foo", "bar", "baz"],
            "int": list(range(1, 4)),
            "uint": np.arange(3, 6).astype("u1"),
            "float": np.arange(4.0, 7.0, dtype="float64"),
            "float_with_nan": [2.0, np.nan, 3.0],
            "bool": [True, False, True],
            "datetime": pd.date_range("20130101", periods=3),
            "datetime_with_nat": [
                pd.Timestamp("20130101"),
                pd.NaT,
                pd.Timestamp("20130103"),
            ],
        }
    )


def check_round_trip(
    df,
    engine=None,
    path=None,
    write_kwargs=None,
    read_kwargs=None,
    expected=None,
    check_names=True,
    check_like=False,
    repeat=2,
):
    """Verify parquet serializer and deserializer produce the same results.

    Performs a pandas to disk and disk to pandas round trip,
    then compares the 2 resulting DataFrames to verify equality.

    Parameters
    ----------
    df: Dataframe
    engine: str, optional
        'pyarrow' or 'fastparquet'
    path: str, optional
    write_kwargs: dict of str:str, optional
    read_kwargs: dict of str:str, optional
    expected: DataFrame, optional
        Expected deserialization result, otherwise will be equal to `df`
    check_names: list of str, optional
        Closed set of column names to be compared
    check_like: bool, optional
        If True, ignore the order of index & columns.
    repeat: int, optional
        How many times to repeat the test
    """
    write_kwargs = write_kwargs or {"compression": None}
    read_kwargs = read_kwargs or {}

    if expected is None:
        expected = df

    if engine:
        write_kwargs["engine"] = engine
        read_kwargs["engine"] = engine

    def compare(repeat):
        for _ in range(repeat):
            df.to_parquet(path, **write_kwargs)
            with catch_warnings(record=True):
                actual = read_parquet(path, **read_kwargs)

            tm.assert_frame_equal(
                expected, actual, check_names=check_names, check_like=check_like
            )

    if path is None:
        with tm.ensure_clean() as path:
            compare(repeat)
    else:
        compare(repeat)


def test_invalid_engine(df_compat):
    with pytest.raises(ValueError):
        check_round_trip(df_compat, "foo", "bar")


def test_options_py(df_compat, pa):
    # use the set option

    with pd.option_context("io.parquet.engine", "pyarrow"):
        check_round_trip(df_compat)


def test_options_fp(df_compat, fp):
    # use the set option

    with pd.option_context("io.parquet.engine", "fastparquet"):
        check_round_trip(df_compat)


def test_options_auto(df_compat, fp, pa):
    # use the set option

    with pd.option_context("io.parquet.engine", "auto"):
        check_round_trip(df_compat)


def test_options_get_engine(fp, pa):
    assert isinstance(get_engine("pyarrow"), PyArrowImpl)
    assert isinstance(get_engine("fastparquet"), FastParquetImpl)

    with pd.option_context("io.parquet.engine", "pyarrow"):
        assert isinstance(get_engine("auto"), PyArrowImpl)
        assert isinstance(get_engine("pyarrow"), PyArrowImpl)
        assert isinstance(get_engine("fastparquet"), FastParquetImpl)

    with pd.option_context("io.parquet.engine", "fastparquet"):
        assert isinstance(get_engine("auto"), FastParquetImpl)
        assert isinstance(get_engine("pyarrow"), PyArrowImpl)
        assert isinstance(get_engine("fastparquet"), FastParquetImpl)

    with pd.option_context("io.parquet.engine", "auto"):
        assert isinstance(get_engine("auto"), PyArrowImpl)
        assert isinstance(get_engine("pyarrow"), PyArrowImpl)
        assert isinstance(get_engine("fastparquet"), FastParquetImpl)


def test_get_engine_auto_error_message():
    # Expect different error messages from get_engine(engine="auto")
    # if engines aren't installed vs. are installed but bad version
    from pandas.compat._optional import VERSIONS

    # Do we have engines installed, but a bad version of them?
    pa_min_ver = VERSIONS.get("pyarrow")
    fp_min_ver = VERSIONS.get("fastparquet")
    have_pa_bad_version = (
        False
        if not _HAVE_PYARROW
        else LooseVersion(pyarrow.__version__) < LooseVersion(pa_min_ver)
    )
    have_fp_bad_version = (
        False
        if not _HAVE_FASTPARQUET
        else LooseVersion(fastparquet.__version__) < LooseVersion(fp_min_ver)
    )
    # Do we have usable engines installed?
    have_usable_pa = _HAVE_PYARROW and not have_pa_bad_version
    have_usable_fp = _HAVE_FASTPARQUET and not have_fp_bad_version

    if not have_usable_pa and not have_usable_fp:
        # No usable engines found.
        if have_pa_bad_version:
            match = f"Pandas requires version .{pa_min_ver}. or newer of .pyarrow."
            with pytest.raises(ImportError, match=match):
                get_engine("auto")
        else:
            match = "Missing optional dependency .pyarrow."
            with pytest.raises(ImportError, match=match):
                get_engine("auto")

        if have_fp_bad_version:
            match = f"Pandas requires version .{fp_min_ver}. or newer of .fastparquet."
            with pytest.raises(ImportError, match=match):
                get_engine("auto")
        else:
            match = "Missing optional dependency .fastparquet."
            with pytest.raises(ImportError, match=match):
                get_engine("auto")


def test_cross_engine_pa_fp(df_cross_compat, pa, fp):
    # cross-compat with differing reading/writing engines

    df = df_cross_compat
    with tm.ensure_clean() as path:
        df.to_parquet(path, engine=pa, compression=None)

        result = read_parquet(path, engine=fp)
        tm.assert_frame_equal(result, df)

        result = read_parquet(path, engine=fp, columns=["a", "d"])
        tm.assert_frame_equal(result, df[["a", "d"]])


def test_cross_engine_fp_pa(df_cross_compat, pa, fp):
    # cross-compat with differing reading/writing engines

    if (
        LooseVersion(pyarrow.__version__) < "0.15"
        and LooseVersion(pyarrow.__version__) >= "0.13"
    ):
        pytest.xfail(
            "Reading fastparquet with pyarrow in 0.14 fails: "
            "https://issues.apache.org/jira/browse/ARROW-6492"
        )

    df = df_cross_compat
    with tm.ensure_clean() as path:
        df.to_parquet(path, engine=fp, compression=None)

        with catch_warnings(record=True):
            result = read_parquet(path, engine=pa)
            tm.assert_frame_equal(result, df)

            result = read_parquet(path, engine=pa, columns=["a", "d"])
            tm.assert_frame_equal(result, df[["a", "d"]])


class Base:
    def check_error_on_write(self, df, engine, exc):
        # check that we are raising the exception on writing
        with tm.ensure_clean() as path:
            with pytest.raises(exc):
                to_parquet(df, path, engine, compression=None)


class TestBasic(Base):
    def test_error(self, engine):
        for obj in [
            pd.Series([1, 2, 3]),
            1,
            "foo",
            pd.Timestamp("20130101"),
            np.array([1, 2, 3]),
        ]:
            self.check_error_on_write(obj, engine, ValueError)

    def test_columns_dtypes(self, engine):
        df = pd.DataFrame({"string": list("abc"), "int": list(range(1, 4))})

        # unicode
        df.columns = ["foo", "bar"]
        check_round_trip(df, engine)

    def test_columns_dtypes_invalid(self, engine):
        df = pd.DataFrame({"string": list("abc"), "int": list(range(1, 4))})

        # numeric
        df.columns = [0, 1]
        self.check_error_on_write(df, engine, ValueError)

        # bytes
        df.columns = [b"foo", b"bar"]
        self.check_error_on_write(df, engine, ValueError)

        # python object
        df.columns = [
            datetime.datetime(2011, 1, 1, 0, 0),
            datetime.datetime(2011, 1, 1, 1, 1),
        ]
        self.check_error_on_write(df, engine, ValueError)

    @pytest.mark.parametrize("compression", [None, "gzip", "snappy", "brotli"])
    def test_compression(self, engine, compression):

        if compression == "snappy":
            pytest.importorskip("snappy")

        elif compression == "brotli":
            pytest.importorskip("brotli")

        df = pd.DataFrame({"A": [1, 2, 3]})
        check_round_trip(df, engine, write_kwargs={"compression": compression})

    def test_read_columns(self, engine):
        # GH18154
        df = pd.DataFrame({"string": list("abc"), "int": list(range(1, 4))})

        expected = pd.DataFrame({"string": list("abc")})
        check_round_trip(
            df, engine, expected=expected, read_kwargs={"columns": ["string"]}
        )

    def test_write_index(self, engine):
        check_names = engine != "fastparquet"

        df = pd.DataFrame({"A": [1, 2, 3]})
        check_round_trip(df, engine)

        indexes = [
            [2, 3, 4],
            pd.date_range("20130101", periods=3),
            list("abc"),
            [1, 3, 4],
        ]
        # non-default index
        for index in indexes:
            df.index = index
            if isinstance(index, pd.DatetimeIndex):
                df.index = df.index._with_freq(None)  # freq doesnt round-trip
            check_round_trip(df, engine, check_names=check_names)

        # index with meta-data
        df.index = [0, 1, 2]
        df.index.name = "foo"
        check_round_trip(df, engine)

    def test_write_multiindex(self, pa):
        # Not supported in fastparquet as of 0.1.3 or older pyarrow version
        engine = pa

        df = pd.DataFrame({"A": [1, 2, 3]})
        index = pd.MultiIndex.from_tuples([("a", 1), ("a", 2), ("b", 1)])
        df.index = index
        check_round_trip(df, engine)

    def test_write_column_multiindex(self, engine):
        # column multi-index
        mi_columns = pd.MultiIndex.from_tuples([("a", 1), ("a", 2), ("b", 1)])
        df = pd.DataFrame(np.random.randn(4, 3), columns=mi_columns)
        self.check_error_on_write(df, engine, ValueError)

    def test_multiindex_with_columns(self, pa):
        engine = pa
        dates = pd.date_range("01-Jan-2018", "01-Dec-2018", freq="MS")
        df = pd.DataFrame(np.random.randn(2 * len(dates), 3), columns=list("ABC"))
        index1 = pd.MultiIndex.from_product(
            [["Level1", "Level2"], dates], names=["level", "date"]
        )
        index2 = index1.copy(names=None)
        for index in [index1, index2]:
            df.index = index

            check_round_trip(df, engine)
            check_round_trip(
                df, engine, read_kwargs={"columns": ["A", "B"]}, expected=df[["A", "B"]]
            )

    def test_write_ignoring_index(self, engine):
        # ENH 20768
        # Ensure index=False omits the index from the written Parquet file.
        df = pd.DataFrame({"a": [1, 2, 3], "b": ["q", "r", "s"]})

        write_kwargs = {"compression": None, "index": False}

        # Because we're dropping the index, we expect the loaded dataframe to
        # have the default integer index.
        expected = df.reset_index(drop=True)

        check_round_trip(df, engine, write_kwargs=write_kwargs, expected=expected)

        # Ignore custom index
        df = pd.DataFrame(
            {"a": [1, 2, 3], "b": ["q", "r", "s"]}, index=["zyx", "wvu", "tsr"]
        )

        check_round_trip(df, engine, write_kwargs=write_kwargs, expected=expected)

        # Ignore multi-indexes as well.
        arrays = [
            ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
            ["one", "two", "one", "two", "one", "two", "one", "two"],
        ]
        df = pd.DataFrame(
            {"one": list(range(8)), "two": [-i for i in range(8)]}, index=arrays
        )

        expected = df.reset_index(drop=True)
        check_round_trip(df, engine, write_kwargs=write_kwargs, expected=expected)


class TestParquetPyArrow(Base):
    def test_basic(self, pa, df_full):

        df = df_full

        # additional supported types for pyarrow
        dti = pd.date_range("20130101", periods=3, tz="Europe/Brussels")
        dti = dti._with_freq(None)  # freq doesnt round-trip
        df["datetime_tz"] = dti
        df["bool_with_none"] = [True, None, True]

        check_round_trip(df, pa)

    def test_basic_subset_columns(self, pa, df_full):
        # GH18628

        df = df_full
        # additional supported types for pyarrow
        df["datetime_tz"] = pd.date_range("20130101", periods=3, tz="Europe/Brussels")

        check_round_trip(
            df,
            pa,
            expected=df[["string", "int"]],
            read_kwargs={"columns": ["string", "int"]},
        )

    def test_duplicate_columns(self, pa):
        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3), columns=list("aaa")).copy()
        self.check_error_on_write(df, pa, ValueError)

    def test_unsupported(self, pa):
        if LooseVersion(pyarrow.__version__) < LooseVersion("0.15.1.dev"):
            # period - will be supported using an extension type with pyarrow 1.0
            df = pd.DataFrame({"a": pd.period_range("2013", freq="M", periods=3)})
            # pyarrow 0.11 raises ArrowTypeError
            # older pyarrows raise ArrowInvalid
            self.check_error_on_write(df, pa, Exception)

        # timedelta
        df = pd.DataFrame({"a": pd.timedelta_range("1 day", periods=3)})
        self.check_error_on_write(df, pa, NotImplementedError)

        # mixed python objects
        df = pd.DataFrame({"a": ["a", 1, 2.0]})
        # pyarrow 0.11 raises ArrowTypeError
        # older pyarrows raise ArrowInvalid
        self.check_error_on_write(df, pa, Exception)

    def test_categorical(self, pa):

        # supported in >= 0.7.0
        df = pd.DataFrame()
        df["a"] = pd.Categorical(list("abcdef"))

        # test for null, out-of-order values, and unobserved category
        df["b"] = pd.Categorical(
            ["bar", "foo", "foo", "bar", None, "bar"],
            dtype=pd.CategoricalDtype(["foo", "bar", "baz"]),
        )

        # test for ordered flag
        df["c"] = pd.Categorical(
            ["a", "b", "c", "a", "c", "b"], categories=["b", "c", "d"], ordered=True
        )

        if LooseVersion(pyarrow.__version__) >= LooseVersion("0.15.0"):
            check_round_trip(df, pa)
        else:
            # de-serialized as object for pyarrow < 0.15
            expected = df.astype(object)
            check_round_trip(df, pa, expected=expected)

    def test_s3_roundtrip_explicit_fs(self, df_compat, s3_resource, pa):
        s3fs = pytest.importorskip("s3fs")
        s3 = s3fs.S3FileSystem()
        kw = dict(filesystem=s3)
        check_round_trip(
            df_compat,
            pa,
            path="pandas-test/pyarrow.parquet",
            read_kwargs=kw,
            write_kwargs=kw,
        )

    def test_s3_roundtrip(self, df_compat, s3_resource, pa):
        # GH #19134
        check_round_trip(df_compat, pa, path="s3://pandas-test/pyarrow.parquet")

    @td.skip_if_no("s3fs")
    @pytest.mark.parametrize("partition_col", [["A"], []])
    def test_s3_roundtrip_for_dir(self, df_compat, s3_resource, pa, partition_col):
        # GH #26388
        expected_df = df_compat.copy()

        # GH #35791
        # read_table uses the new Arrow Datasets API since pyarrow 1.0.0
        # Previous behaviour was pyarrow partitioned columns become 'category' dtypes
        # These are added to back of dataframe on read. In new API category dtype is
        # only used if partition field is string, but this changed again to use
        # category dtype for all types (not only strings) in pyarrow 2.0.0
        pa10 = (LooseVersion(pyarrow.__version__) >= LooseVersion("1.0.0")) and (
            LooseVersion(pyarrow.__version__) < LooseVersion("2.0.0")
        )
        if partition_col:
            if pa10:
                partition_col_type = "int32"
            else:
                partition_col_type = "category"

            expected_df[partition_col] = expected_df[partition_col].astype(
                partition_col_type
            )

        check_round_trip(
            df_compat,
            pa,
            expected=expected_df,
            path="s3://pandas-test/parquet_dir",
            write_kwargs={"partition_cols": partition_col, "compression": None},
            check_like=True,
            repeat=1,
        )

    @tm.network
    @td.skip_if_no("pyarrow")
    def test_parquet_read_from_url(self, df_compat):
        url = (
            "https://raw.githubusercontent.com/pandas-dev/pandas/"
            "master/pandas/tests/io/data/parquet/simple.parquet"
        )
        df = pd.read_parquet(url)
        tm.assert_frame_equal(df, df_compat)

    @td.skip_if_no("pyarrow")
    def test_read_file_like_obj_support(self, df_compat):
        buffer = BytesIO()
        df_compat.to_parquet(buffer)
        df_from_buf = pd.read_parquet(buffer)
        tm.assert_frame_equal(df_compat, df_from_buf)

    @td.skip_if_no("pyarrow")
    def test_expand_user(self, df_compat, monkeypatch):
        monkeypatch.setenv("HOME", "TestingUser")
        monkeypatch.setenv("USERPROFILE", "TestingUser")
        with pytest.raises(OSError, match=r".*TestingUser.*"):
            pd.read_parquet("~/file.parquet")
        with pytest.raises(OSError, match=r".*TestingUser.*"):
            df_compat.to_parquet("~/file.parquet")

    def test_partition_cols_supported(self, pa, df_full):
        # GH #23283
        partition_cols = ["bool", "int"]
        df = df_full
        with tm.ensure_clean_dir() as path:
            df.to_parquet(path, partition_cols=partition_cols, compression=None)
            import pyarrow.parquet as pq

            dataset = pq.ParquetDataset(path, validate_schema=False)
            assert len(dataset.partitions.partition_names) == 2
            assert dataset.partitions.partition_names == set(partition_cols)

    def test_partition_cols_string(self, pa, df_full):
        # GH #27117
        partition_cols = "bool"
        partition_cols_list = [partition_cols]
        df = df_full
        with tm.ensure_clean_dir() as path:
            df.to_parquet(path, partition_cols=partition_cols, compression=None)
            import pyarrow.parquet as pq

            dataset = pq.ParquetDataset(path, validate_schema=False)
            assert len(dataset.partitions.partition_names) == 1
            assert dataset.partitions.partition_names == set(partition_cols_list)

    def test_empty_dataframe(self, pa):
        # GH #27339
        df = pd.DataFrame()
        check_round_trip(df, pa)

    def test_write_with_schema(self, pa):
        import pyarrow

        df = pd.DataFrame({"x": [0, 1]})
        schema = pyarrow.schema([pyarrow.field("x", type=pyarrow.bool_())])
        out_df = df.astype(bool)
        check_round_trip(df, pa, write_kwargs={"schema": schema}, expected=out_df)

    @td.skip_if_no("pyarrow", min_version="0.15.0")
    def test_additional_extension_arrays(self, pa):
        # test additional ExtensionArrays that are supported through the
        # __arrow_array__ protocol
        df = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype="Int64"),
                "b": pd.Series([1, 2, 3], dtype="UInt32"),
                "c": pd.Series(["a", None, "c"], dtype="string"),
            }
        )
        if LooseVersion(pyarrow.__version__) >= LooseVersion("0.16.0"):
            expected = df
        else:
            # de-serialized as plain int / object
            expected = df.assign(
                a=df.a.astype("int64"), b=df.b.astype("int64"), c=df.c.astype("object")
            )
        check_round_trip(df, pa, expected=expected)

        df = pd.DataFrame({"a": pd.Series([1, 2, 3, None], dtype="Int64")})
        if LooseVersion(pyarrow.__version__) >= LooseVersion("0.16.0"):
            expected = df
        else:
            # if missing values in integer, currently de-serialized as float
            expected = df.assign(a=df.a.astype("float64"))
        check_round_trip(df, pa, expected=expected)

    @td.skip_if_no("pyarrow", min_version="0.16.0")
    def test_additional_extension_types(self, pa):
        # test additional ExtensionArrays that are supported through the
        # __arrow_array__ protocol + by defining a custom ExtensionType
        df = pd.DataFrame(
            {
                # Arrow does not yet support struct in writing to Parquet (ARROW-1644)
                # "c": pd.arrays.IntervalArray.from_tuples([(0, 1), (1, 2), (3, 4)]),
                "d": pd.period_range("2012-01-01", periods=3, freq="D"),
            }
        )
        check_round_trip(df, pa)

    @td.skip_if_no("pyarrow", min_version="0.14")
    def test_timestamp_nanoseconds(self, pa):
        # with version 2.0, pyarrow defaults to writing the nanoseconds, so
        # this should work without error
        df = pd.DataFrame({"a": pd.date_range("2017-01-01", freq="1n", periods=10)})
        check_round_trip(df, pa, write_kwargs={"version": "2.0"})

    @td.skip_if_no("pyarrow", min_version="0.17")
    def test_filter_row_groups(self, pa):
        # https://github.com/pandas-dev/pandas/issues/26551
        df = pd.DataFrame({"a": list(range(0, 3))})
        with tm.ensure_clean() as path:
            df.to_parquet(path, pa)
            result = read_parquet(
                path, pa, filters=[("a", "==", 0)], use_legacy_dataset=False
            )
        assert len(result) == 1


class TestParquetFastParquet(Base):
    @td.skip_if_no("fastparquet", min_version="0.3.2")
    def test_basic(self, fp, df_full):
        df = df_full

        dti = pd.date_range("20130101", periods=3, tz="US/Eastern")
        dti = dti._with_freq(None)  # freq doesnt round-trip
        df["datetime_tz"] = dti
        df["timedelta"] = pd.timedelta_range("1 day", periods=3)
        check_round_trip(df, fp)

    @pytest.mark.skip(reason="not supported")
    def test_duplicate_columns(self, fp):

        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3), columns=list("aaa")).copy()
        self.check_error_on_write(df, fp, ValueError)

    def test_bool_with_none(self, fp):
        df = pd.DataFrame({"a": [True, None, False]})
        expected = pd.DataFrame({"a": [1.0, np.nan, 0.0]}, dtype="float16")
        check_round_trip(df, fp, expected=expected)

    def test_unsupported(self, fp):

        # period
        df = pd.DataFrame({"a": pd.period_range("2013", freq="M", periods=3)})
        self.check_error_on_write(df, fp, ValueError)

        # mixed
        df = pd.DataFrame({"a": ["a", 1, 2.0]})
        self.check_error_on_write(df, fp, ValueError)

    def test_categorical(self, fp):
        df = pd.DataFrame({"a": pd.Categorical(list("abc"))})
        check_round_trip(df, fp)

    def test_filter_row_groups(self, fp):
        d = {"a": list(range(0, 3))}
        df = pd.DataFrame(d)
        with tm.ensure_clean() as path:
            df.to_parquet(path, fp, compression=None, row_group_offsets=1)
            result = read_parquet(path, fp, filters=[("a", "==", 0)])
        assert len(result) == 1

    def test_s3_roundtrip(self, df_compat, s3_resource, fp):
        # GH #19134
        check_round_trip(df_compat, fp, path="s3://pandas-test/fastparquet.parquet")

    def test_partition_cols_supported(self, fp, df_full):
        # GH #23283
        partition_cols = ["bool", "int"]
        df = df_full
        with tm.ensure_clean_dir() as path:
            df.to_parquet(
                path,
                engine="fastparquet",
                partition_cols=partition_cols,
                compression=None,
            )
            assert os.path.exists(path)
            import fastparquet  # noqa: F811

            actual_partition_cols = fastparquet.ParquetFile(path, False).cats
            assert len(actual_partition_cols) == 2

    def test_partition_cols_string(self, fp, df_full):
        # GH #27117
        partition_cols = "bool"
        df = df_full
        with tm.ensure_clean_dir() as path:
            df.to_parquet(
                path,
                engine="fastparquet",
                partition_cols=partition_cols,
                compression=None,
            )
            assert os.path.exists(path)
            import fastparquet  # noqa: F811

            actual_partition_cols = fastparquet.ParquetFile(path, False).cats
            assert len(actual_partition_cols) == 1

    def test_partition_on_supported(self, fp, df_full):
        # GH #23283
        partition_cols = ["bool", "int"]
        df = df_full
        with tm.ensure_clean_dir() as path:
            df.to_parquet(
                path,
                engine="fastparquet",
                compression=None,
                partition_on=partition_cols,
            )
            assert os.path.exists(path)
            import fastparquet  # noqa: F811

            actual_partition_cols = fastparquet.ParquetFile(path, False).cats
            assert len(actual_partition_cols) == 2

    def test_error_on_using_partition_cols_and_partition_on(self, fp, df_full):
        # GH #23283
        partition_cols = ["bool", "int"]
        df = df_full
        with pytest.raises(ValueError):
            with tm.ensure_clean_dir() as path:
                df.to_parquet(
                    path,
                    engine="fastparquet",
                    compression=None,
                    partition_on=partition_cols,
                    partition_cols=partition_cols,
                )

    def test_empty_dataframe(self, fp):
        # GH #27339
        df = pd.DataFrame()
        expected = df.copy()
        expected.index.name = "index"
        check_round_trip(df, fp, expected=expected)
