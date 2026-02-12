"""test feather-format compat"""

from datetime import datetime
import zoneinfo

import numpy as np
import pytest

from pandas.compat.pyarrow import (
    pa_version_under18p0,
    pa_version_under19p0,
)

import pandas as pd
import pandas._testing as tm

from pandas.io.feather_format import read_feather, to_feather  # isort:skip

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)


pa = pytest.importorskip("pyarrow")


@pytest.mark.single_cpu
class TestFeather:
    def check_error_on_write(self, df, exc, err_msg, temp_file):
        # check that we are raising the exception
        # on writing

        with pytest.raises(exc, match=err_msg):
            to_feather(df, temp_file)

    def check_external_error_on_write(self, df, temp_file):
        # check that we are raising the exception
        # on writing

        with tm.external_error_raised(Exception):
            to_feather(df, temp_file)

    def check_round_trip(
        self, df, temp_file, expected=None, write_kwargs=None, **read_kwargs
    ):
        if write_kwargs is None:
            write_kwargs = {}
        if expected is None:
            expected = df.copy()

        to_feather(df, temp_file, **write_kwargs)

        result = read_feather(temp_file, **read_kwargs)

        tm.assert_frame_equal(result, expected)

    def test_error(self, temp_file):
        msg = "feather only support IO with DataFrames"
        for obj in [
            pd.Series([1, 2, 3]),
            1,
            "foo",
            pd.Timestamp("20130101"),
            np.array([1, 2, 3]),
        ]:
            self.check_error_on_write(obj, ValueError, msg, temp_file)

    def test_basic(self, temp_file):
        tz = zoneinfo.ZoneInfo("US/Eastern")
        df = pd.DataFrame(
            {
                "string": list("abc"),
                "int": list(range(1, 4)),
                "uint": np.arange(3, 6).astype("u1"),
                "float": np.arange(4.0, 7.0, dtype="float64"),
                "float_with_null": [1.0, np.nan, 3],
                "bool": [True, False, True],
                "bool_with_null": [True, np.nan, False],
                "cat": pd.Categorical(list("abc")),
                "dt": pd.DatetimeIndex(
                    list(pd.date_range("20130101", periods=3)), freq=None
                ),
                "dttz": pd.DatetimeIndex(
                    list(pd.date_range("20130101", periods=3, tz=tz)),
                    freq=None,
                ),
                "dt_with_null": [
                    pd.Timestamp("20130101"),
                    pd.NaT,
                    pd.Timestamp("20130103"),
                ],
                "dtns": pd.DatetimeIndex(
                    list(pd.date_range("20130101", periods=3, freq="ns")), freq=None
                ),
            }
        )
        df["periods"] = pd.period_range("2013", freq="M", periods=3)
        df["timedeltas"] = pd.timedelta_range("1 day", periods=3)
        df["intervals"] = pd.interval_range(0, 3, 3)

        assert df.dttz.dtype.tz.key == "US/Eastern"

        expected = df.copy()
        expected.loc[1, "bool_with_null"] = None
        self.check_round_trip(df, temp_file, expected=expected)

    def test_duplicate_columns(self, temp_file):
        # https://github.com/wesm/feather/issues/53
        # not currently able to handle duplicate columns
        df = pd.DataFrame(np.arange(12).reshape(4, 3), columns=list("aaa")).copy()
        self.check_external_error_on_write(df, temp_file)

    def test_read_columns(self, temp_file):
        # GH 24025
        df = pd.DataFrame(
            {
                "col1": list("abc"),
                "col2": list(range(1, 4)),
                "col3": list("xyz"),
                "col4": list(range(4, 7)),
            }
        )
        columns = ["col1", "col3"]
        self.check_round_trip(df, temp_file, expected=df[columns], columns=columns)

    def test_read_columns_different_order(self, temp_file):
        # GH 33878
        df = pd.DataFrame({"A": [1, 2], "B": ["x", "y"], "C": [True, False]})
        expected = df[["B", "A"]]
        self.check_round_trip(df, temp_file, expected, columns=["B", "A"])

    def test_unsupported_other(self, temp_file):
        # mixed python objects
        df = pd.DataFrame({"a": ["a", 1, 2.0]})
        self.check_external_error_on_write(df, temp_file)

    def test_rw_use_threads(self, temp_file):
        df = pd.DataFrame({"A": np.arange(100000)})
        self.check_round_trip(df, temp_file, use_threads=True)
        self.check_round_trip(df, temp_file, use_threads=False)

    def test_path_pathlib(self, temp_file):
        df = pd.DataFrame(
            1.1 * np.arange(120).reshape((30, 4)),
            columns=pd.Index(list("ABCD")),
            index=pd.Index([f"i-{i}" for i in range(30)]),
        ).reset_index()
        result = tm.round_trip_pathlib(df.to_feather, read_feather, temp_file)
        tm.assert_frame_equal(df, result)

    def test_passthrough_keywords(self, temp_file):
        df = pd.DataFrame(
            1.1 * np.arange(120).reshape((30, 4)),
            columns=pd.Index(list("ABCD")),
            index=pd.Index([f"i-{i}" for i in range(30)]),
        ).reset_index()
        self.check_round_trip(df, temp_file, write_kwargs={"version": 1})

    @pytest.mark.network
    @pytest.mark.single_cpu
    def test_http_path(self, feather_file, httpserver):
        # GH 29055
        expected = read_feather(feather_file)
        with open(feather_file, "rb") as f:
            httpserver.serve_content(content=f.read())
            res = read_feather(httpserver.url)
        tm.assert_frame_equal(expected, res)

    def test_read_feather_dtype_backend(
        self, string_storage, dtype_backend, using_infer_string, temp_file
    ):
        # GH#50765
        df = pd.DataFrame(
            {
                "a": pd.Series([1, pd.NA, 3], dtype="Int64"),
                "b": pd.Series([1, 2, 3], dtype="Int64"),
                "c": pd.Series([1.5, pd.NA, 2.5], dtype="Float64"),
                "d": pd.Series([1.5, 2.0, 2.5], dtype="Float64"),
                "e": [True, False, None],
                "f": [True, False, True],
                "g": ["a", "b", "c"],
                "h": ["a", "b", None],
            }
        )

        to_feather(df, temp_file)
        with pd.option_context("mode.string_storage", string_storage):
            result = read_feather(temp_file, dtype_backend=dtype_backend)

        if dtype_backend == "pyarrow":
            pa = pytest.importorskip("pyarrow")
            if using_infer_string:
                string_dtype = pd.ArrowDtype(pa.large_string())
            else:
                string_dtype = pd.ArrowDtype(pa.string())
        else:
            string_dtype = pd.StringDtype(string_storage)

        expected = pd.DataFrame(
            {
                "a": pd.Series([1, pd.NA, 3], dtype="Int64"),
                "b": pd.Series([1, 2, 3], dtype="Int64"),
                "c": pd.Series([1.5, pd.NA, 2.5], dtype="Float64"),
                "d": pd.Series([1.5, 2.0, 2.5], dtype="Float64"),
                "e": pd.Series([True, False, pd.NA], dtype="boolean"),
                "f": pd.Series([True, False, True], dtype="boolean"),
                "g": pd.Series(["a", "b", "c"], dtype=string_dtype),
                "h": pd.Series(["a", "b", None], dtype=string_dtype),
            }
        )

        if dtype_backend == "pyarrow":
            from pandas.arrays import ArrowExtensionArray

            expected = pd.DataFrame(
                {
                    col: ArrowExtensionArray(pa.array(expected[col], from_pandas=True))
                    for col in expected.columns
                }
            )

        if using_infer_string:
            expected.columns = expected.columns.astype(
                pd.StringDtype(string_storage, na_value=np.nan)
            )
        tm.assert_frame_equal(result, expected)

    def test_int_columns_and_index(self, temp_file):
        df = pd.DataFrame({"a": [1, 2, 3]}, index=pd.Index([3, 4, 5], name="test"))
        self.check_round_trip(df, temp_file)

    def test_invalid_dtype_backend(self, temp_file):
        msg = (
            "dtype_backend numpy is invalid, only 'numpy_nullable' and "
            "'pyarrow' are allowed."
        )
        df = pd.DataFrame({"int": list(range(1, 4))})
        df.to_feather(temp_file)
        with pytest.raises(ValueError, match=msg):
            read_feather(temp_file, dtype_backend="numpy")

    def test_string_inference(self, temp_file, using_infer_string):
        # GH#54431
        df = pd.DataFrame(data={"a": ["x", "y"]})
        df.to_feather(temp_file)
        with pd.option_context("future.infer_string", True):
            result = read_feather(temp_file)
        dtype = pd.StringDtype(na_value=np.nan)
        expected = pd.DataFrame(
            data={"a": ["x", "y"]}, dtype=pd.StringDtype(na_value=np.nan)
        )
        expected = pd.DataFrame(
            data={"a": ["x", "y"]},
            dtype=dtype,
            columns=pd.Index(
                ["a"],
                dtype=object
                if pa_version_under19p0 and not using_infer_string
                else dtype,
            ),
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(pa_version_under18p0, reason="not supported before 18.0")
    def test_string_inference_string_view_type(self, temp_file):
        # GH#54798
        import pyarrow as pa
        from pyarrow import feather

        table = pa.table({"a": pa.array([None, "b", "c"], pa.string_view())})
        feather.write_feather(table, temp_file)

        with pd.option_context("future.infer_string", True):
            result = read_feather(temp_file)

            expected = pd.DataFrame(
                data={"a": [None, "b", "c"]}, dtype=pd.StringDtype(na_value=np.nan)
            )
        tm.assert_frame_equal(result, expected)

    def test_out_of_bounds_datetime_to_feather(self, temp_file):
        # GH#47832
        df = pd.DataFrame(
            {
                "date": [
                    datetime.fromisoformat("1654-01-01"),
                    datetime.fromisoformat("1920-01-01"),
                ],
            }
        )
        self.check_round_trip(df, temp_file)
