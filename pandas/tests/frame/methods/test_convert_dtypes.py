import datetime

import numpy as np
import pytest

from pandas.errors import Pandas4Warning
import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


class TestConvertDtypes:
    @pytest.mark.parametrize(
        "convert_integer, expected", [(False, np.dtype("int32")), (True, "Int32")]
    )
    def test_convert_dtypes(self, convert_integer, expected, string_storage):
        # Specific types are tested in tests/series/test_dtypes.py
        # Just check that it works for DataFrame here
        df = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=np.dtype("int32")),
                "b": pd.Series(["x", "y", "z"], dtype=np.dtype("O")),
            }
        )
        with pd.option_context("string_storage", string_storage):
            msg = "keyword in DataFrame.convert_dtypes is deprecated"
            with tm.assert_produces_warning(Pandas4Warning, match=msg):
                result = df.convert_dtypes(True, True, convert_integer, False)
        expected = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=expected),
                "b": pd.Series(["x", "y", "z"], dtype=f"string[{string_storage}]"),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_convert_empty(self):
        # Empty DataFrame can pass convert_dtypes, see GH#40393
        empty_df = pd.DataFrame()
        tm.assert_frame_equal(empty_df, empty_df.convert_dtypes())

    @td.skip_if_no("pyarrow")
    def test_convert_empty_categorical_to_pyarrow(self):
        # GH#59934
        df = pd.DataFrame(
            {
                "A": pd.Categorical([None] * 5),
                "B": pd.Categorical([None] * 5, categories=["B1", "B2"]),
            }
        )
        converted = df.convert_dtypes(dtype_backend="pyarrow")
        expected = df
        tm.assert_frame_equal(converted, expected)

    def test_convert_dtypes_retain_column_names(self):
        # GH#41435
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        df.columns.name = "cols"

        result = df.convert_dtypes()
        tm.assert_index_equal(result.columns, df.columns)
        assert result.columns.name == "cols"

    def test_pyarrow_dtype_backend(self, using_nan_is_na):
        pa = pytest.importorskip("pyarrow")
        df = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=np.dtype("int32")),
                "b": pd.Series(["x", "y", None], dtype=np.dtype("O")),
                "c": pd.Series([True, False, None], dtype=np.dtype("O")),
                "d": pd.Series([np.nan, 100.5, 200], dtype=np.dtype("float")),
                "e": pd.Series(pd.date_range("2022", periods=3, unit="ns")),
                "f": pd.Series(pd.date_range("2022", periods=3, tz="UTC").as_unit("s")),
                "g": pd.Series(pd.timedelta_range("1D", periods=3)),
            }
        )
        result = df.convert_dtypes(dtype_backend="pyarrow")

        item = None if using_nan_is_na else np.nan
        expected = pd.DataFrame(
            {
                "a": pd.arrays.ArrowExtensionArray(
                    pa.array([1, 2, 3], type=pa.int32())
                ),
                "b": pd.arrays.ArrowExtensionArray(pa.array(["x", "y", None])),
                "c": pd.arrays.ArrowExtensionArray(pa.array([True, False, None])),
                "d": pd.arrays.ArrowExtensionArray(pa.array([item, 100.5, 200.0])),
                "e": pd.arrays.ArrowExtensionArray(
                    pa.array(
                        [
                            datetime.datetime(2022, 1, 1),
                            datetime.datetime(2022, 1, 2),
                            datetime.datetime(2022, 1, 3),
                        ],
                        type=pa.timestamp(unit="ns"),
                    )
                ),
                "f": pd.arrays.ArrowExtensionArray(
                    pa.array(
                        [
                            datetime.datetime(2022, 1, 1),
                            datetime.datetime(2022, 1, 2),
                            datetime.datetime(2022, 1, 3),
                        ],
                        type=pa.timestamp(unit="s", tz="UTC"),
                    )
                ),
                "g": pd.arrays.ArrowExtensionArray(
                    pa.array(
                        [
                            datetime.timedelta(1),
                            datetime.timedelta(2),
                            datetime.timedelta(3),
                        ],
                        type=pa.duration("us"),
                    )
                ),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_pyarrow_dtype_backend_already_pyarrow(self):
        pytest.importorskip("pyarrow")
        expected = pd.DataFrame([1, 2, 3], dtype="int64[pyarrow]")
        result = expected.convert_dtypes(dtype_backend="pyarrow")
        tm.assert_frame_equal(result, expected)

    def test_pyarrow_dtype_backend_from_pandas_nullable(self):
        pa = pytest.importorskip("pyarrow")
        df = pd.DataFrame(
            {
                "a": pd.Series([1, 2, None], dtype="Int32"),
                "b": pd.Series(["x", "y", None], dtype="string[python]"),
                "c": pd.Series([True, False, None], dtype="boolean"),
                "d": pd.Series([None, 100.5, 200], dtype="Float64"),
            }
        )
        result = df.convert_dtypes(dtype_backend="pyarrow")
        expected = pd.DataFrame(
            {
                "a": pd.arrays.ArrowExtensionArray(
                    pa.array([1, 2, None], type=pa.int32())
                ),
                "b": pd.arrays.ArrowExtensionArray(pa.array(["x", "y", None])),
                "c": pd.arrays.ArrowExtensionArray(pa.array([True, False, None])),
                "d": pd.arrays.ArrowExtensionArray(pa.array([None, 100.5, 200.0])),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_pyarrow_dtype_empty_object(self):
        # GH 50970
        pytest.importorskip("pyarrow")
        expected = pd.DataFrame(columns=[0])
        result = expected.convert_dtypes(dtype_backend="pyarrow")
        tm.assert_frame_equal(result, expected)

    def test_pyarrow_engine_lines_false(self):
        # GH 48893
        df = pd.DataFrame({"a": [1, 2, 3]})
        msg = (
            "dtype_backend numpy is invalid, only 'numpy_nullable', "
            "'pyarrow' and None are allowed."
        )
        with pytest.raises(ValueError, match=msg):
            df.convert_dtypes(dtype_backend="numpy")

    def test_pyarrow_backend_no_conversion(self):
        # GH#52872
        pytest.importorskip("pyarrow")
        df = pd.DataFrame({"a": [1, 2], "b": 1.5, "c": True, "d": "x"})
        expected = df.copy()
        msg = "keyword in DataFrame.convert_dtypes is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.convert_dtypes(
                convert_floating=False,
                convert_integer=False,
                convert_boolean=False,
                convert_string=False,
                dtype_backend="pyarrow",
            )
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_pyarrow_to_np_nullable(self):
        # GH 53648
        pytest.importorskip("pyarrow")
        ser = pd.DataFrame(range(2), dtype="int32[pyarrow]")
        result = ser.convert_dtypes(dtype_backend="numpy_nullable")
        expected = pd.DataFrame(range(2), dtype="Int32")
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_pyarrow_timestamp(self):
        # GH 54191
        pytest.importorskip("pyarrow")
        ser = pd.Series(pd.date_range("2020-01-01", "2020-01-02", freq="1min"))
        expected = ser.astype("timestamp[ms][pyarrow]")
        result = expected.convert_dtypes(dtype_backend="pyarrow")
        tm.assert_series_equal(result, expected)

    def test_convert_dtypes_avoid_block_splitting(self):
        # GH#55341
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": "a"})
        msg = "The convert_integer keyword in DataFrame.convert_dtypes is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.convert_dtypes(convert_integer=False)
        expected = pd.DataFrame(
            {
                "a": [1, 2, 3],
                "b": [4, 5, 6],
                "c": pd.Series(["a"] * 3, dtype="string"),
            }
        )
        tm.assert_frame_equal(result, expected)
        assert result._mgr.nblocks == 2

    def test_convert_dtypes_from_arrow(self):
        # GH#56581
        df = pd.DataFrame([["a", datetime.time(18, 12)]], columns=["a", "b"])
        result = df.convert_dtypes()
        expected = df.astype({"a": "string"})
        tm.assert_frame_equal(result, expected)

    def test_convert_dtype_pyarrow_timezone_preserve(self):
        # GH 60237
        pytest.importorskip("pyarrow")
        df = pd.DataFrame(
            {
                "timestamps": pd.Series(
                    pd.to_datetime(range(5), utc=True, unit="h"),
                    dtype="timestamp[ns, tz=UTC][pyarrow]",
                )
            }
        )
        result = df.convert_dtypes(dtype_backend="pyarrow")
        expected = df.copy()
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_complex(self):
        # GH 60129
        df = pd.DataFrame({"a": [1.0 + 5.0j, 1.5 - 3.0j], "b": [1, 2]})
        expected = pd.DataFrame(
            {
                "a": pd.array([1.0 + 5.0j, 1.5 - 3.0j], dtype="complex128"),
                "b": pd.array([1, 2], dtype="Int64"),
            }
        )
        result = df.convert_dtypes()
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_mixed_column_after_slice(self):
        # GH#64702
        df = pd.DataFrame(data=[[1, "a"], [2, "b"], ["c", 3]], columns=["col1", "col2"])
        df = df.loc[[0, 1]].copy()
        result = df.convert_dtypes()
        expected = pd.DataFrame(
            {
                "col1": pd.array([1, 2], dtype="Int64"),
                "col2": pd.array(["a", "b"], dtype="string"),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_int_out_of_range(self):
        # GH#66517 only the column that overflows int64/uint64 stays object
        df = pd.DataFrame({"a": [2**64, 1], "b": [1, 2]}, dtype=object)
        result = df.convert_dtypes()
        expected = pd.DataFrame(
            {
                "a": pd.Series([2**64, 1], dtype=object),
                "b": pd.Series([1, 2], dtype="Int64"),
            }
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "kwarg",
        [
            "infer_objects",
            "convert_string",
            "convert_integer",
            "convert_boolean",
            "convert_floating",
        ],
    )
    def test_convert_dtypes_deprecated_kwargs(self, kwarg):
        # GH#62022
        df = pd.DataFrame({"a": [1, 2, 3]})
        msg = f"The {kwarg} keyword in DataFrame.convert_dtypes"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            df.convert_dtypes(**{kwarg: True})

    def test_convert_dtypes_numpy_backend(self):
        # GH#35694
        df = pd.DataFrame(
            {
                "int": pd.array([1, 2, 3], dtype="Int32"),
                "int_na": pd.array([1, None, 3], dtype="UInt8"),
                "bool": pd.array([True, False, True], dtype="boolean"),
                "bool_na": pd.array([True, None, False], dtype="boolean"),
                "float_na": pd.array([1.5, None, 3.0], dtype="Float32"),
                "string": pd.array(["x", None, "z"], dtype="string"),
                "numpy": np.array([1, 2, 3], dtype=np.int16),
                "cat": pd.Categorical(["a", "b", "a"]),
            }
        )
        result = df.convert_dtypes(dtype_backend=None)
        expected = pd.DataFrame(
            {
                "int": np.array([1, 2, 3], dtype=np.int32),
                "int_na": np.array([1.0, np.nan, 3.0]),
                "bool": np.array([True, False, True]),
                "bool_na": np.array([True, np.nan, False], dtype=object),
                "float_na": np.array([1.5, np.nan, 3.0], dtype=np.float32),
                "string": pd.Series(["x", np.nan, "z"], dtype="str"),
                "numpy": np.array([1, 2, 3], dtype=np.int16),
                "cat": pd.Categorical(["a", "b", "a"]),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_numpy_backend_pyarrow(self):
        # GH#35694
        pa = pytest.importorskip("pyarrow")
        ts = pd.Timestamp("2020-01-01")
        date = datetime.date(2020, 1, 1)
        df = pd.DataFrame(
            {
                "int": pd.array([1, 2], dtype="int8[pyarrow]"),
                "int_na": pd.array([1, None], dtype="int64[pyarrow]"),
                "bool_na": pd.array([True, None], dtype="bool[pyarrow]"),
                "string": pd.array(["x", None], dtype=pd.ArrowDtype(pa.string())),
                "string_pyarrow": pd.array(["x", None], dtype="string[pyarrow]"),
                "ts": pd.array([ts, None], dtype="timestamp[us][pyarrow]"),
                "ts_tz": pd.array(
                    [ts.tz_localize("UTC"), None], dtype="timestamp[s, tz=UTC][pyarrow]"
                ),
                "duration": pd.array(
                    [pd.Timedelta(1, "s"), None], dtype="duration[ms][pyarrow]"
                ),
                "date": pd.array([date, None], dtype="date32[pyarrow]"),
            }
        )
        result = df.convert_dtypes(dtype_backend=None)
        expected = pd.DataFrame(
            {
                "int": np.array([1, 2], dtype=np.int8),
                "int_na": np.array([1.0, np.nan]),
                "bool_na": np.array([True, np.nan], dtype=object),
                "string": pd.Series(["x", np.nan], dtype="str"),
                "string_pyarrow": pd.Series(["x", np.nan], dtype="str"),
                "ts": pd.Series([ts, None], dtype="M8[us]"),
                "ts_tz": pd.Series([ts, None], dtype="M8[s]").dt.tz_localize("UTC"),
                "duration": pd.Series([pd.Timedelta(1, "s"), None], dtype="m8[ms]"),
                "date": df["date"],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_convert_dtypes_numpy_backend_roundtrip(self):
        # GH#35694
        df = pd.DataFrame(
            {"a": np.array([1, 2], dtype=np.int32), "b": [1.5, np.nan], "c": ["x", "y"]}
        )
        result = df.convert_dtypes().convert_dtypes(dtype_backend=None)
        tm.assert_frame_equal(result, df)

    def test_convert_dtypes_numpy_backend_deprecated_kwargs_raise(self):
        # GH#35694
        df = pd.DataFrame({"a": [1, 2, 3]})
        msg = "Cannot pass convert_integer with dtype_backend=None"
        with pytest.raises(ValueError, match=msg):
            df.convert_dtypes(convert_integer=False, dtype_backend=None)
