from collections import OrderedDict
from datetime import timedelta
from io import StringIO
import json
import os

import numpy as np
import pytest

from pandas.compat import PY35, is_platform_32bit
import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, DatetimeIndex, Series, Timestamp, read_json
import pandas.util.testing as tm
from pandas.util.testing import (
    assert_frame_equal,
    assert_series_equal,
    ensure_clean,
    network,
)

_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()

_frame = DataFrame(_seriesd)
_frame2 = DataFrame(_seriesd, columns=["D", "C", "B", "A"])
_intframe = DataFrame({k: v.astype(np.int64) for k, v in _seriesd.items()})

_tsframe = DataFrame(_tsd)
_cat_frame = _frame.copy()
cat = ["bah"] * 5 + ["bar"] * 5 + ["baz"] * 5 + ["foo"] * (len(_cat_frame) - 15)
_cat_frame.index = pd.CategoricalIndex(cat, name="E")
_cat_frame["E"] = list(reversed(cat))
_cat_frame["sort"] = np.arange(len(_cat_frame), dtype="int64")

_mixed_frame = _frame.copy()


class TestPandasContainer:
    @pytest.fixture(scope="function", autouse=True)
    def setup(self, datapath):
        self.dirpath = datapath("io", "json", "data")

        self.ts = tm.makeTimeSeries()
        self.ts.name = "ts"

        self.series = tm.makeStringSeries()
        self.series.name = "series"

        self.objSeries = tm.makeObjectSeries()
        self.objSeries.name = "objects"

        self.empty_series = Series([], index=[])
        self.empty_frame = DataFrame()

        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.intframe = _intframe.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()
        self.categorical = _cat_frame.copy()

        yield

        del self.dirpath

        del self.ts

        del self.series

        del self.objSeries

        del self.empty_series
        del self.empty_frame

        del self.frame
        del self.frame2
        del self.intframe
        del self.tsframe
        del self.mixed_frame

    def test_frame_double_encoded_labels(self, orient):
        df = DataFrame(
            [["a", "b"], ["c", "d"]],
            index=['index " 1', "index / 2"],
            columns=["a \\ b", "y / z"],
        )

        result = read_json(df.to_json(orient=orient), orient=orient)
        expected = df.copy()

        if orient == "records" or orient == "values":
            expected = expected.reset_index(drop=True)
        if orient == "values":
            expected.columns = range(len(expected.columns))

        assert_frame_equal(result, expected)

    @pytest.mark.parametrize("orient", ["split", "records", "values"])
    def test_frame_non_unique_index(self, orient):
        df = DataFrame([["a", "b"], ["c", "d"]], index=[1, 1], columns=["x", "y"])
        result = read_json(df.to_json(orient=orient), orient=orient)
        expected = df.copy()

        if orient == "records" or orient == "values":
            expected = expected.reset_index(drop=True)
        if orient == "values":
            expected.columns = range(len(expected.columns))

        assert_frame_equal(result, expected)

    @pytest.mark.parametrize("orient", ["index", "columns"])
    def test_frame_non_unique_index_raises(self, orient):
        df = DataFrame([["a", "b"], ["c", "d"]], index=[1, 1], columns=["x", "y"])
        msg = "DataFrame index must be unique for orient='{}'".format(orient)
        with pytest.raises(ValueError, match=msg):
            df.to_json(orient=orient)

    @pytest.mark.parametrize("orient", ["split", "values"])
    @pytest.mark.parametrize(
        "data",
        [
            [["a", "b"], ["c", "d"]],
            [[1.5, 2.5], [3.5, 4.5]],
            [[1, 2.5], [3, 4.5]],
            [[Timestamp("20130101"), 3.5], [Timestamp("20130102"), 4.5]],
        ],
    )
    def test_frame_non_unique_columns(self, orient, data):
        df = DataFrame(data, index=[1, 2], columns=["x", "x"])

        result = read_json(
            df.to_json(orient=orient), orient=orient, convert_dates=["x"]
        )
        if orient == "values":
            expected = pd.DataFrame(data)
            if expected.iloc[:, 0].dtype == "datetime64[ns]":
                # orient == "values" by default will write Timestamp objects out
                # in milliseconds; these are internally stored in nanosecond,
                # so divide to get where we need
                # TODO: a to_epoch method would also solve; see GH 14772
                expected.iloc[:, 0] = expected.iloc[:, 0].astype(np.int64) // 1000000
        elif orient == "split":
            expected = df

        assert_frame_equal(result, expected)

    @pytest.mark.parametrize("orient", ["index", "columns", "records"])
    def test_frame_non_unique_columns_raises(self, orient):
        df = DataFrame([["a", "b"], ["c", "d"]], index=[1, 2], columns=["x", "x"])

        msg = "DataFrame columns must be unique for orient='{}'".format(orient)
        with pytest.raises(ValueError, match=msg):
            df.to_json(orient=orient)

    def test_frame_from_json_to_json(self):
        def _check_orient(
            df,
            orient,
            dtype=None,
            numpy=False,
            convert_axes=True,
            check_dtype=True,
            raise_ok=None,
            sort=None,
            check_index_type=True,
            check_column_type=True,
            check_numpy_dtype=False,
        ):
            if sort is not None:
                df = df.sort_values(sort)
            else:
                df = df.sort_index()

            # if we are not unique, then check that we are raising ValueError
            # for the appropriate orients
            if not df.index.is_unique and orient in ["index", "columns"]:
                msg = "DataFrame index must be unique for orient='{}'".format(orient)
                with pytest.raises(ValueError, match=msg):
                    df.to_json(orient=orient)
                return
            if not df.columns.is_unique and orient in ["index", "columns", "records"]:
                # TODO: not executed. fix this.
                with pytest.raises(ValueError, match="ksjkajksfjksjfkjs"):
                    df.to_json(orient=orient)
                return

            dfjson = df.to_json(orient=orient)

            try:
                unser = read_json(
                    dfjson,
                    orient=orient,
                    dtype=dtype,
                    numpy=numpy,
                    convert_axes=convert_axes,
                )
            except Exception as detail:
                if raise_ok is not None:
                    if isinstance(detail, raise_ok):
                        return
                raise

            if sort is not None and sort in unser.columns:
                unser = unser.sort_values(sort)
            else:
                unser = unser.sort_index()

            if not dtype:
                check_dtype = False

            if not convert_axes and df.index.dtype.type == np.datetime64:
                unser.index = DatetimeIndex(unser.index.values.astype("i8") * 1e6)
            if orient == "records":
                # index is not captured in this orientation
                tm.assert_almost_equal(
                    df.values, unser.values, check_dtype=check_numpy_dtype
                )
                tm.assert_index_equal(
                    df.columns, unser.columns, exact=check_column_type
                )
            elif orient == "values":
                # index and cols are not captured in this orientation
                if numpy is True and df.shape == (0, 0):
                    assert unser.shape[0] == 0
                else:
                    tm.assert_almost_equal(
                        df.values, unser.values, check_dtype=check_numpy_dtype
                    )
            elif orient == "split":
                # index and col labels might not be strings
                unser.index = [str(i) for i in unser.index]
                unser.columns = [str(i) for i in unser.columns]

                if sort is None:
                    unser = unser.sort_index()
                tm.assert_almost_equal(
                    df.values, unser.values, check_dtype=check_numpy_dtype
                )
            else:
                if convert_axes:
                    tm.assert_frame_equal(
                        df,
                        unser,
                        check_dtype=check_dtype,
                        check_index_type=check_index_type,
                        check_column_type=check_column_type,
                    )
                else:
                    tm.assert_frame_equal(
                        df, unser, check_less_precise=False, check_dtype=check_dtype
                    )

        def _check_all_orients(
            df,
            dtype=None,
            convert_axes=True,
            raise_ok=None,
            sort=None,
            check_index_type=True,
            check_column_type=True,
        ):

            # numpy=False
            if convert_axes:
                _check_orient(
                    df,
                    "columns",
                    dtype=dtype,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "records",
                    dtype=dtype,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "split",
                    dtype=dtype,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "index",
                    dtype=dtype,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "values",
                    dtype=dtype,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )

            _check_orient(df, "columns", dtype=dtype, convert_axes=False, sort=sort)
            _check_orient(df, "records", dtype=dtype, convert_axes=False, sort=sort)
            _check_orient(df, "split", dtype=dtype, convert_axes=False, sort=sort)
            _check_orient(df, "index", dtype=dtype, convert_axes=False, sort=sort)
            _check_orient(df, "values", dtype=dtype, convert_axes=False, sort=sort)

            # numpy=True and raise_ok might be not None, so ignore the error
            if convert_axes:
                _check_orient(
                    df,
                    "columns",
                    dtype=dtype,
                    numpy=True,
                    raise_ok=raise_ok,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "records",
                    dtype=dtype,
                    numpy=True,
                    raise_ok=raise_ok,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "split",
                    dtype=dtype,
                    numpy=True,
                    raise_ok=raise_ok,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "index",
                    dtype=dtype,
                    numpy=True,
                    raise_ok=raise_ok,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )
                _check_orient(
                    df,
                    "values",
                    dtype=dtype,
                    numpy=True,
                    raise_ok=raise_ok,
                    sort=sort,
                    check_index_type=False,
                    check_column_type=False,
                )

            _check_orient(
                df,
                "columns",
                dtype=dtype,
                numpy=True,
                convert_axes=False,
                raise_ok=raise_ok,
                sort=sort,
            )
            _check_orient(
                df,
                "records",
                dtype=dtype,
                numpy=True,
                convert_axes=False,
                raise_ok=raise_ok,
                sort=sort,
            )
            _check_orient(
                df,
                "split",
                dtype=dtype,
                numpy=True,
                convert_axes=False,
                raise_ok=raise_ok,
                sort=sort,
            )
            _check_orient(
                df,
                "index",
                dtype=dtype,
                numpy=True,
                convert_axes=False,
                raise_ok=raise_ok,
                sort=sort,
            )
            _check_orient(
                df,
                "values",
                dtype=dtype,
                numpy=True,
                convert_axes=False,
                raise_ok=raise_ok,
                sort=sort,
            )

        # basic
        _check_all_orients(self.frame)
        assert self.frame.to_json() == self.frame.to_json(orient="columns")

        _check_all_orients(self.intframe, dtype=self.intframe.values.dtype)
        _check_all_orients(self.intframe, dtype=False)

        # big one
        # index and columns are strings as all unserialised JSON object keys
        # are assumed to be strings
        biggie = DataFrame(
            np.zeros((200, 4)),
            columns=[str(i) for i in range(4)],
            index=[str(i) for i in range(200)],
        )
        _check_all_orients(biggie, dtype=False, convert_axes=False)

        # dtypes
        _check_all_orients(
            DataFrame(biggie, dtype=np.float64), dtype=np.float64, convert_axes=False
        )
        _check_all_orients(
            DataFrame(biggie, dtype=np.int), dtype=np.int, convert_axes=False
        )
        _check_all_orients(
            DataFrame(biggie, dtype="U3"),
            dtype="U3",
            convert_axes=False,
            raise_ok=ValueError,
        )

        # categorical
        _check_all_orients(self.categorical, sort="sort", raise_ok=ValueError)

        # empty
        _check_all_orients(
            self.empty_frame, check_index_type=False, check_column_type=False
        )

        # time series data
        _check_all_orients(self.tsframe)

        # mixed data
        index = pd.Index(["a", "b", "c", "d", "e"])
        data = {
            "A": [0.0, 1.0, 2.0, 3.0, 4.0],
            "B": [0.0, 1.0, 0.0, 1.0, 0.0],
            "C": ["foo1", "foo2", "foo3", "foo4", "foo5"],
            "D": [True, False, True, False, True],
        }
        df = DataFrame(data=data, index=index)
        _check_orient(df, "split", check_dtype=False)
        _check_orient(df, "records", check_dtype=False)
        _check_orient(df, "values", check_dtype=False)
        _check_orient(df, "columns", check_dtype=False)
        # index oriented is problematic as it is read back in in a transposed
        # state, so the columns are interpreted as having mixed data and
        # given object dtypes.
        # force everything to have object dtype beforehand
        _check_orient(df.transpose().transpose(), "index", dtype=False)

    @pytest.mark.parametrize(
        "data,msg,orient",
        [
            ('{"key":b:a:d}', "Expected object or value", "columns"),
            # too few indices
            (
                '{"columns":["A","B"],'
                '"index":["2","3"],'
                '"data":[[1.0,"1"],[2.0,"2"],[null,"3"]]}',
                r"Shape of passed values is \(3, 2\), indices imply \(2, 2\)",
                "split",
            ),
            # too many columns
            (
                '{"columns":["A","B","C"],'
                '"index":["1","2","3"],'
                '"data":[[1.0,"1"],[2.0,"2"],[null,"3"]]}',
                "3 columns passed, passed data had 2 columns",
                "split",
            ),
            # bad key
            (
                '{"badkey":["A","B"],'
                '"index":["2","3"],'
                '"data":[[1.0,"1"],[2.0,"2"],[null,"3"]]}',
                r"unexpected key\(s\): badkey",
                "split",
            ),
        ],
    )
    def test_frame_from_json_bad_data_raises(self, data, msg, orient):
        with pytest.raises(ValueError, match=msg):
            read_json(StringIO(data), orient=orient)

    @pytest.mark.parametrize("dtype", [True, False])
    @pytest.mark.parametrize("convert_axes", [True, False])
    @pytest.mark.parametrize("numpy", [True, False])
    def test_frame_from_json_missing_data(self, orient, convert_axes, numpy, dtype):
        num_df = DataFrame([[1, 2], [4, 5, 6]])
        result = read_json(
            num_df.to_json(orient=orient),
            orient=orient,
            convert_axes=convert_axes,
            dtype=dtype,
        )
        assert np.isnan(result.iloc[0, 2])

        obj_df = DataFrame([["1", "2"], ["4", "5", "6"]])
        result = read_json(
            obj_df.to_json(orient=orient),
            orient=orient,
            convert_axes=convert_axes,
            dtype=dtype,
        )
        if not dtype:  # TODO: Special case for object data; maybe a bug?
            assert result.iloc[0, 2] is None
        else:
            assert np.isnan(result.iloc[0, 2])

    @pytest.mark.parametrize("inf", [np.inf, np.NINF])
    @pytest.mark.parametrize("dtype", [True, False])
    def test_frame_infinity(self, orient, inf, dtype):
        # infinities get mapped to nulls which get mapped to NaNs during
        # deserialisation
        df = DataFrame([[1, 2], [4, 5, 6]])
        df.loc[0, 2] = inf
        result = read_json(df.to_json(), dtype=dtype)
        assert np.isnan(result.iloc[0, 2])

    @pytest.mark.skipif(
        is_platform_32bit(), reason="not compliant on 32-bit, xref #15865"
    )
    @pytest.mark.parametrize(
        "value,precision,expected_val",
        [
            (0.95, 1, 1.0),
            (1.95, 1, 2.0),
            (-1.95, 1, -2.0),
            (0.995, 2, 1.0),
            (0.9995, 3, 1.0),
            (0.99999999999999944, 15, 1.0),
        ],
    )
    def test_frame_to_json_float_precision(self, value, precision, expected_val):
        df = pd.DataFrame([dict(a_float=value)])
        encoded = df.to_json(double_precision=precision)
        assert encoded == '{{"a_float":{{"0":{}}}}}'.format(expected_val)

    def test_frame_to_json_except(self):
        df = DataFrame([1, 2, 3])
        msg = "Invalid value 'garbage' for option 'orient'"
        with pytest.raises(ValueError, match=msg):
            df.to_json(orient="garbage")

    def test_frame_empty(self):
        df = DataFrame(columns=["jim", "joe"])
        assert not df._is_mixed_type
        assert_frame_equal(
            read_json(df.to_json(), dtype=dict(df.dtypes)), df, check_index_type=False
        )
        # GH 7445
        result = pd.DataFrame({"test": []}, index=[]).to_json(orient="columns")
        expected = '{"test":{}}'
        assert result == expected

    def test_frame_empty_mixedtype(self):
        # mixed type
        df = DataFrame(columns=["jim", "joe"])
        df["joe"] = df["joe"].astype("i8")
        assert df._is_mixed_type
        assert_frame_equal(
            read_json(df.to_json(), dtype=dict(df.dtypes)), df, check_index_type=False
        )

    def test_frame_mixedtype_orient(self):  # GH10289
        vals = [
            [10, 1, "foo", 0.1, 0.01],
            [20, 2, "bar", 0.2, 0.02],
            [30, 3, "baz", 0.3, 0.03],
            [40, 4, "qux", 0.4, 0.04],
        ]

        df = DataFrame(
            vals, index=list("abcd"), columns=["1st", "2nd", "3rd", "4th", "5th"]
        )

        assert df._is_mixed_type
        right = df.copy()

        for orient in ["split", "index", "columns"]:
            inp = df.to_json(orient=orient)
            left = read_json(inp, orient=orient, convert_axes=False)
            assert_frame_equal(left, right)

        right.index = np.arange(len(df))
        inp = df.to_json(orient="records")
        left = read_json(inp, orient="records", convert_axes=False)
        assert_frame_equal(left, right)

        right.columns = np.arange(df.shape[1])
        inp = df.to_json(orient="values")
        left = read_json(inp, orient="values", convert_axes=False)
        assert_frame_equal(left, right)

    def test_v12_compat(self):
        df = DataFrame(
            [
                [1.56808523, 0.65727391, 1.81021139, -0.17251653],
                [-0.2550111, -0.08072427, -0.03202878, -0.17581665],
                [1.51493992, 0.11805825, 1.629455, -1.31506612],
                [-0.02765498, 0.44679743, 0.33192641, -0.27885413],
                [0.05951614, -2.69652057, 1.28163262, 0.34703478],
            ],
            columns=["A", "B", "C", "D"],
            index=pd.date_range("2000-01-03", "2000-01-07"),
        )
        df["date"] = pd.Timestamp("19920106 18:21:32.12")
        df.iloc[3, df.columns.get_loc("date")] = pd.Timestamp("20130101")
        df["modified"] = df["date"]
        df.iloc[1, df.columns.get_loc("modified")] = pd.NaT

        v12_json = os.path.join(self.dirpath, "tsframe_v012.json")
        df_unser = pd.read_json(v12_json)
        assert_frame_equal(df, df_unser)

        df_iso = df.drop(["modified"], axis=1)
        v12_iso_json = os.path.join(self.dirpath, "tsframe_iso_v012.json")
        df_unser_iso = pd.read_json(v12_iso_json)
        assert_frame_equal(df_iso, df_unser_iso)

    def test_blocks_compat_GH9037(self):
        index = pd.date_range("20000101", periods=10, freq="H")
        df_mixed = DataFrame(
            OrderedDict(
                float_1=[
                    -0.92077639,
                    0.77434435,
                    1.25234727,
                    0.61485564,
                    -0.60316077,
                    0.24653374,
                    0.28668979,
                    -2.51969012,
                    0.95748401,
                    -1.02970536,
                ],
                int_1=[
                    19680418,
                    75337055,
                    99973684,
                    65103179,
                    79373900,
                    40314334,
                    21290235,
                    4991321,
                    41903419,
                    16008365,
                ],
                str_1=[
                    "78c608f1",
                    "64a99743",
                    "13d2ff52",
                    "ca7f4af2",
                    "97236474",
                    "bde7e214",
                    "1a6bde47",
                    "b1190be5",
                    "7a669144",
                    "8d64d068",
                ],
                float_2=[
                    -0.0428278,
                    -1.80872357,
                    3.36042349,
                    -0.7573685,
                    -0.48217572,
                    0.86229683,
                    1.08935819,
                    0.93898739,
                    -0.03030452,
                    1.43366348,
                ],
                str_2=[
                    "14f04af9",
                    "d085da90",
                    "4bcfac83",
                    "81504caf",
                    "2ffef4a9",
                    "08e2f5c4",
                    "07e1af03",
                    "addbd4a7",
                    "1f6a09ba",
                    "4bfc4d87",
                ],
                int_2=[
                    86967717,
                    98098830,
                    51927505,
                    20372254,
                    12601730,
                    20884027,
                    34193846,
                    10561746,
                    24867120,
                    76131025,
                ],
            ),
            index=index,
        )

        # JSON deserialisation always creates unicode strings
        df_mixed.columns = df_mixed.columns.astype("unicode")

        df_roundtrip = pd.read_json(df_mixed.to_json(orient="split"), orient="split")
        assert_frame_equal(
            df_mixed,
            df_roundtrip,
            check_index_type=True,
            check_column_type=True,
            check_frame_type=True,
            by_blocks=True,
            check_exact=True,
        )

    def test_frame_nonprintable_bytes(self):
        # GH14256: failing column caused segfaults, if it is not the last one

        class BinaryThing:
            def __init__(self, hexed):
                self.hexed = hexed
                self.binary = bytes.fromhex(hexed)

            def __str__(self):
                return self.hexed

        hexed = "574b4454ba8c5eb4f98a8f45"
        binthing = BinaryThing(hexed)

        # verify the proper conversion of printable content
        df_printable = DataFrame({"A": [binthing.hexed]})
        assert df_printable.to_json() == '{{"A":{{"0":"{hex}"}}}}'.format(hex=hexed)

        # check if non-printable content throws appropriate Exception
        df_nonprintable = DataFrame({"A": [binthing]})
        msg = "Unsupported UTF-8 sequence length when encoding string"
        with pytest.raises(OverflowError, match=msg):
            df_nonprintable.to_json()

        # the same with multiple columns threw segfaults
        df_mixed = DataFrame({"A": [binthing], "B": [1]}, columns=["A", "B"])
        with pytest.raises(OverflowError):
            df_mixed.to_json()

        # default_handler should resolve exceptions for non-string types
        assert df_nonprintable.to_json(
            default_handler=str
        ) == '{{"A":{{"0":"{hex}"}}}}'.format(hex=hexed)
        assert df_mixed.to_json(
            default_handler=str
        ) == '{{"A":{{"0":"{hex}"}},"B":{{"0":1}}}}'.format(hex=hexed)

    def test_label_overflow(self):
        # GH14256: buffer length not checked when writing label
        df = pd.DataFrame({"bar" * 100000: [1], "foo": [1337]})
        assert df.to_json() == '{{"{bar}":{{"0":1}},"foo":{{"0":1337}}}}'.format(
            bar=("bar" * 100000)
        )

    def test_series_non_unique_index(self):
        s = Series(["a", "b"], index=[1, 1])

        msg = "Series index must be unique for orient='index'"
        with pytest.raises(ValueError, match=msg):
            s.to_json(orient="index")

        assert_series_equal(
            s, read_json(s.to_json(orient="split"), orient="split", typ="series")
        )
        unser = read_json(s.to_json(orient="records"), orient="records", typ="series")
        tm.assert_numpy_array_equal(s.values, unser.values)

    def test_series_default_orient(self):
        assert self.series.to_json() == self.series.to_json(orient="index")

    @pytest.mark.parametrize("numpy", [True, False])
    def test_series_roundtrip_simple(self, orient, numpy):
        data = self.series.to_json(orient=orient)
        result = pd.read_json(data, typ="series", orient=orient, numpy=numpy)
        expected = self.series.copy()

        if not numpy and PY35 and orient in ("index", "columns"):
            expected = expected.sort_index()
        if orient in ("values", "records"):
            expected = expected.reset_index(drop=True)
        if orient != "split":
            expected.name = None

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [False, None])
    @pytest.mark.parametrize("numpy", [True, False])
    def test_series_roundtrip_object(self, orient, numpy, dtype):
        data = self.objSeries.to_json(orient=orient)
        result = pd.read_json(
            data, typ="series", orient=orient, numpy=numpy, dtype=dtype
        )
        expected = self.objSeries.copy()

        if not numpy and PY35 and orient in ("index", "columns"):
            expected = expected.sort_index()
        if orient in ("values", "records"):
            expected = expected.reset_index(drop=True)
        if orient != "split":
            expected.name = None

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("numpy", [True, False])
    def test_series_roundtrip_empty(self, orient, numpy):
        data = self.empty_series.to_json(orient=orient)
        result = pd.read_json(data, typ="series", orient=orient, numpy=numpy)
        expected = self.empty_series.copy()

        # TODO: see what causes inconsistency
        if not numpy and PY35 and orient == "index":
            expected = expected.sort_index()
        if orient in ("values", "records"):
            expected = expected.reset_index(drop=True)
        else:
            expected.index = expected.index.astype(float)

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("numpy", [True, False])
    def test_series_roundtrip_timeseries(self, orient, numpy):
        data = self.ts.to_json(orient=orient)
        result = pd.read_json(data, typ="series", orient=orient, numpy=numpy)
        expected = self.ts.copy()

        if orient in ("values", "records"):
            expected = expected.reset_index(drop=True)
        if orient != "split":
            expected.name = None

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [np.float64, np.int])
    @pytest.mark.parametrize("numpy", [True, False])
    def test_series_roundtrip_numeric(self, orient, numpy, dtype):
        s = Series(range(6), index=["a", "b", "c", "d", "e", "f"])
        data = s.to_json(orient=orient)
        result = pd.read_json(data, typ="series", orient=orient, numpy=numpy)

        expected = s.copy()
        if orient in ("values", "records"):
            expected = expected.reset_index(drop=True)

        tm.assert_series_equal(result, expected)

    def test_series_to_json_except(self):
        s = Series([1, 2, 3])
        msg = "Invalid value 'garbage' for option 'orient'"
        with pytest.raises(ValueError, match=msg):
            s.to_json(orient="garbage")

    def test_series_from_json_precise_float(self):
        s = Series([4.56, 4.56, 4.56])
        result = read_json(s.to_json(), typ="series", precise_float=True)
        assert_series_equal(result, s, check_index_type=False)

    def test_series_with_dtype(self):
        # GH 21986
        s = Series([4.56, 4.56, 4.56])
        result = read_json(s.to_json(), typ="series", dtype=np.int64)
        expected = Series([4] * 3)
        assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype,expected",
        [
            (True, Series(["2000-01-01"], dtype="datetime64[ns]")),
            (False, Series([946684800000])),
        ],
    )
    def test_series_with_dtype_datetime(self, dtype, expected):
        s = Series(["2000-01-01"], dtype="datetime64[ns]")
        data = s.to_json()
        result = pd.read_json(data, typ="series", dtype=dtype)
        assert_series_equal(result, expected)

    def test_frame_from_json_precise_float(self):
        df = DataFrame([[4.56, 4.56, 4.56], [4.56, 4.56, 4.56]])
        result = read_json(df.to_json(), precise_float=True)
        assert_frame_equal(result, df, check_index_type=False, check_column_type=False)

    def test_typ(self):

        s = Series(range(6), index=["a", "b", "c", "d", "e", "f"], dtype="int64")
        result = read_json(s.to_json(), typ=None)
        assert_series_equal(result, s)

    def test_reconstruction_index(self):

        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        result = read_json(df.to_json())

        assert_frame_equal(result, df)

        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=["A", "B", "C"])
        result = read_json(df.to_json())
        assert_frame_equal(result, df)

    def test_path(self):
        with ensure_clean("test.json") as path:
            for df in [
                self.frame,
                self.frame2,
                self.intframe,
                self.tsframe,
                self.mixed_frame,
            ]:
                df.to_json(path)
                read_json(path)

    def test_axis_dates(self):

        # frame
        json = self.tsframe.to_json()
        result = read_json(json)
        assert_frame_equal(result, self.tsframe)

        # series
        json = self.ts.to_json()
        result = read_json(json, typ="series")
        assert_series_equal(result, self.ts, check_names=False)
        assert result.name is None

    def test_convert_dates(self):

        # frame
        df = self.tsframe.copy()
        df["date"] = Timestamp("20130101")

        json = df.to_json()
        result = read_json(json)
        assert_frame_equal(result, df)

        df["foo"] = 1.0
        json = df.to_json(date_unit="ns")

        result = read_json(json, convert_dates=False)
        expected = df.copy()
        expected["date"] = expected["date"].values.view("i8")
        expected["foo"] = expected["foo"].astype("int64")
        assert_frame_equal(result, expected)

        # series
        ts = Series(Timestamp("20130101"), index=self.ts.index)
        json = ts.to_json()
        result = read_json(json, typ="series")
        assert_series_equal(result, ts)

    @pytest.mark.parametrize(
        "infer_word",
        [
            "trade_time",
            "date",
            "datetime",
            "sold_at",
            "modified",
            "timestamp",
            "timestamps",
        ],
    )
    def test_convert_dates_infer(self, infer_word):
        # GH10747
        from pandas.io.json import dumps

        data = [{"id": 1, infer_word: 1036713600000}, {"id": 2}]
        expected = DataFrame(
            [[1, Timestamp("2002-11-08")], [2, pd.NaT]], columns=["id", infer_word]
        )
        result = read_json(dumps(data))[["id", infer_word]]
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "date,date_unit",
        [
            ("20130101 20:43:42.123", None),
            ("20130101 20:43:42", "s"),
            ("20130101 20:43:42.123", "ms"),
            ("20130101 20:43:42.123456", "us"),
            ("20130101 20:43:42.123456789", "ns"),
        ],
    )
    def test_date_format_frame(self, date, date_unit):
        df = self.tsframe.copy()

        df["date"] = Timestamp(date)
        df.iloc[1, df.columns.get_loc("date")] = pd.NaT
        df.iloc[5, df.columns.get_loc("date")] = pd.NaT
        if date_unit:
            json = df.to_json(date_format="iso", date_unit=date_unit)
        else:
            json = df.to_json(date_format="iso")
        result = read_json(json)
        expected = df.copy()
        # expected.index = expected.index.tz_localize("UTC")
        expected["date"] = expected["date"].dt.tz_localize("UTC")
        assert_frame_equal(result, expected)

    def test_date_format_frame_raises(self):
        df = self.tsframe.copy()
        msg = "Invalid value 'foo' for option 'date_unit'"
        with pytest.raises(ValueError, match=msg):
            df.to_json(date_format="iso", date_unit="foo")

    @pytest.mark.parametrize(
        "date,date_unit",
        [
            ("20130101 20:43:42.123", None),
            ("20130101 20:43:42", "s"),
            ("20130101 20:43:42.123", "ms"),
            ("20130101 20:43:42.123456", "us"),
            ("20130101 20:43:42.123456789", "ns"),
        ],
    )
    def test_date_format_series(self, date, date_unit):
        ts = Series(Timestamp(date), index=self.ts.index)
        ts.iloc[1] = pd.NaT
        ts.iloc[5] = pd.NaT
        if date_unit:
            json = ts.to_json(date_format="iso", date_unit=date_unit)
        else:
            json = ts.to_json(date_format="iso")
        result = read_json(json, typ="series")
        expected = ts.copy()
        # expected.index = expected.index.tz_localize("UTC")
        expected = expected.dt.tz_localize("UTC")
        assert_series_equal(result, expected)

    def test_date_format_series_raises(self):
        ts = Series(Timestamp("20130101 20:43:42.123"), index=self.ts.index)
        msg = "Invalid value 'foo' for option 'date_unit'"
        with pytest.raises(ValueError, match=msg):
            ts.to_json(date_format="iso", date_unit="foo")

    @pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
    def test_date_unit(self, unit):
        df = self.tsframe.copy()
        df["date"] = Timestamp("20130101 20:43:42")
        dl = df.columns.get_loc("date")
        df.iloc[1, dl] = Timestamp("19710101 20:43:42")
        df.iloc[2, dl] = Timestamp("21460101 20:43:42")
        df.iloc[4, dl] = pd.NaT

        json = df.to_json(date_format="epoch", date_unit=unit)

        # force date unit
        result = read_json(json, date_unit=unit)
        assert_frame_equal(result, df)

        # detect date unit
        result = read_json(json, date_unit=None)
        assert_frame_equal(result, df)

    def test_weird_nested_json(self):
        # this used to core dump the parser
        s = r"""{
        "status": "success",
        "data": {
        "posts": [
            {
            "id": 1,
            "title": "A blog post",
            "body": "Some useful content"
            },
            {
            "id": 2,
            "title": "Another blog post",
            "body": "More content"
            }
           ]
          }
        }"""

        read_json(s)

    def test_doc_example(self):
        dfj2 = DataFrame(np.random.randn(5, 2), columns=list("AB"))
        dfj2["date"] = Timestamp("20130101")
        dfj2["ints"] = range(5)
        dfj2["bools"] = True
        dfj2.index = pd.date_range("20130101", periods=5)

        json = dfj2.to_json()
        result = read_json(json, dtype={"ints": np.int64, "bools": np.bool_})
        assert_frame_equal(result, result)

    def test_misc_example(self):

        # parsing unordered input fails
        result = read_json('[{"a": 1, "b": 2}, {"b":2, "a" :1}]', numpy=True)
        expected = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])

        error_msg = """DataFrame\\.index are different

DataFrame\\.index values are different \\(100\\.0 %\\)
\\[left\\]:  Index\\(\\['a', 'b'\\], dtype='object'\\)
\\[right\\]: RangeIndex\\(start=0, stop=2, step=1\\)"""
        with pytest.raises(AssertionError, match=error_msg):
            assert_frame_equal(result, expected, check_index_type=False)

        result = read_json('[{"a": 1, "b": 2}, {"b":2, "a" :1}]')
        expected = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
        assert_frame_equal(result, expected)

    @network
    @pytest.mark.single
    def test_round_trip_exception_(self):
        # GH 3867
        csv = "https://raw.github.com/hayd/lahman2012/master/csvs/Teams.csv"
        df = pd.read_csv(csv)
        s = df.to_json()
        result = pd.read_json(s)
        assert_frame_equal(result.reindex(index=df.index, columns=df.columns), df)

    @network
    @pytest.mark.single
    @pytest.mark.parametrize(
        "field,dtype",
        [
            ["created_at", pd.DatetimeTZDtype(tz="UTC")],
            ["closed_at", "datetime64[ns]"],
            ["updated_at", pd.DatetimeTZDtype(tz="UTC")],
        ],
    )
    def test_url(self, field, dtype):
        url = "https://api.github.com/repos/pandas-dev/pandas/issues?per_page=5"  # noqa
        result = read_json(url, convert_dates=True)
        assert result[field].dtype == dtype

    def test_timedelta(self):
        converter = lambda x: pd.to_timedelta(x, unit="ms")

        s = Series([timedelta(23), timedelta(seconds=5)])
        assert s.dtype == "timedelta64[ns]"

        result = pd.read_json(s.to_json(), typ="series").apply(converter)
        assert_series_equal(result, s)

        s = Series([timedelta(23), timedelta(seconds=5)], index=pd.Index([0, 1]))
        assert s.dtype == "timedelta64[ns]"
        result = pd.read_json(s.to_json(), typ="series").apply(converter)
        assert_series_equal(result, s)

        frame = DataFrame([timedelta(23), timedelta(seconds=5)])
        assert frame[0].dtype == "timedelta64[ns]"
        assert_frame_equal(frame, pd.read_json(frame.to_json()).apply(converter))

        frame = DataFrame(
            {
                "a": [timedelta(days=23), timedelta(seconds=5)],
                "b": [1, 2],
                "c": pd.date_range(start="20130101", periods=2),
            }
        )

        result = pd.read_json(frame.to_json(date_unit="ns"))
        result["a"] = pd.to_timedelta(result.a, unit="ns")
        result["c"] = pd.to_datetime(result.c)
        assert_frame_equal(frame, result)

    def test_mixed_timedelta_datetime(self):
        frame = DataFrame(
            {"a": [timedelta(23), pd.Timestamp("20130101")]}, dtype=object
        )

        expected = DataFrame(
            {"a": [pd.Timedelta(frame.a[0]).value, pd.Timestamp(frame.a[1]).value]}
        )
        result = pd.read_json(frame.to_json(date_unit="ns"), dtype={"a": "int64"})
        assert_frame_equal(result, expected, check_index_type=False)

    def test_default_handler(self):
        value = object()
        frame = DataFrame({"a": [7, value]})
        expected = DataFrame({"a": [7, str(value)]})
        result = pd.read_json(frame.to_json(default_handler=str))
        assert_frame_equal(expected, result, check_index_type=False)

    def test_default_handler_indirect(self):
        from pandas.io.json import dumps

        def default(obj):
            if isinstance(obj, complex):
                return [("mathjs", "Complex"), ("re", obj.real), ("im", obj.imag)]
            return str(obj)

        df_list = [
            9,
            DataFrame(
                {"a": [1, "STR", complex(4, -5)], "b": [float("nan"), None, "N/A"]},
                columns=["a", "b"],
            ),
        ]
        expected = (
            '[9,[[1,null],["STR",null],[[["mathjs","Complex"],'
            '["re",4.0],["im",-5.0]],"N\\/A"]]]'
        )
        assert dumps(df_list, default_handler=default, orient="values") == expected

    def test_default_handler_numpy_unsupported_dtype(self):
        # GH12554 to_json raises 'Unhandled numpy dtype 15'
        df = DataFrame(
            {"a": [1, 2.3, complex(4, -5)], "b": [float("nan"), None, complex(1.2, 0)]},
            columns=["a", "b"],
        )
        expected = (
            '[["(1+0j)","(nan+0j)"],'
            '["(2.3+0j)","(nan+0j)"],'
            '["(4-5j)","(1.2+0j)"]]'
        )
        assert df.to_json(default_handler=str, orient="values") == expected

    def test_default_handler_raises(self):
        msg = "raisin"

        def my_handler_raises(obj):
            raise TypeError(msg)

        with pytest.raises(TypeError, match=msg):
            DataFrame({"a": [1, 2, object()]}).to_json(
                default_handler=my_handler_raises
            )
        with pytest.raises(TypeError, match=msg):
            DataFrame({"a": [1, 2, complex(4, -5)]}).to_json(
                default_handler=my_handler_raises
            )

    def test_categorical(self):
        # GH4377 df.to_json segfaults with non-ndarray blocks
        df = DataFrame({"A": ["a", "b", "c", "a", "b", "b", "a"]})
        df["B"] = df["A"]
        expected = df.to_json()

        df["B"] = df["A"].astype("category")
        assert expected == df.to_json()

        s = df["A"]
        sc = df["B"]
        assert s.to_json() == sc.to_json()

    def test_datetime_tz(self):
        # GH4377 df.to_json segfaults with non-ndarray blocks
        tz_range = pd.date_range("20130101", periods=3, tz="US/Eastern")
        tz_naive = tz_range.tz_convert("utc").tz_localize(None)

        df = DataFrame({"A": tz_range, "B": pd.date_range("20130101", periods=3)})

        df_naive = df.copy()
        df_naive["A"] = tz_naive
        expected = df_naive.to_json()
        assert expected == df.to_json()

        stz = Series(tz_range)
        s_naive = Series(tz_naive)
        assert stz.to_json() == s_naive.to_json()

    def test_sparse(self):
        # GH4377 df.to_json segfaults with non-ndarray blocks
        df = pd.DataFrame(np.random.randn(10, 4))
        df.loc[:8] = np.nan

        sdf = df.astype("Sparse")
        expected = df.to_json()
        assert expected == sdf.to_json()

        s = pd.Series(np.random.randn(10))
        s.loc[:8] = np.nan
        ss = s.astype("Sparse")

        expected = s.to_json()
        assert expected == ss.to_json()

    @pytest.mark.parametrize(
        "ts",
        [
            Timestamp("2013-01-10 05:00:00Z"),
            Timestamp("2013-01-10 00:00:00", tz="US/Eastern"),
            Timestamp("2013-01-10 00:00:00-0500"),
        ],
    )
    def test_tz_is_utc(self, ts):
        from pandas.io.json import dumps

        exp = '"2013-01-10T05:00:00.000Z"'

        assert dumps(ts, iso_dates=True) == exp
        dt = ts.to_pydatetime()
        assert dumps(dt, iso_dates=True) == exp

    @pytest.mark.parametrize(
        "tz_range",
        [
            pd.date_range("2013-01-01 05:00:00Z", periods=2),
            pd.date_range("2013-01-01 00:00:00", periods=2, tz="US/Eastern"),
            pd.date_range("2013-01-01 00:00:00-0500", periods=2),
        ],
    )
    def test_tz_range_is_utc(self, tz_range):
        from pandas.io.json import dumps

        exp = '["2013-01-01T05:00:00.000Z","2013-01-02T05:00:00.000Z"]'
        dfexp = (
            '{"DT":{'
            '"0":"2013-01-01T05:00:00.000Z",'
            '"1":"2013-01-02T05:00:00.000Z"}}'
        )

        assert dumps(tz_range, iso_dates=True) == exp
        dti = pd.DatetimeIndex(tz_range)
        assert dumps(dti, iso_dates=True) == exp
        df = DataFrame({"DT": dti})
        result = dumps(df, iso_dates=True)
        assert result == dfexp

    def test_read_inline_jsonl(self):
        # GH9180
        result = read_json('{"a": 1, "b": 2}\n{"b":2, "a" :1}\n', lines=True)
        expected = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
        assert_frame_equal(result, expected)

    @td.skip_if_not_us_locale
    def test_read_s3_jsonl(self, s3_resource):
        # GH17200

        result = read_json("s3n://pandas-test/items.jsonl", lines=True)
        expected = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
        assert_frame_equal(result, expected)

    def test_read_local_jsonl(self):
        # GH17200
        with ensure_clean("tmp_items.json") as path:
            with open(path, "w") as infile:
                infile.write('{"a": 1, "b": 2}\n{"b":2, "a" :1}\n')
            result = read_json(path, lines=True)
            expected = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
            assert_frame_equal(result, expected)

    def test_read_jsonl_unicode_chars(self):
        # GH15132: non-ascii unicode characters
        # \u201d == RIGHT DOUBLE QUOTATION MARK

        # simulate file handle
        json = '{"a": "foo”", "b": "bar"}\n{"a": "foo", "b": "bar"}\n'
        json = StringIO(json)
        result = read_json(json, lines=True)
        expected = DataFrame([["foo\u201d", "bar"], ["foo", "bar"]], columns=["a", "b"])
        assert_frame_equal(result, expected)

        # simulate string
        json = '{"a": "foo”", "b": "bar"}\n{"a": "foo", "b": "bar"}\n'
        result = read_json(json, lines=True)
        expected = DataFrame([["foo\u201d", "bar"], ["foo", "bar"]], columns=["a", "b"])
        assert_frame_equal(result, expected)

    def test_read_json_large_numbers(self):
        # GH18842
        json = '{"articleId": "1404366058080022500245"}'
        json = StringIO(json)
        result = read_json(json, typ="series")
        expected = Series(1.404366e21, index=["articleId"])
        assert_series_equal(result, expected)

        json = '{"0": {"articleId": "1404366058080022500245"}}'
        json = StringIO(json)
        result = read_json(json)
        expected = DataFrame(1.404366e21, index=["articleId"], columns=[0])
        assert_frame_equal(result, expected)

    def test_to_jsonl(self):
        # GH9180
        df = DataFrame([[1, 2], [1, 2]], columns=["a", "b"])
        result = df.to_json(orient="records", lines=True)
        expected = '{"a":1,"b":2}\n{"a":1,"b":2}'
        assert result == expected

        df = DataFrame([["foo}", "bar"], ['foo"', "bar"]], columns=["a", "b"])
        result = df.to_json(orient="records", lines=True)
        expected = '{"a":"foo}","b":"bar"}\n{"a":"foo\\"","b":"bar"}'
        assert result == expected
        assert_frame_equal(pd.read_json(result, lines=True), df)

        # GH15096: escaped characters in columns and data
        df = DataFrame([["foo\\", "bar"], ['foo"', "bar"]], columns=["a\\", "b"])
        result = df.to_json(orient="records", lines=True)
        expected = '{"a\\\\":"foo\\\\","b":"bar"}\n' '{"a\\\\":"foo\\"","b":"bar"}'
        assert result == expected
        assert_frame_equal(pd.read_json(result, lines=True), df)

    # TODO: there is a near-identical test for pytables; can we share?
    def test_latin_encoding(self):
        # GH 13774
        pytest.skip("encoding not implemented in .to_json(), xref #13774")

        values = [
            [b"E\xc9, 17", b"", b"a", b"b", b"c"],
            [b"E\xc9, 17", b"a", b"b", b"c"],
            [b"EE, 17", b"", b"a", b"b", b"c"],
            [b"E\xc9, 17", b"\xf8\xfc", b"a", b"b", b"c"],
            [b"", b"a", b"b", b"c"],
            [b"\xf8\xfc", b"a", b"b", b"c"],
            [b"A\xf8\xfc", b"", b"a", b"b", b"c"],
            [np.nan, b"", b"b", b"c"],
            [b"A\xf8\xfc", np.nan, b"", b"b", b"c"],
        ]

        values = [
            [x.decode("latin-1") if isinstance(x, bytes) else x for x in y]
            for y in values
        ]

        examples = []
        for dtype in ["category", object]:
            for val in values:
                examples.append(Series(val, dtype=dtype))

        def roundtrip(s, encoding="latin-1"):
            with ensure_clean("test.json") as path:
                s.to_json(path, encoding=encoding)
                retr = read_json(path, encoding=encoding)
                assert_series_equal(s, retr, check_categorical=False)

        for s in examples:
            roundtrip(s)

    def test_data_frame_size_after_to_json(self):
        # GH15344
        df = DataFrame({"a": [str(1)]})

        size_before = df.memory_usage(index=True, deep=True).sum()
        df.to_json()
        size_after = df.memory_usage(index=True, deep=True).sum()

        assert size_before == size_after

    @pytest.mark.parametrize(
        "index", [None, [1, 2], [1.0, 2.0], ["a", "b"], ["1", "2"], ["1.", "2."]]
    )
    @pytest.mark.parametrize("columns", [["a", "b"], ["1", "2"], ["1.", "2."]])
    def test_from_json_to_json_table_index_and_columns(self, index, columns):
        # GH25433 GH25435
        expected = DataFrame([[1, 2], [3, 4]], index=index, columns=columns)
        dfjson = expected.to_json(orient="table")
        result = pd.read_json(dfjson, orient="table")
        assert_frame_equal(result, expected)

    def test_from_json_to_json_table_dtypes(self):
        # GH21345
        expected = pd.DataFrame({"a": [1, 2], "b": [3.0, 4.0], "c": ["5", "6"]})
        dfjson = expected.to_json(orient="table")
        result = pd.read_json(dfjson, orient="table")
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize("dtype", [True, {"b": int, "c": int}])
    def test_read_json_table_dtype_raises(self, dtype):
        # GH21345
        df = pd.DataFrame({"a": [1, 2], "b": [3.0, 4.0], "c": ["5", "6"]})
        dfjson = df.to_json(orient="table")
        msg = "cannot pass both dtype and orient='table'"
        with pytest.raises(ValueError, match=msg):
            pd.read_json(dfjson, orient="table", dtype=dtype)

    def test_read_json_table_convert_axes_raises(self):
        # GH25433 GH25435
        df = DataFrame([[1, 2], [3, 4]], index=[1.0, 2.0], columns=["1.", "2."])
        dfjson = df.to_json(orient="table")
        msg = "cannot pass both convert_axes and orient='table'"
        with pytest.raises(ValueError, match=msg):
            pd.read_json(dfjson, orient="table", convert_axes=True)

    @pytest.mark.parametrize(
        "data, expected",
        [
            (
                DataFrame([[1, 2], [4, 5]], columns=["a", "b"]),
                {"columns": ["a", "b"], "data": [[1, 2], [4, 5]]},
            ),
            (
                DataFrame([[1, 2], [4, 5]], columns=["a", "b"]).rename_axis("foo"),
                {"columns": ["a", "b"], "data": [[1, 2], [4, 5]]},
            ),
            (
                DataFrame(
                    [[1, 2], [4, 5]], columns=["a", "b"], index=[["a", "b"], ["c", "d"]]
                ),
                {"columns": ["a", "b"], "data": [[1, 2], [4, 5]]},
            ),
            (Series([1, 2, 3], name="A"), {"name": "A", "data": [1, 2, 3]}),
            (
                Series([1, 2, 3], name="A").rename_axis("foo"),
                {"name": "A", "data": [1, 2, 3]},
            ),
            (
                Series([1, 2], name="A", index=[["a", "b"], ["c", "d"]]),
                {"name": "A", "data": [1, 2]},
            ),
        ],
    )
    def test_index_false_to_json_split(self, data, expected):
        # GH 17394
        # Testing index=False in to_json with orient='split'

        result = data.to_json(orient="split", index=False)
        result = json.loads(result)

        assert result == expected

    @pytest.mark.parametrize(
        "data",
        [
            (DataFrame([[1, 2], [4, 5]], columns=["a", "b"])),
            (DataFrame([[1, 2], [4, 5]], columns=["a", "b"]).rename_axis("foo")),
            (
                DataFrame(
                    [[1, 2], [4, 5]], columns=["a", "b"], index=[["a", "b"], ["c", "d"]]
                )
            ),
            (Series([1, 2, 3], name="A")),
            (Series([1, 2, 3], name="A").rename_axis("foo")),
            (Series([1, 2], name="A", index=[["a", "b"], ["c", "d"]])),
        ],
    )
    def test_index_false_to_json_table(self, data):
        # GH 17394
        # Testing index=False in to_json with orient='table'

        result = data.to_json(orient="table", index=False)
        result = json.loads(result)

        expected = {
            "schema": pd.io.json.build_table_schema(data, index=False),
            "data": DataFrame(data).to_dict(orient="records"),
        }

        assert result == expected

    @pytest.mark.parametrize("orient", ["records", "index", "columns", "values"])
    def test_index_false_error_to_json(self, orient):
        # GH 17394
        # Testing error message from to_json with index=False

        df = pd.DataFrame([[1, 2], [4, 5]], columns=["a", "b"])

        msg = "'index=False' is only valid when 'orient' is 'split' or 'table'"
        with pytest.raises(ValueError, match=msg):
            df.to_json(orient=orient, index=False)

    @pytest.mark.parametrize("orient", ["split", "table"])
    @pytest.mark.parametrize("index", [True, False])
    def test_index_false_from_json_to_json(self, orient, index):
        # GH25170
        # Test index=False in from_json to_json
        expected = DataFrame({"a": [1, 2], "b": [3, 4]})
        dfjson = expected.to_json(orient=orient, index=index)
        result = read_json(dfjson, orient=orient)
        assert_frame_equal(result, expected)

    def test_read_timezone_information(self):
        # GH 25546
        result = read_json(
            '{"2019-01-01T11:00:00.000Z":88}', typ="series", orient="index"
        )
        expected = Series([88], index=DatetimeIndex(["2019-01-01 11:00:00"], tz="UTC"))
        assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "date_format,key", [("epoch", 86400000), ("iso", "P1DT0H0M0S")]
    )
    def test_timedelta_as_label(self, date_format, key):
        df = pd.DataFrame([[1]], columns=[pd.Timedelta("1D")])
        expected = '{{"{key}":{{"0":1}}}}'.format(key=key)
        result = df.to_json(date_format=date_format)

        assert result == expected

    @pytest.mark.parametrize(
        "orient,expected",
        [
            ("index", "{\"('a', 'b')\":{\"('c', 'd')\":1}}"),
            ("columns", "{\"('c', 'd')\":{\"('a', 'b')\":1}}"),
            # TODO: the below have separate encoding procedures
            # They produce JSON but not in a consistent manner
            pytest.param("split", "", marks=pytest.mark.skip),
            pytest.param("table", "", marks=pytest.mark.skip),
        ],
    )
    def test_tuple_labels(self, orient, expected):
        # GH 20500
        df = pd.DataFrame([[1]], index=[("a", "b")], columns=[("c", "d")])
        result = df.to_json(orient=orient)
        assert result == expected

    @pytest.mark.parametrize("indent", [1, 2, 4])
    def test_to_json_indent(self, indent):
        # GH 12004
        df = pd.DataFrame([["foo", "bar"], ["baz", "qux"]], columns=["a", "b"])

        result = df.to_json(indent=indent)
        spaces = " " * indent
        expected = """{{
{spaces}"a":{{
{spaces}{spaces}"0":"foo",
{spaces}{spaces}"1":"baz"
{spaces}}},
{spaces}"b":{{
{spaces}{spaces}"0":"bar",
{spaces}{spaces}"1":"qux"
{spaces}}}
}}""".format(
            spaces=spaces
        )

        assert result == expected

    @pytest.mark.parametrize(
        "orient,expected",
        [
            (
                "split",
                """{
    "columns":[
        "a",
        "b"
    ],
    "index":[
        0,
        1
    ],
    "data":[
        [
            "foo",
            "bar"
        ],
        [
            "baz",
            "qux"
        ]
    ]
}""",
            ),
            (
                "records",
                """[
    {
        "a":"foo",
        "b":"bar"
    },
    {
        "a":"baz",
        "b":"qux"
    }
]""",
            ),
            (
                "index",
                """{
    "0":{
        "a":"foo",
        "b":"bar"
    },
    "1":{
        "a":"baz",
        "b":"qux"
    }
}""",
            ),
            (
                "columns",
                """{
    "a":{
        "0":"foo",
        "1":"baz"
    },
    "b":{
        "0":"bar",
        "1":"qux"
    }
}""",
            ),
            (
                "values",
                """[
    [
        "foo",
        "bar"
    ],
    [
        "baz",
        "qux"
    ]
]""",
            ),
            (
                "table",
                """{
    "schema":{
        "fields":[
            {
                "name":"index",
                "type":"integer"
            },
            {
                "name":"a",
                "type":"string"
            },
            {
                "name":"b",
                "type":"string"
            }
        ],
        "primaryKey":[
            "index"
        ],
        "pandas_version":"0.20.0"
    },
    "data":[
        {
            "index":0,
            "a":"foo",
            "b":"bar"
        },
        {
            "index":1,
            "a":"baz",
            "b":"qux"
        }
    ]
}""",
            ),
        ],
    )
    def test_json_indent_all_orients(self, orient, expected):
        # GH 12004
        df = pd.DataFrame([["foo", "bar"], ["baz", "qux"]], columns=["a", "b"])
        result = df.to_json(orient=orient, indent=4)

        if PY35:
            assert json.loads(result) == json.loads(expected)
        else:
            assert result == expected

    def test_json_negative_indent_raises(self):
        with pytest.raises(ValueError, match="must be a nonnegative integer"):
            pd.DataFrame().to_json(indent=-1)
