from collections import OrderedDict, abc
from datetime import date, datetime, timedelta
import functools
import itertools

import numpy as np
import numpy.ma as ma
import numpy.ma.mrecords as mrecords
import pytest

from pandas.compat import is_platform_little_endian

from pandas.core.dtypes.common import is_integer_dtype

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    Index,
    MultiIndex,
    RangeIndex,
    Series,
    Timedelta,
    Timestamp,
    date_range,
    isna,
)
import pandas._testing as tm
from pandas.arrays import IntervalArray, PeriodArray, SparseArray
from pandas.core.construction import create_series_with_explicit_dtype

MIXED_FLOAT_DTYPES = ["float16", "float32", "float64"]
MIXED_INT_DTYPES = [
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    "int8",
    "int16",
    "int32",
    "int64",
]


class TestDataFrameConstructors:
    def test_series_with_name_not_matching_column(self):
        # GH#9232
        x = pd.Series(range(5), name=1)
        y = pd.Series(range(5), name=0)

        result = pd.DataFrame(x, columns=[0])
        expected = pd.DataFrame([], columns=[0])
        tm.assert_frame_equal(result, expected)

        result = pd.DataFrame(y, columns=[1])
        expected = pd.DataFrame([], columns=[1])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "constructor",
        [
            lambda: DataFrame(),
            lambda: DataFrame(None),
            lambda: DataFrame({}),
            lambda: DataFrame(()),
            lambda: DataFrame([]),
            lambda: DataFrame((_ for _ in [])),
            lambda: DataFrame(range(0)),
            lambda: DataFrame(data=None),
            lambda: DataFrame(data={}),
            lambda: DataFrame(data=()),
            lambda: DataFrame(data=[]),
            lambda: DataFrame(data=(_ for _ in [])),
            lambda: DataFrame(data=range(0)),
        ],
    )
    def test_empty_constructor(self, constructor):
        expected = DataFrame()
        result = constructor()
        assert len(result.index) == 0
        assert len(result.columns) == 0
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "emptylike,expected_index,expected_columns",
        [
            ([[]], RangeIndex(1), RangeIndex(0)),
            ([[], []], RangeIndex(2), RangeIndex(0)),
            ([(_ for _ in [])], RangeIndex(1), RangeIndex(0)),
        ],
    )
    def test_emptylike_constructor(self, emptylike, expected_index, expected_columns):
        expected = DataFrame(index=expected_index, columns=expected_columns)
        result = DataFrame(emptylike)
        tm.assert_frame_equal(result, expected)

    def test_constructor_mixed(self, float_string_frame):
        index, data = tm.getMixedTypeDict()

        # TODO(wesm), incomplete test?
        indexed_frame = DataFrame(data, index=index)  # noqa
        unindexed_frame = DataFrame(data)  # noqa

        assert float_string_frame["foo"].dtype == np.object_

    def test_constructor_cast_failure(self):
        foo = DataFrame({"a": ["a", "b", "c"]}, dtype=np.float64)
        assert foo["a"].dtype == object

        # GH 3010, constructing with odd arrays
        df = DataFrame(np.ones((4, 2)))

        # this is ok
        df["foo"] = np.ones((4, 2)).tolist()

        # this is not ok
        msg = "Wrong number of items passed 2, placement implies 1"
        with pytest.raises(ValueError, match=msg):
            df["test"] = np.ones((4, 2))

        # this is ok
        df["foo2"] = np.ones((4, 2)).tolist()

    def test_constructor_dtype_copy(self):
        orig_df = DataFrame({"col1": [1.0], "col2": [2.0], "col3": [3.0]})

        new_df = pd.DataFrame(orig_df, dtype=float, copy=True)

        new_df["col1"] = 200.0
        assert orig_df["col1"][0] == 1.0

    def test_constructor_dtype_nocast_view(self):
        df = DataFrame([[1, 2]])
        should_be_view = DataFrame(df, dtype=df[0].dtype)
        should_be_view[0][0] = 99
        assert df.values[0, 0] == 99

        should_be_view = DataFrame(df.values, dtype=df[0].dtype)
        should_be_view[0][0] = 97
        assert df.values[0, 0] == 97

    def test_constructor_dtype_list_data(self):
        df = DataFrame([[1, "2"], [None, "a"]], dtype=object)
        assert df.loc[1, 0] is None
        assert df.loc[0, 1] == "2"

    def test_constructor_list_frames(self):
        # see gh-3243
        result = DataFrame([DataFrame()])
        assert result.shape == (1, 0)

        result = DataFrame([DataFrame(dict(A=np.arange(5)))])
        assert isinstance(result.iloc[0, 0], DataFrame)

    def test_constructor_mixed_dtypes(self):
        def _make_mixed_dtypes_df(typ, ad=None):

            if typ == "int":
                dtypes = MIXED_INT_DTYPES
                arrays = [np.array(np.random.rand(10), dtype=d) for d in dtypes]
            elif typ == "float":
                dtypes = MIXED_FLOAT_DTYPES
                arrays = [
                    np.array(np.random.randint(10, size=10), dtype=d) for d in dtypes
                ]

            for d, a in zip(dtypes, arrays):
                assert a.dtype == d
            if ad is None:
                ad = dict()
            ad.update({d: a for d, a in zip(dtypes, arrays)})
            return DataFrame(ad)

        def _check_mixed_dtypes(df, dtypes=None):
            if dtypes is None:
                dtypes = MIXED_FLOAT_DTYPES + MIXED_INT_DTYPES
            for d in dtypes:
                if d in df:
                    assert df.dtypes[d] == d

        # mixed floating and integer coexist in the same frame
        df = _make_mixed_dtypes_df("float")
        _check_mixed_dtypes(df)

        # add lots of types
        df = _make_mixed_dtypes_df("float", dict(A=1, B="foo", C="bar"))
        _check_mixed_dtypes(df)

        # GH 622
        df = _make_mixed_dtypes_df("int")
        _check_mixed_dtypes(df)

    def test_constructor_complex_dtypes(self):
        # GH10952
        a = np.random.rand(10).astype(np.complex64)
        b = np.random.rand(10).astype(np.complex128)

        df = DataFrame({"a": a, "b": b})
        assert a.dtype == df.a.dtype
        assert b.dtype == df.b.dtype

    def test_constructor_dtype_str_na_values(self, string_dtype):
        # https://github.com/pandas-dev/pandas/issues/21083
        df = DataFrame({"A": ["x", None]}, dtype=string_dtype)
        result = df.isna()
        expected = DataFrame({"A": [False, True]})
        tm.assert_frame_equal(result, expected)
        assert df.iloc[1, 0] is None

        df = DataFrame({"A": ["x", np.nan]}, dtype=string_dtype)
        assert np.isnan(df.iloc[1, 0])

    def test_constructor_rec(self, float_frame):
        rec = float_frame.to_records(index=False)
        rec.dtype.names = list(rec.dtype.names)[::-1]

        index = float_frame.index

        df = DataFrame(rec)
        tm.assert_index_equal(df.columns, pd.Index(rec.dtype.names))

        df2 = DataFrame(rec, index=index)
        tm.assert_index_equal(df2.columns, pd.Index(rec.dtype.names))
        tm.assert_index_equal(df2.index, index)

        rng = np.arange(len(rec))[::-1]
        df3 = DataFrame(rec, index=rng, columns=["C", "B"])
        expected = DataFrame(rec, index=rng).reindex(columns=["C", "B"])
        tm.assert_frame_equal(df3, expected)

    def test_constructor_bool(self):
        df = DataFrame({0: np.ones(10, dtype=bool), 1: np.zeros(10, dtype=bool)})
        assert df.values.dtype == np.bool_

    def test_constructor_overflow_int64(self):
        # see gh-14881
        values = np.array([2 ** 64 - i for i in range(1, 10)], dtype=np.uint64)

        result = DataFrame({"a": values})
        assert result["a"].dtype == np.uint64

        # see gh-2355
        data_scores = [
            (6311132704823138710, 273),
            (2685045978526272070, 23),
            (8921811264899370420, 45),
            (17019687244989530680, 270),
            (9930107427299601010, 273),
        ]
        dtype = [("uid", "u8"), ("score", "u8")]
        data = np.zeros((len(data_scores),), dtype=dtype)
        data[:] = data_scores
        df_crawls = DataFrame(data)
        assert df_crawls["uid"].dtype == np.uint64

    @pytest.mark.parametrize(
        "values",
        [
            np.array([2 ** 64], dtype=object),
            np.array([2 ** 65]),
            [2 ** 64 + 1],
            np.array([-(2 ** 63) - 4], dtype=object),
            np.array([-(2 ** 64) - 1]),
            [-(2 ** 65) - 2],
        ],
    )
    def test_constructor_int_overflow(self, values):
        # see gh-18584
        value = values[0]
        result = DataFrame(values)

        assert result[0].dtype == object
        assert result[0][0] == value

    def test_constructor_ordereddict(self):
        import random

        nitems = 100
        nums = list(range(nitems))
        random.shuffle(nums)
        expected = ["A{i:d}".format(i=i) for i in nums]
        df = DataFrame(OrderedDict(zip(expected, [[0]] * nitems)))
        assert expected == list(df.columns)

    def test_constructor_dict(self):
        datetime_series = tm.makeTimeSeries(nper=30)
        # test expects index shifted by 5
        datetime_series_short = tm.makeTimeSeries(nper=30)[5:]

        frame = DataFrame({"col1": datetime_series, "col2": datetime_series_short})

        # col2 is padded with NaN
        assert len(datetime_series) == 30
        assert len(datetime_series_short) == 25

        tm.assert_series_equal(frame["col1"], datetime_series.rename("col1"))

        exp = pd.Series(
            np.concatenate([[np.nan] * 5, datetime_series_short.values]),
            index=datetime_series.index,
            name="col2",
        )
        tm.assert_series_equal(exp, frame["col2"])

        frame = DataFrame(
            {"col1": datetime_series, "col2": datetime_series_short},
            columns=["col2", "col3", "col4"],
        )

        assert len(frame) == len(datetime_series_short)
        assert "col1" not in frame
        assert isna(frame["col3"]).all()

        # Corner cases
        assert len(DataFrame()) == 0

        # mix dict and array, wrong size - no spec for which error should raise
        # first
        with pytest.raises(ValueError):
            DataFrame({"A": {"a": "a", "b": "b"}, "B": ["a", "b", "c"]})

        # Length-one dict micro-optimization
        frame = DataFrame({"A": {"1": 1, "2": 2}})
        tm.assert_index_equal(frame.index, pd.Index(["1", "2"]))

        # empty dict plus index
        idx = Index([0, 1, 2])
        frame = DataFrame({}, index=idx)
        assert frame.index is idx

        # empty dict with index and columns
        idx = Index([0, 1, 2])
        frame = DataFrame({}, index=idx, columns=idx)
        assert frame.index is idx
        assert frame.columns is idx
        assert len(frame._series) == 3

        # with dict of empty list and Series
        frame = DataFrame({"A": [], "B": []}, columns=["A", "B"])
        tm.assert_index_equal(frame.index, Index([], dtype=np.int64))

        # GH 14381
        # Dict with None value
        frame_none = DataFrame(dict(a=None), index=[0])
        frame_none_list = DataFrame(dict(a=[None]), index=[0])
        assert frame_none._get_value(0, "a") is None
        assert frame_none_list._get_value(0, "a") is None
        tm.assert_frame_equal(frame_none, frame_none_list)

        # GH10856
        # dict with scalar values should raise error, even if columns passed
        msg = "If using all scalar values, you must pass an index"
        with pytest.raises(ValueError, match=msg):
            DataFrame({"a": 0.7})

        with pytest.raises(ValueError, match=msg):
            DataFrame({"a": 0.7}, columns=["a"])

    @pytest.mark.parametrize("scalar", [2, np.nan, None, "D"])
    def test_constructor_invalid_items_unused(self, scalar):
        # No error if invalid (scalar) value is in fact not used:
        result = DataFrame({"a": scalar}, columns=["b"])
        expected = DataFrame(columns=["b"])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("value", [2, np.nan, None, float("nan")])
    def test_constructor_dict_nan_key(self, value):
        # GH 18455
        cols = [1, value, 3]
        idx = ["a", value]
        values = [[0, 3], [1, 4], [2, 5]]
        data = {cols[c]: Series(values[c], index=idx) for c in range(3)}
        result = DataFrame(data).sort_values(1).sort_values("a", axis=1)
        expected = DataFrame(
            np.arange(6, dtype="int64").reshape(2, 3), index=idx, columns=cols
        )
        tm.assert_frame_equal(result, expected)

        result = DataFrame(data, index=idx).sort_values("a", axis=1)
        tm.assert_frame_equal(result, expected)

        result = DataFrame(data, index=idx, columns=cols)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("value", [np.nan, None, float("nan")])
    def test_constructor_dict_nan_tuple_key(self, value):
        # GH 18455
        cols = Index([(11, 21), (value, 22), (13, value)])
        idx = Index([("a", value), (value, 2)])
        values = [[0, 3], [1, 4], [2, 5]]
        data = {cols[c]: Series(values[c], index=idx) for c in range(3)}
        result = DataFrame(data).sort_values((11, 21)).sort_values(("a", value), axis=1)
        expected = DataFrame(
            np.arange(6, dtype="int64").reshape(2, 3), index=idx, columns=cols
        )
        tm.assert_frame_equal(result, expected)

        result = DataFrame(data, index=idx).sort_values(("a", value), axis=1)
        tm.assert_frame_equal(result, expected)

        result = DataFrame(data, index=idx, columns=cols)
        tm.assert_frame_equal(result, expected)

    def test_constructor_dict_order_insertion(self):
        datetime_series = tm.makeTimeSeries(nper=30)
        datetime_series_short = tm.makeTimeSeries(nper=25)

        # GH19018
        # initialization ordering: by insertion order if python>= 3.6
        d = {"b": datetime_series_short, "a": datetime_series}
        frame = DataFrame(data=d)
        expected = DataFrame(data=d, columns=list("ba"))
        tm.assert_frame_equal(frame, expected)

    def test_constructor_multi_index(self):
        # GH 4078
        # construction error with mi and all-nan frame
        tuples = [(2, 3), (3, 3), (3, 3)]
        mi = MultiIndex.from_tuples(tuples)
        df = DataFrame(index=mi, columns=mi)
        assert pd.isna(df).values.ravel().all()

        tuples = [(3, 3), (2, 3), (3, 3)]
        mi = MultiIndex.from_tuples(tuples)
        df = DataFrame(index=mi, columns=mi)
        assert pd.isna(df).values.ravel().all()

    def test_constructor_2d_index(self):
        # GH 25416
        # handling of 2d index in construction
        df = pd.DataFrame([[1]], columns=[[1]], index=[1, 2])
        expected = pd.DataFrame(
            [1, 1],
            index=pd.Int64Index([1, 2], dtype="int64"),
            columns=pd.MultiIndex(levels=[[1]], codes=[[0]]),
        )
        tm.assert_frame_equal(df, expected)

        df = pd.DataFrame([[1]], columns=[[1]], index=[[1, 2]])
        expected = pd.DataFrame(
            [1, 1],
            index=pd.MultiIndex(levels=[[1, 2]], codes=[[0, 1]]),
            columns=pd.MultiIndex(levels=[[1]], codes=[[0]]),
        )
        tm.assert_frame_equal(df, expected)

    def test_constructor_error_msgs(self):
        msg = "Empty data passed with indices specified."
        # passing an empty array with columns specified.
        with pytest.raises(ValueError, match=msg):
            DataFrame(np.empty(0), columns=list("abc"))

        msg = "Mixing dicts with non-Series may lead to ambiguous ordering."
        # mix dict and array, wrong size
        with pytest.raises(ValueError, match=msg):
            DataFrame({"A": {"a": "a", "b": "b"}, "B": ["a", "b", "c"]})

        # wrong size ndarray, GH 3105
        msg = r"Shape of passed values is \(4, 3\), indices imply \(3, 3\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(
                np.arange(12).reshape((4, 3)),
                columns=["foo", "bar", "baz"],
                index=pd.date_range("2000-01-01", periods=3),
            )

        arr = np.array([[4, 5, 6]])
        msg = r"Shape of passed values is \(1, 3\), indices imply \(1, 4\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(index=[0], columns=range(0, 4), data=arr)

        arr = np.array([4, 5, 6])
        msg = r"Shape of passed values is \(3, 1\), indices imply \(1, 4\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(index=[0], columns=range(0, 4), data=arr)

        # higher dim raise exception
        with pytest.raises(ValueError, match="Must pass 2-d input"):
            DataFrame(np.zeros((3, 3, 3)), columns=["A", "B", "C"], index=[1])

        # wrong size axis labels
        msg = r"Shape of passed values is \(2, 3\), indices imply \(1, 3\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(np.random.rand(2, 3), columns=["A", "B", "C"], index=[1])

        msg = r"Shape of passed values is \(2, 3\), indices imply \(2, 2\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(np.random.rand(2, 3), columns=["A", "B"], index=[1, 2])

        # gh-26429
        msg = "2 columns passed, passed data had 10 columns"
        with pytest.raises(ValueError, match=msg):
            DataFrame((range(10), range(10, 20)), columns=("ones", "twos"))

        msg = "If using all scalar values, you must pass an index"
        with pytest.raises(ValueError, match=msg):
            DataFrame({"a": False, "b": True})

    def test_constructor_with_embedded_frames(self):

        # embedded data frames
        df1 = DataFrame({"a": [1, 2, 3], "b": [3, 4, 5]})
        df2 = DataFrame([df1, df1 + 10])

        df2.dtypes
        str(df2)

        result = df2.loc[0, 0]
        tm.assert_frame_equal(result, df1)

        result = df2.loc[1, 0]
        tm.assert_frame_equal(result, df1 + 10)

    def test_constructor_subclass_dict(self, float_frame, dict_subclass):
        # Test for passing dict subclass to constructor
        data = {
            "col1": dict_subclass((x, 10.0 * x) for x in range(10)),
            "col2": dict_subclass((x, 20.0 * x) for x in range(10)),
        }
        df = DataFrame(data)
        refdf = DataFrame({col: dict(val.items()) for col, val in data.items()})
        tm.assert_frame_equal(refdf, df)

        data = dict_subclass(data.items())
        df = DataFrame(data)
        tm.assert_frame_equal(refdf, df)

        # try with defaultdict
        from collections import defaultdict

        data = {}
        float_frame["B"][:10] = np.nan
        for k, v in float_frame.items():
            dct = defaultdict(dict)
            dct.update(v.to_dict())
            data[k] = dct
        frame = DataFrame(data)
        expected = frame.reindex(index=float_frame.index)
        tm.assert_frame_equal(float_frame, expected)

    def test_constructor_dict_block(self):
        expected = np.array([[4.0, 3.0, 2.0, 1.0]])
        df = DataFrame(
            {"d": [4.0], "c": [3.0], "b": [2.0], "a": [1.0]},
            columns=["d", "c", "b", "a"],
        )
        tm.assert_numpy_array_equal(df.values, expected)

    def test_constructor_dict_cast(self):
        # cast float tests
        test_data = {"A": {"1": 1, "2": 2}, "B": {"1": "1", "2": "2", "3": "3"}}
        frame = DataFrame(test_data, dtype=float)
        assert len(frame) == 3
        assert frame["B"].dtype == np.float64
        assert frame["A"].dtype == np.float64

        frame = DataFrame(test_data)
        assert len(frame) == 3
        assert frame["B"].dtype == np.object_
        assert frame["A"].dtype == np.float64

        # can't cast to float
        test_data = {
            "A": dict(zip(range(20), tm.makeStringIndex(20))),
            "B": dict(zip(range(15), np.random.randn(15))),
        }
        frame = DataFrame(test_data, dtype=float)
        assert len(frame) == 20
        assert frame["A"].dtype == np.object_
        assert frame["B"].dtype == np.float64

    def test_constructor_dict_dont_upcast(self):
        d = {"Col1": {"Row1": "A String", "Row2": np.nan}}
        df = DataFrame(d)
        assert isinstance(df["Col1"]["Row2"], float)

        dm = DataFrame([[1, 2], ["a", "b"]], index=[1, 2], columns=[1, 2])
        assert isinstance(dm[1][1], int)

    def test_constructor_dict_of_tuples(self):
        # GH #1491
        data = {"a": (1, 2, 3), "b": (4, 5, 6)}

        result = DataFrame(data)
        expected = DataFrame({k: list(v) for k, v in data.items()})
        tm.assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_dict_of_ranges(self):
        # GH 26356
        data = {"a": range(3), "b": range(3, 6)}

        result = DataFrame(data)
        expected = DataFrame({"a": [0, 1, 2], "b": [3, 4, 5]})
        tm.assert_frame_equal(result, expected)

    def test_constructor_dict_of_iterators(self):
        # GH 26349
        data = {"a": iter(range(3)), "b": reversed(range(3))}

        result = DataFrame(data)
        expected = DataFrame({"a": [0, 1, 2], "b": [2, 1, 0]})
        tm.assert_frame_equal(result, expected)

    def test_constructor_dict_of_generators(self):
        # GH 26349
        data = {"a": (i for i in (range(3))), "b": (i for i in reversed(range(3)))}
        result = DataFrame(data)
        expected = DataFrame({"a": [0, 1, 2], "b": [2, 1, 0]})
        tm.assert_frame_equal(result, expected)

    def test_constructor_dict_multiindex(self):
        def check(result, expected):
            return tm.assert_frame_equal(
                result,
                expected,
                check_dtype=True,
                check_index_type=True,
                check_column_type=True,
                check_names=True,
            )

        d = {
            ("a", "a"): {("i", "i"): 0, ("i", "j"): 1, ("j", "i"): 2},
            ("b", "a"): {("i", "i"): 6, ("i", "j"): 5, ("j", "i"): 4},
            ("b", "c"): {("i", "i"): 7, ("i", "j"): 8, ("j", "i"): 9},
        }
        _d = sorted(d.items())
        df = DataFrame(d)
        expected = DataFrame(
            [x[1] for x in _d], index=MultiIndex.from_tuples([x[0] for x in _d])
        ).T
        expected.index = MultiIndex.from_tuples(expected.index)
        check(df, expected)

        d["z"] = {"y": 123.0, ("i", "i"): 111, ("i", "j"): 111, ("j", "i"): 111}
        _d.insert(0, ("z", d["z"]))
        expected = DataFrame(
            [x[1] for x in _d], index=Index([x[0] for x in _d], tupleize_cols=False)
        ).T
        expected.index = Index(expected.index, tupleize_cols=False)
        df = DataFrame(d)
        df = df.reindex(columns=expected.columns, index=expected.index)
        check(df, expected)

    def test_constructor_dict_datetime64_index(self):
        # GH 10160
        dates_as_str = ["1984-02-19", "1988-11-06", "1989-12-03", "1990-03-15"]

        def create_data(constructor):
            return {i: {constructor(s): 2 * i} for i, s in enumerate(dates_as_str)}

        data_datetime64 = create_data(np.datetime64)
        data_datetime = create_data(lambda x: datetime.strptime(x, "%Y-%m-%d"))
        data_Timestamp = create_data(Timestamp)

        expected = DataFrame(
            [
                {0: 0, 1: None, 2: None, 3: None},
                {0: None, 1: 2, 2: None, 3: None},
                {0: None, 1: None, 2: 4, 3: None},
                {0: None, 1: None, 2: None, 3: 6},
            ],
            index=[Timestamp(dt) for dt in dates_as_str],
        )

        result_datetime64 = DataFrame(data_datetime64)
        result_datetime = DataFrame(data_datetime)
        result_Timestamp = DataFrame(data_Timestamp)
        tm.assert_frame_equal(result_datetime64, expected)
        tm.assert_frame_equal(result_datetime, expected)
        tm.assert_frame_equal(result_Timestamp, expected)

    def test_constructor_dict_timedelta64_index(self):
        # GH 10160
        td_as_int = [1, 2, 3, 4]

        def create_data(constructor):
            return {i: {constructor(s): 2 * i} for i, s in enumerate(td_as_int)}

        data_timedelta64 = create_data(lambda x: np.timedelta64(x, "D"))
        data_timedelta = create_data(lambda x: timedelta(days=x))
        data_Timedelta = create_data(lambda x: Timedelta(x, "D"))

        expected = DataFrame(
            [
                {0: 0, 1: None, 2: None, 3: None},
                {0: None, 1: 2, 2: None, 3: None},
                {0: None, 1: None, 2: 4, 3: None},
                {0: None, 1: None, 2: None, 3: 6},
            ],
            index=[Timedelta(td, "D") for td in td_as_int],
        )

        result_timedelta64 = DataFrame(data_timedelta64)
        result_timedelta = DataFrame(data_timedelta)
        result_Timedelta = DataFrame(data_Timedelta)
        tm.assert_frame_equal(result_timedelta64, expected)
        tm.assert_frame_equal(result_timedelta, expected)
        tm.assert_frame_equal(result_Timedelta, expected)

    def test_constructor_period(self):
        # PeriodIndex
        a = pd.PeriodIndex(["2012-01", "NaT", "2012-04"], freq="M")
        b = pd.PeriodIndex(["2012-02-01", "2012-03-01", "NaT"], freq="D")
        df = pd.DataFrame({"a": a, "b": b})
        assert df["a"].dtype == a.dtype
        assert df["b"].dtype == b.dtype

        # list of periods
        df = pd.DataFrame(
            {"a": a.astype(object).tolist(), "b": b.astype(object).tolist()}
        )
        assert df["a"].dtype == a.dtype
        assert df["b"].dtype == b.dtype

    def test_nested_dict_frame_constructor(self):
        rng = pd.period_range("1/1/2000", periods=5)
        df = DataFrame(np.random.randn(10, 5), columns=rng)

        data = {}
        for col in df.columns:
            for row in df.index:
                data.setdefault(col, {})[row] = df._get_value(row, col)

        result = DataFrame(data, columns=rng)
        tm.assert_frame_equal(result, df)

        data = {}
        for col in df.columns:
            for row in df.index:
                data.setdefault(row, {})[col] = df._get_value(row, col)

        result = DataFrame(data, index=rng).T
        tm.assert_frame_equal(result, df)

    def _check_basic_constructor(self, empty):
        # mat: 2d matrix with shape (3, 2) to input. empty - makes sized
        # objects
        mat = empty((2, 3), dtype=float)
        # 2-D input
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2])

        assert len(frame.index) == 2
        assert len(frame.columns) == 3

        # 1-D input
        frame = DataFrame(empty((3,)), columns=["A"], index=[1, 2, 3])
        assert len(frame.index) == 3
        assert len(frame.columns) == 1

        # cast type
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2], dtype=np.int64)
        assert frame.values.dtype == np.int64

        # wrong size axis labels
        msg = r"Shape of passed values is \(2, 3\), indices imply \(1, 3\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(mat, columns=["A", "B", "C"], index=[1])
        msg = r"Shape of passed values is \(2, 3\), indices imply \(2, 2\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(mat, columns=["A", "B"], index=[1, 2])

        # higher dim raise exception
        with pytest.raises(ValueError, match="Must pass 2-d input"):
            DataFrame(empty((3, 3, 3)), columns=["A", "B", "C"], index=[1])

        # automatic labeling
        frame = DataFrame(mat)
        tm.assert_index_equal(frame.index, pd.Int64Index(range(2)))
        tm.assert_index_equal(frame.columns, pd.Int64Index(range(3)))

        frame = DataFrame(mat, index=[1, 2])
        tm.assert_index_equal(frame.columns, pd.Int64Index(range(3)))

        frame = DataFrame(mat, columns=["A", "B", "C"])
        tm.assert_index_equal(frame.index, pd.Int64Index(range(2)))

        # 0-length axis
        frame = DataFrame(empty((0, 3)))
        assert len(frame.index) == 0

        frame = DataFrame(empty((3, 0)))
        assert len(frame.columns) == 0

    def test_constructor_ndarray(self):
        self._check_basic_constructor(np.ones)

        frame = DataFrame(["foo", "bar"], index=[0, 1], columns=["A"])
        assert len(frame) == 2

    def test_constructor_maskedarray(self):
        self._check_basic_constructor(ma.masked_all)

        # Check non-masked values
        mat = ma.masked_all((2, 3), dtype=float)
        mat[0, 0] = 1.0
        mat[1, 2] = 2.0
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2])
        assert 1.0 == frame["A"][1]
        assert 2.0 == frame["C"][2]

        # what is this even checking??
        mat = ma.masked_all((2, 3), dtype=float)
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2])
        assert np.all(~np.asarray(frame == frame))

    def test_constructor_maskedarray_nonfloat(self):
        # masked int promoted to float
        mat = ma.masked_all((2, 3), dtype=int)
        # 2-D input
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2])

        assert len(frame.index) == 2
        assert len(frame.columns) == 3
        assert np.all(~np.asarray(frame == frame))

        # cast type
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2], dtype=np.float64)
        assert frame.values.dtype == np.float64

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1
        mat2[1, 2] = 2
        frame = DataFrame(mat2, columns=["A", "B", "C"], index=[1, 2])
        assert 1 == frame["A"][1]
        assert 2 == frame["C"][2]

        # masked np.datetime64 stays (use NaT as null)
        mat = ma.masked_all((2, 3), dtype="M8[ns]")
        # 2-D input
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2])

        assert len(frame.index) == 2
        assert len(frame.columns) == 3
        assert isna(frame).values.all()

        # cast type
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2], dtype=np.int64)
        assert frame.values.dtype == np.int64

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1
        mat2[1, 2] = 2
        frame = DataFrame(mat2, columns=["A", "B", "C"], index=[1, 2])
        assert 1 == frame["A"].view("i8")[1]
        assert 2 == frame["C"].view("i8")[2]

        # masked bool promoted to object
        mat = ma.masked_all((2, 3), dtype=bool)
        # 2-D input
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2])

        assert len(frame.index) == 2
        assert len(frame.columns) == 3
        assert np.all(~np.asarray(frame == frame))

        # cast type
        frame = DataFrame(mat, columns=["A", "B", "C"], index=[1, 2], dtype=object)
        assert frame.values.dtype == object

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = True
        mat2[1, 2] = False
        frame = DataFrame(mat2, columns=["A", "B", "C"], index=[1, 2])
        assert frame["A"][1] is True
        assert frame["C"][2] is False

    def test_constructor_maskedarray_hardened(self):
        # Check numpy masked arrays with hard masks -- from GH24574
        mat_hard = ma.masked_all((2, 2), dtype=float).harden_mask()
        result = pd.DataFrame(mat_hard, columns=["A", "B"], index=[1, 2])
        expected = pd.DataFrame(
            {"A": [np.nan, np.nan], "B": [np.nan, np.nan]},
            columns=["A", "B"],
            index=[1, 2],
            dtype=float,
        )
        tm.assert_frame_equal(result, expected)
        # Check case where mask is hard but no data are masked
        mat_hard = ma.ones((2, 2), dtype=float).harden_mask()
        result = pd.DataFrame(mat_hard, columns=["A", "B"], index=[1, 2])
        expected = pd.DataFrame(
            {"A": [1.0, 1.0], "B": [1.0, 1.0]},
            columns=["A", "B"],
            index=[1, 2],
            dtype=float,
        )
        tm.assert_frame_equal(result, expected)

    def test_constructor_maskedrecarray_dtype(self):
        # Ensure constructor honors dtype
        data = np.ma.array(
            np.ma.zeros(5, dtype=[("date", "<f8"), ("price", "<f8")]), mask=[False] * 5
        )
        data = data.view(mrecords.mrecarray)
        result = pd.DataFrame(data, dtype=int)
        expected = pd.DataFrame(np.zeros((5, 2), dtype=int), columns=["date", "price"])
        tm.assert_frame_equal(result, expected)

    def test_constructor_mrecarray(self):
        # Ensure mrecarray produces frame identical to dict of masked arrays
        # from GH3479

        assert_fr_equal = functools.partial(
            tm.assert_frame_equal,
            check_index_type=True,
            check_column_type=True,
            check_frame_type=True,
        )
        arrays = [
            ("float", np.array([1.5, 2.0])),
            ("int", np.array([1, 2])),
            ("str", np.array(["abc", "def"])),
        ]
        for name, arr in arrays[:]:
            arrays.append(
                ("masked1_" + name, np.ma.masked_array(arr, mask=[False, True]))
            )
        arrays.append(("masked_all", np.ma.masked_all((2,))))
        arrays.append(("masked_none", np.ma.masked_array([1.0, 2.5], mask=False)))

        # call assert_frame_equal for all selections of 3 arrays
        for comb in itertools.combinations(arrays, 3):
            names, data = zip(*comb)
            mrecs = mrecords.fromarrays(data, names=names)

            # fill the comb
            comb = {k: (v.filled() if hasattr(v, "filled") else v) for k, v in comb}

            expected = DataFrame(comb, columns=names)
            result = DataFrame(mrecs)
            assert_fr_equal(result, expected)

            # specify columns
            expected = DataFrame(comb, columns=names[::-1])
            result = DataFrame(mrecs, columns=names[::-1])
            assert_fr_equal(result, expected)

            # specify index
            expected = DataFrame(comb, columns=names, index=[1, 2])
            result = DataFrame(mrecs, index=[1, 2])
            assert_fr_equal(result, expected)

    def test_constructor_corner_shape(self):
        df = DataFrame(index=[])
        assert df.values.shape == (0, 0)

    @pytest.mark.parametrize(
        "data, index, columns, dtype, expected",
        [
            (None, list(range(10)), ["a", "b"], object, np.object_),
            (None, None, ["a", "b"], "int64", np.dtype("int64")),
            (None, list(range(10)), ["a", "b"], int, np.dtype("float64")),
            ({}, None, ["foo", "bar"], None, np.object_),
            ({"b": 1}, list(range(10)), list("abc"), int, np.dtype("float64")),
        ],
    )
    def test_constructor_dtype(self, data, index, columns, dtype, expected):
        df = DataFrame(data, index, columns, dtype)
        assert df.values.dtype == expected

    def test_constructor_scalar_inference(self):
        data = {"int": 1, "bool": True, "float": 3.0, "complex": 4j, "object": "foo"}
        df = DataFrame(data, index=np.arange(10))

        assert df["int"].dtype == np.int64
        assert df["bool"].dtype == np.bool_
        assert df["float"].dtype == np.float64
        assert df["complex"].dtype == np.complex128
        assert df["object"].dtype == np.object_

    def test_constructor_arrays_and_scalars(self):
        df = DataFrame({"a": np.random.randn(10), "b": True})
        exp = DataFrame({"a": df["a"].values, "b": [True] * 10})

        tm.assert_frame_equal(df, exp)
        with pytest.raises(ValueError, match="must pass an index"):
            DataFrame({"a": False, "b": True})

    def test_constructor_DataFrame(self, float_frame):
        df = DataFrame(float_frame)
        tm.assert_frame_equal(df, float_frame)

        df_casted = DataFrame(float_frame, dtype=np.int64)
        assert df_casted.values.dtype == np.int64

    def test_constructor_more(self, float_frame):
        # used to be in test_matrix.py
        arr = np.random.randn(10)
        dm = DataFrame(arr, columns=["A"], index=np.arange(10))
        assert dm.values.ndim == 2

        arr = np.random.randn(0)
        dm = DataFrame(arr)
        assert dm.values.ndim == 2
        assert dm.values.ndim == 2

        # no data specified
        dm = DataFrame(columns=["A", "B"], index=np.arange(10))
        assert dm.values.shape == (10, 2)

        dm = DataFrame(columns=["A", "B"])
        assert dm.values.shape == (0, 2)

        dm = DataFrame(index=np.arange(10))
        assert dm.values.shape == (10, 0)

        # can't cast
        mat = np.array(["foo", "bar"], dtype=object).reshape(2, 1)
        with pytest.raises(ValueError, match="cast"):
            DataFrame(mat, index=[0, 1], columns=[0], dtype=float)

        dm = DataFrame(DataFrame(float_frame._series))
        tm.assert_frame_equal(dm, float_frame)

        # int cast
        dm = DataFrame(
            {"A": np.ones(10, dtype=int), "B": np.ones(10, dtype=np.float64)},
            index=np.arange(10),
        )

        assert len(dm.columns) == 2
        assert dm.values.dtype == np.float64

    def test_constructor_empty_list(self):
        df = DataFrame([], index=[])
        expected = DataFrame(index=[])
        tm.assert_frame_equal(df, expected)

        # GH 9939
        df = DataFrame([], columns=["A", "B"])
        expected = DataFrame({}, columns=["A", "B"])
        tm.assert_frame_equal(df, expected)

        # Empty generator: list(empty_gen()) == []
        def empty_gen():
            return
            yield

        df = DataFrame(empty_gen(), columns=["A", "B"])
        tm.assert_frame_equal(df, expected)

    def test_constructor_list_of_lists(self):
        # GH #484
        df = DataFrame(data=[[1, "a"], [2, "b"]], columns=["num", "str"])
        assert is_integer_dtype(df["num"])
        assert df["str"].dtype == np.object_

        # GH 4851
        # list of 0-dim ndarrays
        expected = DataFrame({0: np.arange(10)})
        data = [np.array(x) for x in range(10)]
        result = DataFrame(data)
        tm.assert_frame_equal(result, expected)

    def test_constructor_sequence_like(self):
        # GH 3783
        # collections.Squence like

        class DummyContainer(abc.Sequence):
            def __init__(self, lst):
                self._lst = lst

            def __getitem__(self, n):
                return self._lst.__getitem__(n)

            def __len__(self, n):
                return self._lst.__len__()

        lst_containers = [DummyContainer([1, "a"]), DummyContainer([2, "b"])]
        columns = ["num", "str"]
        result = DataFrame(lst_containers, columns=columns)
        expected = DataFrame([[1, "a"], [2, "b"]], columns=columns)
        tm.assert_frame_equal(result, expected, check_dtype=False)

        # GH 4297
        # support Array
        import array

        result = DataFrame({"A": array.array("i", range(10))})
        expected = DataFrame({"A": list(range(10))})
        tm.assert_frame_equal(result, expected, check_dtype=False)

        expected = DataFrame([list(range(10)), list(range(10))])
        result = DataFrame([array.array("i", range(10)), array.array("i", range(10))])
        tm.assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_range(self):
        # GH26342
        result = DataFrame(range(10))
        expected = DataFrame(list(range(10)))
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_ranges(self):
        result = DataFrame([range(10), range(10)])
        expected = DataFrame([list(range(10)), list(range(10))])
        tm.assert_frame_equal(result, expected)

    def test_constructor_iterable(self):
        # GH 21987
        class Iter:
            def __iter__(self):
                for i in range(10):
                    yield [1, 2, 3]

        expected = DataFrame([[1, 2, 3]] * 10)
        result = DataFrame(Iter())
        tm.assert_frame_equal(result, expected)

    def test_constructor_iterator(self):
        result = DataFrame(iter(range(10)))
        expected = DataFrame(list(range(10)))
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_iterators(self):
        result = DataFrame([iter(range(10)), iter(range(10))])
        expected = DataFrame([list(range(10)), list(range(10))])
        tm.assert_frame_equal(result, expected)

    def test_constructor_generator(self):
        # related #2305

        gen1 = (i for i in range(10))
        gen2 = (i for i in range(10))

        expected = DataFrame([list(range(10)), list(range(10))])
        result = DataFrame([gen1, gen2])
        tm.assert_frame_equal(result, expected)

        gen = ([i, "a"] for i in range(10))
        result = DataFrame(gen)
        expected = DataFrame({0: range(10), 1: "a"})
        tm.assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_list_of_odicts(self):
        data = [
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4], ["d", 6]]),
            OrderedDict([["a", 1.5], ["b", 3], ["d", 6]]),
            OrderedDict([["a", 1.5], ["d", 6]]),
            OrderedDict(),
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4]]),
            OrderedDict([["b", 3], ["c", 4], ["d", 6]]),
        ]

        result = DataFrame(data)
        expected = DataFrame.from_dict(
            dict(zip(range(len(data)), data)), orient="index"
        )
        tm.assert_frame_equal(result, expected.reindex(result.index))

        result = DataFrame([{}])
        expected = DataFrame(index=[0])
        tm.assert_frame_equal(result, expected)

    def test_constructor_ordered_dict_preserve_order(self):
        # see gh-13304
        expected = DataFrame([[2, 1]], columns=["b", "a"])

        data = OrderedDict()
        data["b"] = [2]
        data["a"] = [1]

        result = DataFrame(data)
        tm.assert_frame_equal(result, expected)

        data = OrderedDict()
        data["b"] = 2
        data["a"] = 1

        result = DataFrame([data])
        tm.assert_frame_equal(result, expected)

    def test_constructor_ordered_dict_conflicting_orders(self):
        # the first dict element sets the ordering for the DataFrame,
        # even if there are conflicting orders from subsequent ones
        row_one = OrderedDict()
        row_one["b"] = 2
        row_one["a"] = 1

        row_two = OrderedDict()
        row_two["a"] = 1
        row_two["b"] = 2

        row_three = {"b": 2, "a": 1}

        expected = DataFrame([[2, 1], [2, 1]], columns=["b", "a"])
        result = DataFrame([row_one, row_two])
        tm.assert_frame_equal(result, expected)

        expected = DataFrame([[2, 1], [2, 1], [2, 1]], columns=["b", "a"])
        result = DataFrame([row_one, row_two, row_three])
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_series(self):
        data = [
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 4.0]]),
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 6.0]]),
        ]
        sdict = OrderedDict(zip(["x", "y"], data))
        idx = Index(["a", "b", "c"])

        # all named
        data2 = [
            Series([1.5, 3, 4], idx, dtype="O", name="x"),
            Series([1.5, 3, 6], idx, name="y"),
        ]
        result = DataFrame(data2)
        expected = DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected)

        # some unnamed
        data2 = [
            Series([1.5, 3, 4], idx, dtype="O", name="x"),
            Series([1.5, 3, 6], idx),
        ]
        result = DataFrame(data2)

        sdict = OrderedDict(zip(["x", "Unnamed 0"], data))
        expected = DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected)

        # none named
        data = [
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4], ["d", 6]]),
            OrderedDict([["a", 1.5], ["b", 3], ["d", 6]]),
            OrderedDict([["a", 1.5], ["d", 6]]),
            OrderedDict(),
            OrderedDict([["a", 1.5], ["b", 3], ["c", 4]]),
            OrderedDict([["b", 3], ["c", 4], ["d", 6]]),
        ]
        data = [
            create_series_with_explicit_dtype(d, dtype_if_empty=object) for d in data
        ]

        result = DataFrame(data)
        sdict = OrderedDict(zip(range(len(data)), data))
        expected = DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected.reindex(result.index))

        result2 = DataFrame(data, index=np.arange(6))
        tm.assert_frame_equal(result, result2)

        result = DataFrame([Series(dtype=object)])
        expected = DataFrame(index=[0])
        tm.assert_frame_equal(result, expected)

        data = [
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 4.0]]),
            OrderedDict([["a", 1.5], ["b", 3.0], ["c", 6.0]]),
        ]
        sdict = OrderedDict(zip(range(len(data)), data))

        idx = Index(["a", "b", "c"])
        data2 = [Series([1.5, 3, 4], idx, dtype="O"), Series([1.5, 3, 6], idx)]
        result = DataFrame(data2)
        expected = DataFrame.from_dict(sdict, orient="index")
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_series_aligned_index(self):
        series = [pd.Series(i, index=["b", "a", "c"], name=str(i)) for i in range(3)]
        result = pd.DataFrame(series)
        expected = pd.DataFrame(
            {"b": [0, 1, 2], "a": [0, 1, 2], "c": [0, 1, 2]},
            columns=["b", "a", "c"],
            index=["0", "1", "2"],
        )
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_derived_dicts(self):
        class CustomDict(dict):
            pass

        d = {"a": 1.5, "b": 3}

        data_custom = [CustomDict(d)]
        data = [d]

        result_custom = DataFrame(data_custom)
        result = DataFrame(data)
        tm.assert_frame_equal(result, result_custom)

    def test_constructor_ragged(self):
        data = {"A": np.random.randn(10), "B": np.random.randn(8)}
        with pytest.raises(ValueError, match="arrays must all be same length"):
            DataFrame(data)

    def test_constructor_scalar(self):
        idx = Index(range(3))
        df = DataFrame({"a": 0}, index=idx)
        expected = DataFrame({"a": [0, 0, 0]}, index=idx)
        tm.assert_frame_equal(df, expected, check_dtype=False)

    def test_constructor_Series_copy_bug(self, float_frame):
        df = DataFrame(float_frame["A"], index=float_frame.index, columns=["A"])
        df.copy()

    def test_constructor_mixed_dict_and_Series(self):
        data = {}
        data["A"] = {"foo": 1, "bar": 2, "baz": 3}
        data["B"] = Series([4, 3, 2, 1], index=["bar", "qux", "baz", "foo"])

        result = DataFrame(data)
        assert result.index.is_monotonic

        # ordering ambiguous, raise exception
        with pytest.raises(ValueError, match="ambiguous ordering"):
            DataFrame({"A": ["a", "b"], "B": {"a": "a", "b": "b"}})

        # this is OK though
        result = DataFrame({"A": ["a", "b"], "B": Series(["a", "b"], index=["a", "b"])})
        expected = DataFrame({"A": ["a", "b"], "B": ["a", "b"]}, index=["a", "b"])
        tm.assert_frame_equal(result, expected)

    def test_constructor_mixed_type_rows(self):
        # Issue 25075
        data = [[1, 2], (3, 4)]
        result = DataFrame(data)
        expected = DataFrame([[1, 2], [3, 4]])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "tuples,lists",
        [
            ((), []),
            ((()), []),
            (((), ()), [(), ()]),
            (((), ()), [[], []]),
            (([], []), [[], []]),
            (([1, 2, 3], [4, 5, 6]), [[1, 2, 3], [4, 5, 6]]),
        ],
    )
    def test_constructor_tuple(self, tuples, lists):
        # GH 25691
        result = DataFrame(tuples)
        expected = DataFrame(lists)
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_tuples(self):
        result = DataFrame({"A": [(1, 2), (3, 4)]})
        expected = DataFrame({"A": Series([(1, 2), (3, 4)])})
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_namedtuples(self):
        # GH11181
        from collections import namedtuple

        named_tuple = namedtuple("Pandas", list("ab"))
        tuples = [named_tuple(1, 3), named_tuple(2, 4)]
        expected = DataFrame({"a": [1, 2], "b": [3, 4]})
        result = DataFrame(tuples)
        tm.assert_frame_equal(result, expected)

        # with columns
        expected = DataFrame({"y": [1, 2], "z": [3, 4]})
        result = DataFrame(tuples, columns=["y", "z"])
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_of_dict_order(self):
        # GH10056
        data = [
            {"First": 1, "Second": 4, "Third": 7, "Fourth": 10},
            {"Second": 5, "First": 2, "Fourth": 11, "Third": 8},
            {"Second": 6, "First": 3, "Fourth": 12, "Third": 9, "YYY": 14, "XXX": 13},
        ]
        expected = DataFrame(
            {
                "First": [1, 2, 3],
                "Second": [4, 5, 6],
                "Third": [7, 8, 9],
                "Fourth": [10, 11, 12],
                "YYY": [None, None, 14],
                "XXX": [None, None, 13],
            }
        )
        result = DataFrame(data)
        tm.assert_frame_equal(result, expected)

    def test_constructor_orient(self, float_string_frame):
        data_dict = float_string_frame.T._series
        recons = DataFrame.from_dict(data_dict, orient="index")
        expected = float_string_frame.reindex(index=recons.index)
        tm.assert_frame_equal(recons, expected)

        # dict of sequence
        a = {"hi": [32, 3, 3], "there": [3, 5, 3]}
        rs = DataFrame.from_dict(a, orient="index")
        xp = DataFrame.from_dict(a).T.reindex(list(a.keys()))
        tm.assert_frame_equal(rs, xp)

    def test_constructor_from_ordered_dict(self):
        # GH8425
        a = OrderedDict(
            [
                ("one", OrderedDict([("col_a", "foo1"), ("col_b", "bar1")])),
                ("two", OrderedDict([("col_a", "foo2"), ("col_b", "bar2")])),
                ("three", OrderedDict([("col_a", "foo3"), ("col_b", "bar3")])),
            ]
        )
        expected = DataFrame.from_dict(a, orient="columns").T
        result = DataFrame.from_dict(a, orient="index")
        tm.assert_frame_equal(result, expected)

    def test_from_dict_columns_parameter(self):
        # GH 18529
        # Test new columns parameter for from_dict that was added to make
        # from_items(..., orient='index', columns=[...]) easier to replicate
        result = DataFrame.from_dict(
            OrderedDict([("A", [1, 2]), ("B", [4, 5])]),
            orient="index",
            columns=["one", "two"],
        )
        expected = DataFrame([[1, 2], [4, 5]], index=["A", "B"], columns=["one", "two"])
        tm.assert_frame_equal(result, expected)

        msg = "cannot use columns parameter with orient='columns'"
        with pytest.raises(ValueError, match=msg):
            DataFrame.from_dict(
                dict([("A", [1, 2]), ("B", [4, 5])]),
                orient="columns",
                columns=["one", "two"],
            )
        with pytest.raises(ValueError, match=msg):
            DataFrame.from_dict(
                dict([("A", [1, 2]), ("B", [4, 5])]), columns=["one", "two"]
            )

    @pytest.mark.parametrize(
        "data_dict, keys",
        [
            ([{("a",): 1}, {("a",): 2}], [("a",)]),
            ([OrderedDict([(("a",), 1), (("b",), 2)])], [("a",), ("b",)]),
            ([{("a", "b"): 1}], [("a", "b")]),
        ],
    )
    def test_constructor_from_dict_tuples(self, data_dict, keys):
        # GH 16769
        df = DataFrame.from_dict(data_dict)

        result = df.columns
        expected = Index(keys, dtype="object", tupleize_cols=False)

        tm.assert_index_equal(result, expected)

    def test_constructor_Series_named(self):
        a = Series([1, 2, 3], index=["a", "b", "c"], name="x")
        df = DataFrame(a)
        assert df.columns[0] == "x"
        tm.assert_index_equal(df.index, a.index)

        # ndarray like
        arr = np.random.randn(10)
        s = Series(arr, name="x")
        df = DataFrame(s)
        expected = DataFrame(dict(x=s))
        tm.assert_frame_equal(df, expected)

        s = Series(arr, index=range(3, 13))
        df = DataFrame(s)
        expected = DataFrame({0: s})
        tm.assert_frame_equal(df, expected)

        msg = r"Shape of passed values is \(10, 1\), indices imply \(10, 2\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame(s, columns=[1, 2])

        # #2234
        a = Series([], name="x", dtype=object)
        df = DataFrame(a)
        assert df.columns[0] == "x"

        # series with name and w/o
        s1 = Series(arr, name="x")
        df = DataFrame([s1, arr]).T
        expected = DataFrame({"x": s1, "Unnamed 0": arr}, columns=["x", "Unnamed 0"])
        tm.assert_frame_equal(df, expected)

        # this is a bit non-intuitive here; the series collapse down to arrays
        df = DataFrame([arr, s1]).T
        expected = DataFrame({1: s1, 0: arr}, columns=[0, 1])
        tm.assert_frame_equal(df, expected)

    def test_constructor_Series_named_and_columns(self):
        # GH 9232 validation

        s0 = Series(range(5), name=0)
        s1 = Series(range(5), name=1)

        # matching name and column gives standard frame
        tm.assert_frame_equal(pd.DataFrame(s0, columns=[0]), s0.to_frame())
        tm.assert_frame_equal(pd.DataFrame(s1, columns=[1]), s1.to_frame())

        # non-matching produces empty frame
        assert pd.DataFrame(s0, columns=[1]).empty
        assert pd.DataFrame(s1, columns=[0]).empty

    def test_constructor_Series_differently_indexed(self):
        # name
        s1 = Series([1, 2, 3], index=["a", "b", "c"], name="x")

        # no name
        s2 = Series([1, 2, 3], index=["a", "b", "c"])

        other_index = Index(["a", "b"])

        df1 = DataFrame(s1, index=other_index)
        exp1 = DataFrame(s1.reindex(other_index))
        assert df1.columns[0] == "x"
        tm.assert_frame_equal(df1, exp1)

        df2 = DataFrame(s2, index=other_index)
        exp2 = DataFrame(s2.reindex(other_index))
        assert df2.columns[0] == 0
        tm.assert_index_equal(df2.index, other_index)
        tm.assert_frame_equal(df2, exp2)

    def test_constructor_manager_resize(self, float_frame):
        index = list(float_frame.index[:5])
        columns = list(float_frame.columns[:3])

        result = DataFrame(float_frame._data, index=index, columns=columns)
        tm.assert_index_equal(result.index, Index(index))
        tm.assert_index_equal(result.columns, Index(columns))

    def test_constructor_mix_series_nonseries(self, float_frame):
        df = DataFrame(
            {"A": float_frame["A"], "B": list(float_frame["B"])}, columns=["A", "B"]
        )
        tm.assert_frame_equal(df, float_frame.loc[:, ["A", "B"]])

        msg = "does not match index length"
        with pytest.raises(ValueError, match=msg):
            DataFrame({"A": float_frame["A"], "B": list(float_frame["B"])[:-2]})

    def test_constructor_miscast_na_int_dtype(self):
        df = DataFrame([[np.nan, 1], [1, 0]], dtype=np.int64)
        expected = DataFrame([[np.nan, 1], [1, 0]])
        tm.assert_frame_equal(df, expected)

    def test_constructor_column_duplicates(self):
        # it works! #2079
        df = DataFrame([[8, 5]], columns=["a", "a"])
        edf = DataFrame([[8, 5]])
        edf.columns = ["a", "a"]

        tm.assert_frame_equal(df, edf)

        idf = DataFrame.from_records([(8, 5)], columns=["a", "a"])

        tm.assert_frame_equal(idf, edf)

        msg = "If using all scalar values, you must pass an index"
        with pytest.raises(ValueError, match=msg):
            DataFrame.from_dict(OrderedDict([("b", 8), ("a", 5), ("a", 6)]))

    def test_constructor_empty_with_string_dtype(self):
        # GH 9428
        expected = DataFrame(index=[0, 1], columns=[0, 1], dtype=object)

        df = DataFrame(index=[0, 1], columns=[0, 1], dtype=str)
        tm.assert_frame_equal(df, expected)
        df = DataFrame(index=[0, 1], columns=[0, 1], dtype=np.str_)
        tm.assert_frame_equal(df, expected)
        df = DataFrame(index=[0, 1], columns=[0, 1], dtype=np.unicode_)
        tm.assert_frame_equal(df, expected)
        df = DataFrame(index=[0, 1], columns=[0, 1], dtype="U5")
        tm.assert_frame_equal(df, expected)

    def test_constructor_single_value(self):
        # expecting single value upcasting here
        df = DataFrame(0.0, index=[1, 2, 3], columns=["a", "b", "c"])
        tm.assert_frame_equal(
            df, DataFrame(np.zeros(df.shape).astype("float64"), df.index, df.columns)
        )

        df = DataFrame(0, index=[1, 2, 3], columns=["a", "b", "c"])
        tm.assert_frame_equal(
            df, DataFrame(np.zeros(df.shape).astype("int64"), df.index, df.columns)
        )

        df = DataFrame("a", index=[1, 2], columns=["a", "c"])
        tm.assert_frame_equal(
            df,
            DataFrame(
                np.array([["a", "a"], ["a", "a"]], dtype=object),
                index=[1, 2],
                columns=["a", "c"],
            ),
        )

        msg = "DataFrame constructor not properly called!"
        with pytest.raises(ValueError, match=msg):
            DataFrame("a", [1, 2])
        with pytest.raises(ValueError, match=msg):
            DataFrame("a", columns=["a", "c"])

        msg = "incompatible data and dtype"
        with pytest.raises(TypeError, match=msg):
            DataFrame("a", [1, 2], ["a", "c"], float)

    def test_constructor_with_datetimes(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype("M8[ns]").name
        objectname = np.dtype(np.object_).name

        # single item
        df = DataFrame(
            {
                "A": 1,
                "B": "foo",
                "C": "bar",
                "D": Timestamp("20010101"),
                "E": datetime(2001, 1, 2, 0, 0),
            },
            index=np.arange(10),
        )
        result = df.dtypes
        expected = Series(
            [np.dtype("int64")]
            + [np.dtype(objectname)] * 2
            + [np.dtype(datetime64name)] * 2,
            index=list("ABCDE"),
        )
        tm.assert_series_equal(result, expected)

        # check with ndarray construction ndim==0 (e.g. we are passing a ndim 0
        # ndarray with a dtype specified)
        df = DataFrame(
            {
                "a": 1.0,
                "b": 2,
                "c": "foo",
                floatname: np.array(1.0, dtype=floatname),
                intname: np.array(1, dtype=intname),
            },
            index=np.arange(10),
        )
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")]
            + [np.dtype("int64")]
            + [np.dtype("object")]
            + [np.dtype("float64")]
            + [np.dtype(intname)],
            index=["a", "b", "c", floatname, intname],
        )
        tm.assert_series_equal(result, expected)

        # check with ndarray construction ndim>0
        df = DataFrame(
            {
                "a": 1.0,
                "b": 2,
                "c": "foo",
                floatname: np.array([1.0] * 10, dtype=floatname),
                intname: np.array([1] * 10, dtype=intname),
            },
            index=np.arange(10),
        )
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")]
            + [np.dtype("int64")]
            + [np.dtype("object")]
            + [np.dtype("float64")]
            + [np.dtype(intname)],
            index=["a", "b", "c", floatname, intname],
        )
        tm.assert_series_equal(result, expected)

        # GH 2809
        ind = date_range(start="2000-01-01", freq="D", periods=10)
        datetimes = [ts.to_pydatetime() for ts in ind]
        datetime_s = Series(datetimes)
        assert datetime_s.dtype == "M8[ns]"

        # GH 2810
        ind = date_range(start="2000-01-01", freq="D", periods=10)
        datetimes = [ts.to_pydatetime() for ts in ind]
        dates = [ts.date() for ts in ind]
        df = DataFrame(datetimes, columns=["datetimes"])
        df["dates"] = dates
        result = df.dtypes
        expected = Series(
            [np.dtype("datetime64[ns]"), np.dtype("object")],
            index=["datetimes", "dates"],
        )
        tm.assert_series_equal(result, expected)

        # GH 7594
        # don't coerce tz-aware
        import pytz

        tz = pytz.timezone("US/Eastern")
        dt = tz.localize(datetime(2012, 1, 1))

        df = DataFrame({"End Date": dt}, index=[0])
        assert df.iat[0, 0] == dt
        tm.assert_series_equal(
            df.dtypes, Series({"End Date": "datetime64[ns, US/Eastern]"})
        )

        df = DataFrame([{"End Date": dt}])
        assert df.iat[0, 0] == dt
        tm.assert_series_equal(
            df.dtypes, Series({"End Date": "datetime64[ns, US/Eastern]"})
        )

        # tz-aware (UTC and other tz's)
        # GH 8411
        dr = date_range("20130101", periods=3)
        df = DataFrame({"value": dr})
        assert df.iat[0, 0].tz is None
        dr = date_range("20130101", periods=3, tz="UTC")
        df = DataFrame({"value": dr})
        assert str(df.iat[0, 0].tz) == "UTC"
        dr = date_range("20130101", periods=3, tz="US/Eastern")
        df = DataFrame({"value": dr})
        assert str(df.iat[0, 0].tz) == "US/Eastern"

        # GH 7822
        # preserver an index with a tz on dict construction
        i = date_range("1/1/2011", periods=5, freq="10s", tz="US/Eastern")

        expected = DataFrame({"a": i.to_series().reset_index(drop=True)})
        df = DataFrame()
        df["a"] = i
        tm.assert_frame_equal(df, expected)

        df = DataFrame({"a": i})
        tm.assert_frame_equal(df, expected)

        # multiples
        i_no_tz = date_range("1/1/2011", periods=5, freq="10s")
        df = DataFrame({"a": i, "b": i_no_tz})
        expected = DataFrame({"a": i.to_series().reset_index(drop=True), "b": i_no_tz})
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "arr",
        [
            np.array([None, None, None, None, datetime.now(), None]),
            np.array([None, None, datetime.now(), None]),
            [[np.datetime64("NaT")], [None]],
            [[np.datetime64("NaT")], [pd.NaT]],
            [[None], [np.datetime64("NaT")]],
            [[None], [pd.NaT]],
            [[pd.NaT], [np.datetime64("NaT")]],
            [[pd.NaT], [None]],
        ],
    )
    def test_constructor_datetimes_with_nulls(self, arr):
        # gh-15869, GH#11220
        result = DataFrame(arr).dtypes
        expected = Series([np.dtype("datetime64[ns]")])
        tm.assert_series_equal(result, expected)

    def test_constructor_for_list_with_dtypes(self):
        # test list of lists/ndarrays
        df = DataFrame([np.arange(5) for x in range(5)])
        result = df.dtypes
        expected = Series([np.dtype("int64")] * 5)
        tm.assert_series_equal(result, expected)

        df = DataFrame([np.array(np.arange(5), dtype="int32") for x in range(5)])
        result = df.dtypes
        expected = Series([np.dtype("int64")] * 5)
        tm.assert_series_equal(result, expected)

        # overflow issue? (we always expecte int64 upcasting here)
        df = DataFrame({"a": [2 ** 31, 2 ** 31 + 1]})
        assert df.dtypes.iloc[0] == np.dtype("int64")

        # GH #2751 (construction with no index specified), make sure we cast to
        # platform values
        df = DataFrame([1, 2])
        assert df.dtypes.iloc[0] == np.dtype("int64")

        df = DataFrame([1.0, 2.0])
        assert df.dtypes.iloc[0] == np.dtype("float64")

        df = DataFrame({"a": [1, 2]})
        assert df.dtypes.iloc[0] == np.dtype("int64")

        df = DataFrame({"a": [1.0, 2.0]})
        assert df.dtypes.iloc[0] == np.dtype("float64")

        df = DataFrame({"a": 1}, index=range(3))
        assert df.dtypes.iloc[0] == np.dtype("int64")

        df = DataFrame({"a": 1.0}, index=range(3))
        assert df.dtypes.iloc[0] == np.dtype("float64")

        # with object list
        df = DataFrame(
            {
                "a": [1, 2, 4, 7],
                "b": [1.2, 2.3, 5.1, 6.3],
                "c": list("abcd"),
                "d": [datetime(2000, 1, 1) for i in range(4)],
                "e": [1.0, 2, 4.0, 7],
            }
        )
        result = df.dtypes
        expected = Series(
            [
                np.dtype("int64"),
                np.dtype("float64"),
                np.dtype("object"),
                np.dtype("datetime64[ns]"),
                np.dtype("float64"),
            ],
            index=list("abcde"),
        )
        tm.assert_series_equal(result, expected)

    def test_constructor_frame_copy(self, float_frame):
        cop = DataFrame(float_frame, copy=True)
        cop["A"] = 5
        assert (cop["A"] == 5).all()
        assert not (float_frame["A"] == 5).all()

    def test_constructor_ndarray_copy(self, float_frame):
        df = DataFrame(float_frame.values)

        float_frame.values[5] = 5
        assert (df.values[5] == 5).all()

        df = DataFrame(float_frame.values, copy=True)
        float_frame.values[6] = 6
        assert not (df.values[6] == 6).all()

    def test_constructor_series_copy(self, float_frame):
        series = float_frame._series

        df = DataFrame({"A": series["A"]})
        df["A"][:] = 5

        assert not (series["A"] == 5).all()

    def test_constructor_with_nas(self):
        # GH 5016
        # na's in indices

        def check(df):
            for i in range(len(df.columns)):
                df.iloc[:, i]

            indexer = np.arange(len(df.columns))[isna(df.columns)]

            # No NaN found -> error
            if len(indexer) == 0:
                msg = (
                    "cannot do label indexing on "
                    r"<class 'pandas\.core\.indexes\.range\.RangeIndex'> "
                    r"with these indexers \[nan\] of <class 'float'>"
                )
                with pytest.raises(TypeError, match=msg):
                    df.loc[:, np.nan]
            # single nan should result in Series
            elif len(indexer) == 1:
                tm.assert_series_equal(df.iloc[:, indexer[0]], df.loc[:, np.nan])
            # multiple nans should result in DataFrame
            else:
                tm.assert_frame_equal(df.iloc[:, indexer], df.loc[:, np.nan])

        df = DataFrame([[1, 2, 3], [4, 5, 6]], index=[1, np.nan])
        check(df)

        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=[1.1, 2.2, np.nan])
        check(df)

        df = DataFrame([[0, 1, 2, 3], [4, 5, 6, 7]], columns=[np.nan, 1.1, 2.2, np.nan])
        check(df)

        df = DataFrame(
            [[0.0, 1, 2, 3.0], [4, 5, 6, 7]], columns=[np.nan, 1.1, 2.2, np.nan]
        )
        check(df)

        # GH 21428 (non-unique columns)
        df = DataFrame([[0.0, 1, 2, 3.0], [4, 5, 6, 7]], columns=[np.nan, 1, 2, 2])
        check(df)

    def test_constructor_lists_to_object_dtype(self):
        # from #1074
        d = DataFrame({"a": [np.nan, False]})
        assert d["a"].dtype == np.object_
        assert not d["a"][1]

    def test_constructor_categorical(self):

        # GH8626

        # dict creation
        df = DataFrame({"A": list("abc")}, dtype="category")
        expected = Series(list("abc"), dtype="category", name="A")
        tm.assert_series_equal(df["A"], expected)

        # to_frame
        s = Series(list("abc"), dtype="category")
        result = s.to_frame()
        expected = Series(list("abc"), dtype="category", name=0)
        tm.assert_series_equal(result[0], expected)
        result = s.to_frame(name="foo")
        expected = Series(list("abc"), dtype="category", name="foo")
        tm.assert_series_equal(result["foo"], expected)

        # list-like creation
        df = DataFrame(list("abc"), dtype="category")
        expected = Series(list("abc"), dtype="category", name=0)
        tm.assert_series_equal(df[0], expected)

        # ndim != 1
        df = DataFrame([Categorical(list("abc"))])
        expected = DataFrame({0: Series(list("abc"), dtype="category")})
        tm.assert_frame_equal(df, expected)

        df = DataFrame([Categorical(list("abc")), Categorical(list("abd"))])
        expected = DataFrame(
            {
                0: Series(list("abc"), dtype="category"),
                1: Series(list("abd"), dtype="category"),
            },
            columns=[0, 1],
        )
        tm.assert_frame_equal(df, expected)

        # mixed
        df = DataFrame([Categorical(list("abc")), list("def")])
        expected = DataFrame(
            {0: Series(list("abc"), dtype="category"), 1: list("def")}, columns=[0, 1]
        )
        tm.assert_frame_equal(df, expected)

        # invalid (shape)
        msg = r"Shape of passed values is \(6, 2\), indices imply \(3, 2\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame([Categorical(list("abc")), Categorical(list("abdefg"))])

        # ndim > 1
        msg = "> 1 ndim Categorical are not supported at this time"
        with pytest.raises(NotImplementedError, match=msg):
            Categorical(np.array([list("abcd")]))

    def test_constructor_categorical_series(self):

        items = [1, 2, 3, 1]
        exp = Series(items).astype("category")
        res = Series(items, dtype="category")
        tm.assert_series_equal(res, exp)

        items = ["a", "b", "c", "a"]
        exp = Series(items).astype("category")
        res = Series(items, dtype="category")
        tm.assert_series_equal(res, exp)

        # insert into frame with different index
        # GH 8076
        index = date_range("20000101", periods=3)
        expected = Series(
            Categorical(values=[np.nan, np.nan, np.nan], categories=["a", "b", "c"])
        )
        expected.index = index

        expected = DataFrame({"x": expected})
        df = DataFrame({"x": Series(["a", "b", "c"], dtype="category")}, index=index)
        tm.assert_frame_equal(df, expected)

    def test_from_records_to_records(self):
        # from numpy documentation
        arr = np.zeros((2,), dtype=("i4,f4,a10"))
        arr[:] = [(1, 2.0, "Hello"), (2, 3.0, "World")]

        # TODO(wesm): unused
        frame = DataFrame.from_records(arr)  # noqa

        index = pd.Index(np.arange(len(arr))[::-1])
        indexed_frame = DataFrame.from_records(arr, index=index)
        tm.assert_index_equal(indexed_frame.index, index)

        # without names, it should go to last ditch
        arr2 = np.zeros((2, 3))
        tm.assert_frame_equal(DataFrame.from_records(arr2), DataFrame(arr2))

        # wrong length
        msg = r"Shape of passed values is \(2, 3\), indices imply \(1, 3\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame.from_records(arr, index=index[:-1])

        indexed_frame = DataFrame.from_records(arr, index="f1")

        # what to do?
        records = indexed_frame.to_records()
        assert len(records.dtype.names) == 3

        records = indexed_frame.to_records(index=False)
        assert len(records.dtype.names) == 2
        assert "index" not in records.dtype.names

    def test_from_records_nones(self):
        tuples = [(1, 2, None, 3), (1, 2, None, 3), (None, 2, 5, 3)]

        df = DataFrame.from_records(tuples, columns=["a", "b", "c", "d"])
        assert np.isnan(df["c"][0])

    def test_from_records_iterator(self):
        arr = np.array(
            [(1.0, 1.0, 2, 2), (3.0, 3.0, 4, 4), (5.0, 5.0, 6, 6), (7.0, 7.0, 8, 8)],
            dtype=[
                ("x", np.float64),
                ("u", np.float32),
                ("y", np.int64),
                ("z", np.int32),
            ],
        )
        df = DataFrame.from_records(iter(arr), nrows=2)
        xp = DataFrame(
            {
                "x": np.array([1.0, 3.0], dtype=np.float64),
                "u": np.array([1.0, 3.0], dtype=np.float32),
                "y": np.array([2, 4], dtype=np.int64),
                "z": np.array([2, 4], dtype=np.int32),
            }
        )
        tm.assert_frame_equal(df.reindex_like(xp), xp)

        # no dtypes specified here, so just compare with the default
        arr = [(1.0, 2), (3.0, 4), (5.0, 6), (7.0, 8)]
        df = DataFrame.from_records(iter(arr), columns=["x", "y"], nrows=2)
        tm.assert_frame_equal(df, xp.reindex(columns=["x", "y"]), check_dtype=False)

    def test_from_records_tuples_generator(self):
        def tuple_generator(length):
            for i in range(length):
                letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                yield (i, letters[i % len(letters)], i / length)

        columns_names = ["Integer", "String", "Float"]
        columns = [
            [i[j] for i in tuple_generator(10)] for j in range(len(columns_names))
        ]
        data = {"Integer": columns[0], "String": columns[1], "Float": columns[2]}
        expected = DataFrame(data, columns=columns_names)

        generator = tuple_generator(10)
        result = DataFrame.from_records(generator, columns=columns_names)
        tm.assert_frame_equal(result, expected)

    def test_from_records_lists_generator(self):
        def list_generator(length):
            for i in range(length):
                letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                yield [i, letters[i % len(letters)], i / length]

        columns_names = ["Integer", "String", "Float"]
        columns = [
            [i[j] for i in list_generator(10)] for j in range(len(columns_names))
        ]
        data = {"Integer": columns[0], "String": columns[1], "Float": columns[2]}
        expected = DataFrame(data, columns=columns_names)

        generator = list_generator(10)
        result = DataFrame.from_records(generator, columns=columns_names)
        tm.assert_frame_equal(result, expected)

    def test_from_records_columns_not_modified(self):
        tuples = [(1, 2, 3), (1, 2, 3), (2, 5, 3)]

        columns = ["a", "b", "c"]
        original_columns = list(columns)

        df = DataFrame.from_records(tuples, columns=columns, index="a")  # noqa

        assert columns == original_columns

    def test_from_records_decimal(self):
        from decimal import Decimal

        tuples = [(Decimal("1.5"),), (Decimal("2.5"),), (None,)]

        df = DataFrame.from_records(tuples, columns=["a"])
        assert df["a"].dtype == object

        df = DataFrame.from_records(tuples, columns=["a"], coerce_float=True)
        assert df["a"].dtype == np.float64
        assert np.isnan(df["a"].values[-1])

    def test_from_records_duplicates(self):
        result = DataFrame.from_records([(1, 2, 3), (4, 5, 6)], columns=["a", "b", "a"])

        expected = DataFrame([(1, 2, 3), (4, 5, 6)], columns=["a", "b", "a"])

        tm.assert_frame_equal(result, expected)

    def test_from_records_set_index_name(self):
        def create_dict(order_id):
            return {
                "order_id": order_id,
                "quantity": np.random.randint(1, 10),
                "price": np.random.randint(1, 10),
            }

        documents = [create_dict(i) for i in range(10)]
        # demo missing data
        documents.append({"order_id": 10, "quantity": 5})

        result = DataFrame.from_records(documents, index="order_id")
        assert result.index.name == "order_id"

        # MultiIndex
        result = DataFrame.from_records(documents, index=["order_id", "quantity"])
        assert result.index.names == ("order_id", "quantity")

    def test_from_records_misc_brokenness(self):
        # #2179

        data = {1: ["foo"], 2: ["bar"]}

        result = DataFrame.from_records(data, columns=["a", "b"])
        exp = DataFrame(data, columns=["a", "b"])
        tm.assert_frame_equal(result, exp)

        # overlap in index/index_names

        data = {"a": [1, 2, 3], "b": [4, 5, 6]}

        result = DataFrame.from_records(data, index=["a", "b", "c"])
        exp = DataFrame(data, index=["a", "b", "c"])
        tm.assert_frame_equal(result, exp)

        # GH 2623
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), "hi"])  # test col upconverts to obj
        df2_obj = DataFrame.from_records(rows, columns=["date", "test"])
        result = df2_obj.dtypes
        expected = Series(
            [np.dtype("datetime64[ns]"), np.dtype("object")], index=["date", "test"]
        )
        tm.assert_series_equal(result, expected)

        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 1])
        df2_obj = DataFrame.from_records(rows, columns=["date", "test"])
        result = df2_obj.dtypes
        expected = Series(
            [np.dtype("datetime64[ns]"), np.dtype("int64")], index=["date", "test"]
        )
        tm.assert_series_equal(result, expected)

    def test_from_records_empty(self):
        # 3562
        result = DataFrame.from_records([], columns=["a", "b", "c"])
        expected = DataFrame(columns=["a", "b", "c"])
        tm.assert_frame_equal(result, expected)

        result = DataFrame.from_records([], columns=["a", "b", "b"])
        expected = DataFrame(columns=["a", "b", "b"])
        tm.assert_frame_equal(result, expected)

    def test_from_records_empty_with_nonempty_fields_gh3682(self):
        a = np.array([(1, 2)], dtype=[("id", np.int64), ("value", np.int64)])
        df = DataFrame.from_records(a, index="id")
        tm.assert_index_equal(df.index, Index([1], name="id"))
        assert df.index.name == "id"
        tm.assert_index_equal(df.columns, Index(["value"]))

        b = np.array([], dtype=[("id", np.int64), ("value", np.int64)])
        df = DataFrame.from_records(b, index="id")
        tm.assert_index_equal(df.index, Index([], name="id"))
        assert df.index.name == "id"

    def test_from_records_with_datetimes(self):

        # this may fail on certain platforms because of a numpy issue
        # related GH6140
        if not is_platform_little_endian():
            pytest.skip("known failure of test on non-little endian")

        # construction with a null in a recarray
        # GH 6140
        expected = DataFrame({"EXPIRY": [datetime(2005, 3, 1, 0, 0), None]})

        arrdata = [np.array([datetime(2005, 3, 1, 0, 0), None])]
        dtypes = [("EXPIRY", "<M8[ns]")]

        try:
            recarray = np.core.records.fromarrays(arrdata, dtype=dtypes)
        except (ValueError):
            pytest.skip("known failure of numpy rec array creation")

        result = DataFrame.from_records(recarray)
        tm.assert_frame_equal(result, expected)

        # coercion should work too
        arrdata = [np.array([datetime(2005, 3, 1, 0, 0), None])]
        dtypes = [("EXPIRY", "<M8[m]")]
        recarray = np.core.records.fromarrays(arrdata, dtype=dtypes)
        result = DataFrame.from_records(recarray)
        tm.assert_frame_equal(result, expected)

    def test_from_records_sequencelike(self):
        df = DataFrame(
            {
                "A": np.array(np.random.randn(6), dtype=np.float64),
                "A1": np.array(np.random.randn(6), dtype=np.float64),
                "B": np.array(np.arange(6), dtype=np.int64),
                "C": ["foo"] * 6,
                "D": np.array([True, False] * 3, dtype=bool),
                "E": np.array(np.random.randn(6), dtype=np.float32),
                "E1": np.array(np.random.randn(6), dtype=np.float32),
                "F": np.array(np.arange(6), dtype=np.int32),
            }
        )

        # this is actually tricky to create the recordlike arrays and
        # have the dtypes be intact
        blocks = df._to_dict_of_blocks()
        tuples = []
        columns = []
        dtypes = []
        for dtype, b in blocks.items():
            columns.extend(b.columns)
            dtypes.extend([(c, np.dtype(dtype).descr[0][1]) for c in b.columns])
        for i in range(len(df.index)):
            tup = []
            for _, b in blocks.items():
                tup.extend(b.iloc[i].values)
            tuples.append(tuple(tup))

        recarray = np.array(tuples, dtype=dtypes).view(np.recarray)
        recarray2 = df.to_records()
        lists = [list(x) for x in tuples]

        # tuples (lose the dtype info)
        result = DataFrame.from_records(tuples, columns=columns).reindex(
            columns=df.columns
        )

        # created recarray and with to_records recarray (have dtype info)
        result2 = DataFrame.from_records(recarray, columns=columns).reindex(
            columns=df.columns
        )
        result3 = DataFrame.from_records(recarray2, columns=columns).reindex(
            columns=df.columns
        )

        # list of tupels (no dtype info)
        result4 = DataFrame.from_records(lists, columns=columns).reindex(
            columns=df.columns
        )

        tm.assert_frame_equal(result, df, check_dtype=False)
        tm.assert_frame_equal(result2, df)
        tm.assert_frame_equal(result3, df)
        tm.assert_frame_equal(result4, df, check_dtype=False)

        # tuples is in the order of the columns
        result = DataFrame.from_records(tuples)
        tm.assert_index_equal(result.columns, pd.RangeIndex(8))

        # test exclude parameter & we are casting the results here (as we don't
        # have dtype info to recover)
        columns_to_test = [columns.index("C"), columns.index("E1")]

        exclude = list(set(range(8)) - set(columns_to_test))
        result = DataFrame.from_records(tuples, exclude=exclude)
        result.columns = [columns[i] for i in sorted(columns_to_test)]
        tm.assert_series_equal(result["C"], df["C"])
        tm.assert_series_equal(result["E1"], df["E1"].astype("float64"))

        # empty case
        result = DataFrame.from_records([], columns=["foo", "bar", "baz"])
        assert len(result) == 0
        tm.assert_index_equal(result.columns, pd.Index(["foo", "bar", "baz"]))

        result = DataFrame.from_records([])
        assert len(result) == 0
        assert len(result.columns) == 0

    def test_from_records_dictlike(self):

        # test the dict methods
        df = DataFrame(
            {
                "A": np.array(np.random.randn(6), dtype=np.float64),
                "A1": np.array(np.random.randn(6), dtype=np.float64),
                "B": np.array(np.arange(6), dtype=np.int64),
                "C": ["foo"] * 6,
                "D": np.array([True, False] * 3, dtype=bool),
                "E": np.array(np.random.randn(6), dtype=np.float32),
                "E1": np.array(np.random.randn(6), dtype=np.float32),
                "F": np.array(np.arange(6), dtype=np.int32),
            }
        )

        # columns is in a different order here than the actual items iterated
        # from the dict
        blocks = df._to_dict_of_blocks()
        columns = []
        for dtype, b in blocks.items():
            columns.extend(b.columns)

        asdict = {x: y for x, y in df.items()}
        asdict2 = {x: y.values for x, y in df.items()}

        # dict of series & dict of ndarrays (have dtype info)
        results = []
        results.append(DataFrame.from_records(asdict).reindex(columns=df.columns))
        results.append(
            DataFrame.from_records(asdict, columns=columns).reindex(columns=df.columns)
        )
        results.append(
            DataFrame.from_records(asdict2, columns=columns).reindex(columns=df.columns)
        )

        for r in results:
            tm.assert_frame_equal(r, df)

    def test_from_records_with_index_data(self):
        df = DataFrame(np.random.randn(10, 3), columns=["A", "B", "C"])

        data = np.random.randn(10)
        df1 = DataFrame.from_records(df, index=data)
        tm.assert_index_equal(df1.index, Index(data))

    def test_from_records_bad_index_column(self):
        df = DataFrame(np.random.randn(10, 3), columns=["A", "B", "C"])

        # should pass
        df1 = DataFrame.from_records(df, index=["C"])
        tm.assert_index_equal(df1.index, Index(df.C))

        df1 = DataFrame.from_records(df, index="C")
        tm.assert_index_equal(df1.index, Index(df.C))

        # should fail
        msg = r"Shape of passed values is \(10, 3\), indices imply \(1, 3\)"
        with pytest.raises(ValueError, match=msg):
            DataFrame.from_records(df, index=[2])
        with pytest.raises(KeyError, match=r"^2$"):
            DataFrame.from_records(df, index=2)

    def test_from_records_non_tuple(self):
        class Record:
            def __init__(self, *args):
                self.args = args

            def __getitem__(self, i):
                return self.args[i]

            def __iter__(self):
                return iter(self.args)

        recs = [Record(1, 2, 3), Record(4, 5, 6), Record(7, 8, 9)]
        tups = [tuple(rec) for rec in recs]

        result = DataFrame.from_records(recs)
        expected = DataFrame.from_records(tups)
        tm.assert_frame_equal(result, expected)

    def test_from_records_len0_with_columns(self):
        # #2633
        result = DataFrame.from_records([], index="foo", columns=["foo", "bar"])
        expected = Index(["bar"])

        assert len(result) == 0
        assert result.index.name == "foo"
        tm.assert_index_equal(result.columns, expected)

    def test_from_records_series_list_dict(self):
        # GH27358
        expected = DataFrame([[{"a": 1, "b": 2}, {"a": 3, "b": 4}]]).T
        data = Series([[{"a": 1, "b": 2}], [{"a": 3, "b": 4}]])
        result = DataFrame.from_records(data)
        tm.assert_frame_equal(result, expected)

    def test_to_frame_with_falsey_names(self):
        # GH 16114
        result = Series(name=0, dtype=object).to_frame().dtypes
        expected = Series({0: object})
        tm.assert_series_equal(result, expected)

        result = DataFrame(Series(name=0, dtype=object)).dtypes
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", [None, "uint8", "category"])
    def test_constructor_range_dtype(self, dtype):
        expected = DataFrame({"A": [0, 1, 2, 3, 4]}, dtype=dtype or "int64")

        # GH 26342
        result = DataFrame(range(5), columns=["A"], dtype=dtype)
        tm.assert_frame_equal(result, expected)

        # GH 16804
        result = DataFrame({"A": range(5)}, dtype=dtype)
        tm.assert_frame_equal(result, expected)

    def test_frame_from_list_subclass(self):
        # GH21226
        class List(list):
            pass

        expected = DataFrame([[1, 2, 3], [4, 5, 6]])
        result = DataFrame(List([List([1, 2, 3]), List([4, 5, 6])]))
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "extension_arr",
        [
            Categorical(list("aabbc")),
            SparseArray([1, np.nan, np.nan, np.nan]),
            IntervalArray([pd.Interval(0, 1), pd.Interval(1, 5)]),
            PeriodArray(pd.period_range(start="1/1/2017", end="1/1/2018", freq="M")),
        ],
    )
    def test_constructor_with_extension_array(self, extension_arr):
        # GH11363
        expected = DataFrame(Series(extension_arr))
        result = DataFrame(extension_arr)
        tm.assert_frame_equal(result, expected)

    def test_datetime_date_tuple_columns_from_dict(self):
        # GH 10863
        v = date.today()
        tup = v, v
        result = DataFrame({tup: Series(range(3), index=range(3))}, columns=[tup])
        expected = DataFrame([0, 1, 2], columns=pd.Index(pd.Series([tup])))
        tm.assert_frame_equal(result, expected)


class TestDataFrameConstructorWithDatetimeTZ:
    def test_from_dict(self):

        # 8260
        # support datetime64 with tz

        idx = Index(date_range("20130101", periods=3, tz="US/Eastern"), name="foo")
        dr = date_range("20130110", periods=3)

        # construction
        df = DataFrame({"A": idx, "B": dr})
        assert df["A"].dtype, "M8[ns, US/Eastern"
        assert df["A"].name == "A"
        tm.assert_series_equal(df["A"], Series(idx, name="A"))
        tm.assert_series_equal(df["B"], Series(dr, name="B"))

    def test_from_index(self):

        # from index
        idx2 = date_range("20130101", periods=3, tz="US/Eastern", name="foo")
        df2 = DataFrame(idx2)
        tm.assert_series_equal(df2["foo"], Series(idx2, name="foo"))
        df2 = DataFrame(Series(idx2))
        tm.assert_series_equal(df2["foo"], Series(idx2, name="foo"))

        idx2 = date_range("20130101", periods=3, tz="US/Eastern")
        df2 = DataFrame(idx2)
        tm.assert_series_equal(df2[0], Series(idx2, name=0))
        df2 = DataFrame(Series(idx2))
        tm.assert_series_equal(df2[0], Series(idx2, name=0))

    def test_frame_dict_constructor_datetime64_1680(self):
        dr = date_range("1/1/2012", periods=10)
        s = Series(dr, index=dr)

        # it works!
        DataFrame({"a": "foo", "b": s}, index=dr)
        DataFrame({"a": "foo", "b": s.values}, index=dr)

    def test_frame_datetime64_mixed_index_ctor_1681(self):
        dr = date_range("2011/1/1", "2012/1/1", freq="W-FRI")
        ts = Series(dr)

        # it works!
        d = DataFrame({"A": "foo", "B": ts}, index=dr)
        assert d["B"].isna().all()

    def test_frame_timeseries_to_records(self):
        index = date_range("1/1/2000", periods=10)
        df = DataFrame(np.random.randn(10, 3), index=index, columns=["a", "b", "c"])

        result = df.to_records()
        result["index"].dtype == "M8[ns]"

        result = df.to_records(index=False)

    def test_frame_timeseries_column(self):
        # GH19157
        dr = date_range(start="20130101T10:00:00", periods=3, freq="T", tz="US/Eastern")
        result = DataFrame(dr, columns=["timestamps"])
        expected = DataFrame(
            {
                "timestamps": [
                    Timestamp("20130101T10:00:00", tz="US/Eastern"),
                    Timestamp("20130101T10:01:00", tz="US/Eastern"),
                    Timestamp("20130101T10:02:00", tz="US/Eastern"),
                ]
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_nested_dict_construction(self):
        # GH22227
        columns = ["Nevada", "Ohio"]
        pop = {
            "Nevada": {2001: 2.4, 2002: 2.9},
            "Ohio": {2000: 1.5, 2001: 1.7, 2002: 3.6},
        }
        result = pd.DataFrame(pop, index=[2001, 2002, 2003], columns=columns)
        expected = pd.DataFrame(
            [(2.4, 1.7), (2.9, 3.6), (np.nan, np.nan)],
            columns=columns,
            index=pd.Index([2001, 2002, 2003]),
        )
        tm.assert_frame_equal(result, expected)

    def test_from_tzaware_object_array(self):
        # GH#26825 2D object array of tzaware timestamps should not raise
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")
        data = dti._data.astype(object).reshape(1, -1)
        df = pd.DataFrame(data)
        assert df.shape == (1, 3)
        assert (df.dtypes == dti.dtype).all()
        assert (df == dti).all().all()

    def test_from_tzaware_mixed_object_array(self):
        # GH#26825
        arr = np.array(
            [
                [
                    Timestamp("2013-01-01 00:00:00"),
                    Timestamp("2013-01-02 00:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
            ],
            dtype=object,
        ).T
        res = DataFrame(arr, columns=["A", "B", "C"])

        expected_dtypes = [
            "datetime64[ns]",
            "datetime64[ns, US/Eastern]",
            "datetime64[ns, CET]",
        ]
        assert (res.dtypes == expected_dtypes).all()

    def test_from_2d_ndarray_with_dtype(self):
        # GH#12513
        array_dim2 = np.arange(10).reshape((5, 2))
        df = pd.DataFrame(array_dim2, dtype="datetime64[ns, UTC]")

        expected = pd.DataFrame(array_dim2).astype("datetime64[ns, UTC]")
        tm.assert_frame_equal(df, expected)
