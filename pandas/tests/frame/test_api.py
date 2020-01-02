from copy import deepcopy
import datetime
import pydoc

import numpy as np
import pytest

from pandas.compat import PY37

import pandas as pd
from pandas import Categorical, DataFrame, Series, compat, date_range, timedelta_range
import pandas.util.testing as tm


class TestDataFrameMisc:
    def test_copy_index_name_checking(self, float_frame):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy
        for attr in ("index", "columns"):
            ind = getattr(float_frame, attr)
            ind.name = None
            cp = float_frame.copy()
            getattr(cp, attr).name = "foo"
            assert getattr(float_frame, attr).name is None

    def test_getitem_pop_assign_name(self, float_frame):
        s = float_frame["A"]
        assert s.name == "A"

        s = float_frame.pop("A")
        assert s.name == "A"

        s = float_frame.loc[:, "B"]
        assert s.name == "B"

        s2 = s.loc[:]
        assert s2.name == "B"

    def test_get_value(self, float_frame):
        for idx in float_frame.index:
            for col in float_frame.columns:
                result = float_frame._get_value(idx, col)
                expected = float_frame[col][idx]
                tm.assert_almost_equal(result, expected)

    def test_add_prefix_suffix(self, float_frame):
        with_prefix = float_frame.add_prefix("foo#")
        expected = pd.Index(["foo#{c}".format(c=c) for c in float_frame.columns])
        tm.assert_index_equal(with_prefix.columns, expected)

        with_suffix = float_frame.add_suffix("#foo")
        expected = pd.Index(["{c}#foo".format(c=c) for c in float_frame.columns])
        tm.assert_index_equal(with_suffix.columns, expected)

        with_pct_prefix = float_frame.add_prefix("%")
        expected = pd.Index(["%{c}".format(c=c) for c in float_frame.columns])
        tm.assert_index_equal(with_pct_prefix.columns, expected)

        with_pct_suffix = float_frame.add_suffix("%")
        expected = pd.Index(["{c}%".format(c=c) for c in float_frame.columns])
        tm.assert_index_equal(with_pct_suffix.columns, expected)

    def test_get_axis(self, float_frame):
        f = float_frame
        assert f._get_axis_number(0) == 0
        assert f._get_axis_number(1) == 1
        assert f._get_axis_number("index") == 0
        assert f._get_axis_number("rows") == 0
        assert f._get_axis_number("columns") == 1

        assert f._get_axis_name(0) == "index"
        assert f._get_axis_name(1) == "columns"
        assert f._get_axis_name("index") == "index"
        assert f._get_axis_name("rows") == "index"
        assert f._get_axis_name("columns") == "columns"

        assert f._get_axis(0) is f.index
        assert f._get_axis(1) is f.columns

        with pytest.raises(ValueError, match="No axis named"):
            f._get_axis_number(2)

        with pytest.raises(ValueError, match="No axis.*foo"):
            f._get_axis_name("foo")

        with pytest.raises(ValueError, match="No axis.*None"):
            f._get_axis_name(None)

        with pytest.raises(ValueError, match="No axis named"):
            f._get_axis_number(None)

    def test_keys(self, float_frame):
        getkeys = float_frame.keys
        assert getkeys() is float_frame.columns

    def test_column_contains_raises(self, float_frame):
        with pytest.raises(TypeError, match="unhashable type: 'Index'"):
            float_frame.columns in float_frame

    def test_tab_completion(self):
        # DataFrame whose columns are identifiers shall have them in __dir__.
        df = pd.DataFrame([list("abcd"), list("efgh")], columns=list("ABCD"))
        for key in list("ABCD"):
            assert key in dir(df)
        assert isinstance(df.__getitem__("A"), pd.Series)

        # DataFrame whose first-level columns are identifiers shall have
        # them in __dir__.
        df = pd.DataFrame(
            [list("abcd"), list("efgh")],
            columns=pd.MultiIndex.from_tuples(list(zip("ABCD", "EFGH"))),
        )
        for key in list("ABCD"):
            assert key in dir(df)
        for key in list("EFGH"):
            assert key not in dir(df)
        assert isinstance(df.__getitem__("A"), pd.DataFrame)

    def test_not_hashable(self):
        empty_frame = DataFrame()

        df = DataFrame([1])
        msg = "'DataFrame' objects are mutable, thus they cannot be hashed"
        with pytest.raises(TypeError, match=msg):
            hash(df)
        with pytest.raises(TypeError, match=msg):
            hash(empty_frame)

    def test_new_empty_index(self):
        df1 = DataFrame(np.random.randn(0, 3))
        df2 = DataFrame(np.random.randn(0, 3))
        df1.index.name = "foo"
        assert df2.index.name is None

    def test_array_interface(self, float_frame):
        with np.errstate(all="ignore"):
            result = np.sqrt(float_frame)
        assert isinstance(result, type(float_frame))
        assert result.index is float_frame.index
        assert result.columns is float_frame.columns

        tm.assert_frame_equal(result, float_frame.apply(np.sqrt))

    def test_get_agg_axis(self, float_frame):
        cols = float_frame._get_agg_axis(0)
        assert cols is float_frame.columns

        idx = float_frame._get_agg_axis(1)
        assert idx is float_frame.index

        msg = r"Axis must be 0 or 1 \(got 2\)"
        with pytest.raises(ValueError, match=msg):
            float_frame._get_agg_axis(2)

    def test_nonzero(self, float_frame, float_string_frame):
        empty_frame = DataFrame()
        assert empty_frame.empty

        assert not float_frame.empty
        assert not float_string_frame.empty

        # corner case
        df = DataFrame({"A": [1.0, 2.0, 3.0], "B": ["a", "b", "c"]}, index=np.arange(3))
        del df["A"]
        assert not df.empty

    def test_iteritems(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=["a", "a", "b"])
        for k, v in df.items():
            assert isinstance(v, DataFrame._constructor_sliced)

    def test_items(self):
        # GH 17213, GH 13918
        cols = ["a", "b", "c"]
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=cols)
        for c, (k, v) in zip(cols, df.items()):
            assert c == k
            assert isinstance(v, Series)
            assert (df[k] == v).all()

    def test_iter(self, float_frame):
        assert tm.equalContents(list(float_frame), float_frame.columns)

    def test_iterrows(self, float_frame, float_string_frame):
        for k, v in float_frame.iterrows():
            exp = float_frame.loc[k]
            tm.assert_series_equal(v, exp)

        for k, v in float_string_frame.iterrows():
            exp = float_string_frame.loc[k]
            tm.assert_series_equal(v, exp)

    def test_iterrows_iso8601(self):
        # GH 19671
        s = DataFrame(
            {
                "non_iso8601": ["M1701", "M1802", "M1903", "M2004"],
                "iso8601": date_range("2000-01-01", periods=4, freq="M"),
            }
        )
        for k, v in s.iterrows():
            exp = s.loc[k]
            tm.assert_series_equal(v, exp)

    def test_iterrows_corner(self):
        # gh-12222
        df = DataFrame(
            {
                "a": [datetime.datetime(2015, 1, 1)],
                "b": [None],
                "c": [None],
                "d": [""],
                "e": [[]],
                "f": [set()],
                "g": [{}],
            }
        )
        expected = Series(
            [datetime.datetime(2015, 1, 1), None, None, "", [], set(), {}],
            index=list("abcdefg"),
            name=0,
            dtype="object",
        )
        _, result = next(df.iterrows())
        tm.assert_series_equal(result, expected)

    def test_itertuples(self, float_frame):
        for i, tup in enumerate(float_frame.itertuples()):
            s = DataFrame._constructor_sliced(tup[1:])
            s.name = tup[0]
            expected = float_frame.iloc[i, :].reset_index(drop=True)
            tm.assert_series_equal(s, expected)

        df = DataFrame(
            {"floats": np.random.randn(5), "ints": range(5)}, columns=["floats", "ints"]
        )

        for tup in df.itertuples(index=False):
            assert isinstance(tup[1], int)

        df = DataFrame(data={"a": [1, 2, 3], "b": [4, 5, 6]})
        dfaa = df[["a", "a"]]

        assert list(dfaa.itertuples()) == [(0, 1, 1), (1, 2, 2), (2, 3, 3)]

        # repr with int on 32-bit/windows
        if not (compat.is_platform_windows() or compat.is_platform_32bit()):
            assert (
                repr(list(df.itertuples(name=None)))
                == "[(0, 1, 4), (1, 2, 5), (2, 3, 6)]"
            )

        tup = next(df.itertuples(name="TestName"))
        assert tup._fields == ("Index", "a", "b")
        assert (tup.Index, tup.a, tup.b) == tup
        assert type(tup).__name__ == "TestName"

        df.columns = ["def", "return"]
        tup2 = next(df.itertuples(name="TestName"))
        assert tup2 == (0, 1, 4)
        assert tup2._fields == ("Index", "_1", "_2")

        df3 = DataFrame({"f" + str(i): [i] for i in range(1024)})
        # will raise SyntaxError if trying to create namedtuple
        tup3 = next(df3.itertuples())
        assert isinstance(tup3, tuple)
        if PY37:
            assert hasattr(tup3, "_fields")
        else:
            assert not hasattr(tup3, "_fields")

        # GH 28282
        df_254_columns = DataFrame([{f"foo_{i}": f"bar_{i}" for i in range(254)}])
        result_254_columns = next(df_254_columns.itertuples(index=False))
        assert isinstance(result_254_columns, tuple)
        assert hasattr(result_254_columns, "_fields")

        df_255_columns = DataFrame([{f"foo_{i}": f"bar_{i}" for i in range(255)}])
        result_255_columns = next(df_255_columns.itertuples(index=False))
        assert isinstance(result_255_columns, tuple)

        # Dataframes with >=255 columns will fallback to regular tuples on python < 3.7
        if PY37:
            assert hasattr(result_255_columns, "_fields")
        else:
            assert not hasattr(result_255_columns, "_fields")

    def test_sequence_like_with_categorical(self):

        # GH 7839
        # make sure can iterate
        df = DataFrame(
            {"id": [1, 2, 3, 4, 5, 6], "raw_grade": ["a", "b", "b", "a", "a", "e"]}
        )
        df["grade"] = Categorical(df["raw_grade"])

        # basic sequencing testing
        result = list(df.grade.values)
        expected = np.array(df.grade.values).tolist()
        tm.assert_almost_equal(result, expected)

        # iteration
        for t in df.itertuples(index=False):
            str(t)

        for row, s in df.iterrows():
            str(s)

        for c, col in df.items():
            str(s)

    def test_len(self, float_frame):
        assert len(float_frame) == len(float_frame.index)

    def test_values_mixed_dtypes(self, float_frame, float_string_frame):
        frame = float_frame
        arr = frame.values

        frame_cols = frame.columns
        for i, row in enumerate(arr):
            for j, value in enumerate(row):
                col = frame_cols[j]
                if np.isnan(value):
                    assert np.isnan(frame[col][i])
                else:
                    assert value == frame[col][i]

        # mixed type
        arr = float_string_frame[["foo", "A"]].values
        assert arr[0, 0] == "bar"

        df = DataFrame({"complex": [1j, 2j, 3j], "real": [1, 2, 3]})
        arr = df.values
        assert arr[0, 0] == 1j

        # single block corner case
        arr = float_frame[["A", "B"]].values
        expected = float_frame.reindex(columns=["A", "B"]).values
        tm.assert_almost_equal(arr, expected)

    def test_to_numpy(self):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4.5]})
        expected = np.array([[1, 3], [2, 4.5]])
        result = df.to_numpy()
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_dtype(self):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4.5]})
        expected = np.array([[1, 3], [2, 4]], dtype="int64")
        result = df.to_numpy(dtype="int64")
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_copy(self):
        arr = np.random.randn(4, 3)
        df = pd.DataFrame(arr)
        assert df.values.base is arr
        assert df.to_numpy(copy=False).base is arr
        assert df.to_numpy(copy=True).base is None

    def test_transpose(self, float_frame):
        frame = float_frame
        dft = frame.T
        for idx, series in dft.items():
            for col, value in series.items():
                if np.isnan(value):
                    assert np.isnan(frame[col][idx])
                else:
                    assert value == frame[col][idx]

        # mixed type
        index, data = tm.getMixedTypeDict()
        mixed = DataFrame(data, index=index)

        mixed_T = mixed.T
        for col, s in mixed_T.items():
            assert s.dtype == np.object_

    def test_swapaxes(self):
        df = DataFrame(np.random.randn(10, 5))
        tm.assert_frame_equal(df.T, df.swapaxes(0, 1))
        tm.assert_frame_equal(df.T, df.swapaxes(1, 0))
        tm.assert_frame_equal(df, df.swapaxes(0, 0))
        msg = (
            "No axis named 2 for object type"
            r" <class 'pandas.core(.sparse)?.frame.(Sparse)?DataFrame'>"
        )
        with pytest.raises(ValueError, match=msg):
            df.swapaxes(2, 5)

    def test_axis_aliases(self, float_frame):
        f = float_frame

        # reg name
        expected = f.sum(axis=0)
        result = f.sum(axis="index")
        tm.assert_series_equal(result, expected)

        expected = f.sum(axis=1)
        result = f.sum(axis="columns")
        tm.assert_series_equal(result, expected)

    def test_class_axis(self):
        # GH 18147
        # no exception and no empty docstring
        assert pydoc.getdoc(DataFrame.index)
        assert pydoc.getdoc(DataFrame.columns)

    def test_more_values(self, float_string_frame):
        values = float_string_frame.values
        assert values.shape[1] == len(float_string_frame.columns)

    def test_repr_with_mi_nat(self, float_string_frame):
        df = DataFrame(
            {"X": [1, 2]}, index=[[pd.NaT, pd.Timestamp("20130101")], ["a", "b"]]
        )
        result = repr(df)
        expected = "              X\nNaT        a  1\n2013-01-01 b  2"
        assert result == expected

    def test_items_names(self, float_string_frame):
        for k, v in float_string_frame.items():
            assert v.name == k

    def test_series_put_names(self, float_string_frame):
        series = float_string_frame._series
        for k, v in series.items():
            assert v.name == k

    def test_empty_nonzero(self):
        df = DataFrame([1, 2, 3])
        assert not df.empty
        df = DataFrame(index=[1], columns=[1])
        assert not df.empty
        df = DataFrame(index=["a", "b"], columns=["c", "d"]).dropna()
        assert df.empty
        assert df.T.empty
        empty_frames = [
            DataFrame(),
            DataFrame(index=[1]),
            DataFrame(columns=[1]),
            DataFrame({1: []}),
        ]
        for df in empty_frames:
            assert df.empty
            assert df.T.empty

    def test_with_datetimelikes(self):

        df = DataFrame(
            {
                "A": date_range("20130101", periods=10),
                "B": timedelta_range("1 day", periods=10),
            }
        )
        t = df.T

        result = t.dtypes.value_counts()
        expected = Series({np.dtype("object"): 10})
        tm.assert_series_equal(result, expected)

    def test_values(self, float_frame):
        float_frame.values[:, 0] = 5.0
        assert (float_frame.values[:, 0] == 5).all()

    def test_deepcopy(self, float_frame):
        cp = deepcopy(float_frame)
        series = cp["A"]
        series[:] = 10
        for idx, value in series.items():
            assert float_frame["A"][idx] != value

    def test_transpose_get_view(self, float_frame):
        dft = float_frame.T
        dft.values[:, 5:10] = 5

        assert (float_frame.values[5:10] == 5).all()

    def test_inplace_return_self(self):
        # GH 1893

        data = DataFrame(
            {"a": ["foo", "bar", "baz", "qux"], "b": [0, 0, 1, 1], "c": [1, 2, 3, 4]}
        )

        def _check_f(base, f):
            result = f(base)
            assert result is None

        # -----DataFrame-----

        # set_index
        f = lambda x: x.set_index("a", inplace=True)
        _check_f(data.copy(), f)

        # reset_index
        f = lambda x: x.reset_index(inplace=True)
        _check_f(data.set_index("a"), f)

        # drop_duplicates
        f = lambda x: x.drop_duplicates(inplace=True)
        _check_f(data.copy(), f)

        # sort
        f = lambda x: x.sort_values("b", inplace=True)
        _check_f(data.copy(), f)

        # sort_index
        f = lambda x: x.sort_index(inplace=True)
        _check_f(data.copy(), f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(data.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(data.copy(), f)

        # rename
        f = lambda x: x.rename({1: "foo"}, inplace=True)
        _check_f(data.copy(), f)

        # -----Series-----
        d = data.copy()["c"]

        # reset_index
        f = lambda x: x.reset_index(inplace=True, drop=True)
        _check_f(data.set_index("a")["c"], f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(d.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(d.copy(), f)

        # rename
        f = lambda x: x.rename({1: "foo"}, inplace=True)
        _check_f(d.copy(), f)

    def test_tab_complete_warning(self, ip):
        # GH 16409
        pytest.importorskip("IPython", minversion="6.0.0")
        from IPython.core.completer import provisionalcompleter

        code = "import pandas as pd; df = pd.DataFrame()"
        ip.run_code(code)
        with tm.assert_produces_warning(None):
            with provisionalcompleter("ignore"):
                list(ip.Completer.completions("df.", 1))
