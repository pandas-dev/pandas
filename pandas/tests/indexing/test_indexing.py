""" test fancy indexing & misc """

from datetime import datetime
import re
import weakref

import numpy as np
import pytest

from pandas.errors import AbstractMethodError

from pandas.core.dtypes.common import is_float_dtype, is_integer_dtype

import pandas as pd
from pandas import DataFrame, Index, NaT, Series
from pandas.core.generic import NDFrame
from pandas.core.indexers import validate_indices
from pandas.core.indexing import _maybe_numeric_slice, _non_reducing_slice
from pandas.tests.indexing.common import Base, _mklbl
import pandas.util.testing as tm

# ------------------------------------------------------------------------
# Indexing test cases


class TestFancy(Base):
    """ pure get/set item & fancy indexing """

    def test_setitem_ndarray_1d(self):
        # GH5508

        # len of indexer vs length of the 1d ndarray
        df = DataFrame(index=Index(np.arange(1, 11)))
        df["foo"] = np.zeros(10, dtype=np.float64)
        df["bar"] = np.zeros(10, dtype=np.complex)

        # invalid
        with pytest.raises(ValueError):
            df.loc[df.index[2:5], "bar"] = np.array([2.33j, 1.23 + 0.1j, 2.2, 1.0])

        # valid
        df.loc[df.index[2:6], "bar"] = np.array([2.33j, 1.23 + 0.1j, 2.2, 1.0])

        result = df.loc[df.index[2:6], "bar"]
        expected = Series(
            [2.33j, 1.23 + 0.1j, 2.2, 1.0], index=[3, 4, 5, 6], name="bar"
        )
        tm.assert_series_equal(result, expected)

        # dtype getting changed?
        df = DataFrame(index=Index(np.arange(1, 11)))
        df["foo"] = np.zeros(10, dtype=np.float64)
        df["bar"] = np.zeros(10, dtype=np.complex)

        with pytest.raises(ValueError):
            df[2:5] = np.arange(1, 4) * 1j

    @pytest.mark.parametrize(
        "index", tm.all_index_generator(5), ids=lambda x: type(x).__name__
    )
    @pytest.mark.parametrize(
        "obj",
        [
            lambda i: Series(np.arange(len(i)), index=i),
            lambda i: DataFrame(np.random.randn(len(i), len(i)), index=i, columns=i),
        ],
        ids=["Series", "DataFrame"],
    )
    @pytest.mark.parametrize(
        "idxr, idxr_id",
        [
            (lambda x: x, "getitem"),
            (lambda x: x.loc, "loc"),
            (lambda x: x.iloc, "iloc"),
        ],
    )
    def test_getitem_ndarray_3d(self, index, obj, idxr, idxr_id):
        # GH 25567
        obj = obj(index)
        idxr = idxr(obj)
        nd3 = np.random.randint(5, size=(2, 2, 2))

        msg = (
            r"Buffer has wrong number of dimensions \(expected 1,"
            r" got 3\)|"
            "The truth value of an array with more than one element is"
            " ambiguous|"
            "Cannot index with multidimensional key|"
            r"Wrong number of dimensions. values.ndim != ndim \[3 != 1\]|"
            "No matching signature found|"  # TypeError
            "unhashable type: 'numpy.ndarray'"  # TypeError
        )

        if (
            isinstance(obj, Series)
            and idxr_id == "getitem"
            and index.inferred_type
            in [
                "string",
                "datetime64",
                "period",
                "timedelta64",
                "boolean",
                "categorical",
            ]
        ):
            idxr[nd3]
        else:
            if (
                isinstance(obj, DataFrame)
                and idxr_id == "getitem"
                and index.inferred_type == "boolean"
            ):
                error = TypeError
            elif idxr_id == "getitem" and index.inferred_type == "interval":
                error = TypeError
            else:
                error = ValueError

            with pytest.raises(error, match=msg):
                idxr[nd3]

    @pytest.mark.parametrize(
        "index", tm.all_index_generator(5), ids=lambda x: type(x).__name__
    )
    @pytest.mark.parametrize(
        "obj",
        [
            lambda i: Series(np.arange(len(i)), index=i),
            lambda i: DataFrame(np.random.randn(len(i), len(i)), index=i, columns=i),
        ],
        ids=["Series", "DataFrame"],
    )
    @pytest.mark.parametrize(
        "idxr, idxr_id",
        [
            (lambda x: x, "setitem"),
            (lambda x: x.loc, "loc"),
            (lambda x: x.iloc, "iloc"),
        ],
    )
    def test_setitem_ndarray_3d(self, index, obj, idxr, idxr_id):
        # GH 25567
        obj = obj(index)
        idxr = idxr(obj)
        nd3 = np.random.randint(5, size=(2, 2, 2))

        msg = (
            r"Buffer has wrong number of dimensions \(expected 1,"
            r" got 3\)|"
            "The truth value of an array with more than one element is"
            " ambiguous|"
            "Only 1-dimensional input arrays are supported|"
            "'pandas._libs.interval.IntervalTree' object has no attribute"
            " 'set_value'|"  # AttributeError
            "unhashable type: 'numpy.ndarray'|"  # TypeError
            "No matching signature found|"  # TypeError
            r"^\[\[\["  # pandas.core.indexing.IndexingError
        )

        if (idxr_id == "iloc") or (
            (
                isinstance(obj, Series)
                and idxr_id == "setitem"
                and index.inferred_type
                in [
                    "floating",
                    "string",
                    "datetime64",
                    "period",
                    "timedelta64",
                    "boolean",
                    "categorical",
                ]
            )
        ):
            idxr[nd3] = 0
        else:
            with pytest.raises(
                (ValueError, AttributeError, TypeError, pd.core.indexing.IndexingError),
                match=msg,
            ):
                idxr[nd3] = 0

    def test_inf_upcast(self):
        # GH 16957
        # We should be able to use np.inf as a key
        # np.inf should cause an index to convert to float

        # Test with np.inf in rows
        df = DataFrame(columns=[0])
        df.loc[1] = 1
        df.loc[2] = 2
        df.loc[np.inf] = 3

        # make sure we can look up the value
        assert df.loc[np.inf, 0] == 3

        result = df.index
        expected = pd.Float64Index([1, 2, np.inf])
        tm.assert_index_equal(result, expected)

        # Test with np.inf in columns
        df = DataFrame()
        df.loc[0, 0] = 1
        df.loc[1, 1] = 2
        df.loc[0, np.inf] = 3

        result = df.columns
        expected = pd.Float64Index([0, 1, np.inf])
        tm.assert_index_equal(result, expected)

    def test_setitem_dtype_upcast(self):

        # GH3216
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df["c"] = np.nan
        assert df["c"].dtype == np.float64

        df.loc[0, "c"] = "foo"
        expected = DataFrame(
            [{"a": 1, "b": np.nan, "c": "foo"}, {"a": 3, "b": 2, "c": np.nan}]
        )
        tm.assert_frame_equal(df, expected)

        # GH10280
        df = DataFrame(
            np.arange(6, dtype="int64").reshape(2, 3),
            index=list("ab"),
            columns=["foo", "bar", "baz"],
        )

        for val in [3.14, "wxyz"]:
            left = df.copy()
            left.loc["a", "bar"] = val
            right = DataFrame(
                [[0, val, 2], [3, 4, 5]],
                index=list("ab"),
                columns=["foo", "bar", "baz"],
            )

            tm.assert_frame_equal(left, right)
            assert is_integer_dtype(left["foo"])
            assert is_integer_dtype(left["baz"])

        left = DataFrame(
            np.arange(6, dtype="int64").reshape(2, 3) / 10.0,
            index=list("ab"),
            columns=["foo", "bar", "baz"],
        )
        left.loc["a", "bar"] = "wxyz"

        right = DataFrame(
            [[0, "wxyz", 0.2], [0.3, 0.4, 0.5]],
            index=list("ab"),
            columns=["foo", "bar", "baz"],
        )

        tm.assert_frame_equal(left, right)
        assert is_float_dtype(left["foo"])
        assert is_float_dtype(left["baz"])

    def test_dups_fancy_indexing(self):

        # GH 3455

        df = tm.makeCustomDataframe(10, 3)
        df.columns = ["a", "a", "b"]
        result = df[["b", "a"]].columns
        expected = Index(["b", "a", "a"])
        tm.assert_index_equal(result, expected)

        # across dtypes
        df = DataFrame([[1, 2, 1.0, 2.0, 3.0, "foo", "bar"]], columns=list("aaaaaaa"))
        df.head()
        str(df)
        result = DataFrame([[1, 2, 1.0, 2.0, 3.0, "foo", "bar"]])
        result.columns = list("aaaaaaa")

        # TODO(wesm): unused?
        df_v = df.iloc[:, 4]  # noqa
        res_v = result.iloc[:, 4]  # noqa

        tm.assert_frame_equal(df, result)

        # GH 3561, dups not in selected order
        df = DataFrame(
            {"test": [5, 7, 9, 11], "test1": [4.0, 5, 6, 7], "other": list("abcd")},
            index=["A", "A", "B", "C"],
        )
        rows = ["C", "B"]
        expected = DataFrame(
            {"test": [11, 9], "test1": [7.0, 6], "other": ["d", "c"]}, index=rows
        )
        result = df.loc[rows]
        tm.assert_frame_equal(result, expected)

        result = df.loc[Index(rows)]
        tm.assert_frame_equal(result, expected)

        rows = ["C", "B", "E"]
        with pytest.raises(KeyError, match="with any missing labels"):
            df.loc[rows]

        # see GH5553, make sure we use the right indexer
        rows = ["F", "G", "H", "C", "B", "E"]
        with pytest.raises(KeyError, match="with any missing labels"):
            df.loc[rows]

        # List containing only missing label
        dfnu = DataFrame(np.random.randn(5, 3), index=list("AABCD"))
        with pytest.raises(
            KeyError,
            match=re.escape(
                "\"None of [Index(['E'], dtype='object')] are in the [index]\""
            ),
        ):
            dfnu.loc[["E"]]

        # ToDo: check_index_type can be True after GH 11497

        # GH 4619; duplicate indexer with missing label
        df = DataFrame({"A": [0, 1, 2]})
        with pytest.raises(KeyError, match="with any missing labels"):
            df.loc[[0, 8, 0]]

        df = DataFrame({"A": list("abc")})
        with pytest.raises(KeyError, match="with any missing labels"):
            df.loc[[0, 8, 0]]

        # non unique with non unique selector
        df = DataFrame({"test": [5, 7, 9, 11]}, index=["A", "A", "B", "C"])
        with pytest.raises(KeyError, match="with any missing labels"):
            df.loc[["A", "A", "E"]]

    def test_dups_fancy_indexing2(self):
        # GH 5835
        # dups on index and missing values
        df = DataFrame(np.random.randn(5, 5), columns=["A", "B", "B", "B", "A"])

        with pytest.raises(KeyError, match="with any missing labels"):
            df.loc[:, ["A", "B", "C"]]

        # GH 6504, multi-axis indexing
        df = DataFrame(
            np.random.randn(9, 2), index=[1, 1, 1, 2, 2, 2, 3, 3, 3], columns=["a", "b"]
        )

        expected = df.iloc[0:6]
        result = df.loc[[1, 2]]
        tm.assert_frame_equal(result, expected)

        expected = df
        result = df.loc[:, ["a", "b"]]
        tm.assert_frame_equal(result, expected)

        expected = df.iloc[0:6, :]
        result = df.loc[[1, 2], ["a", "b"]]
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("case", [lambda s: s, lambda s: s.loc])
    def test_duplicate_int_indexing(self, case):
        # GH 17347
        s = pd.Series(range(3), index=[1, 1, 3])
        expected = s[1]
        result = case(s)[[1]]
        tm.assert_series_equal(result, expected)

    def test_indexing_mixed_frame_bug(self):

        # GH3492
        df = DataFrame(
            {"a": {1: "aaa", 2: "bbb", 3: "ccc"}, "b": {1: 111, 2: 222, 3: 333}}
        )

        # this works, new column is created correctly
        df["test"] = df["a"].apply(lambda x: "_" if x == "aaa" else x)

        # this does not work, ie column test is not changed
        idx = df["test"] == "_"
        temp = df.loc[idx, "a"].apply(lambda x: "-----" if x == "aaa" else x)
        df.loc[idx, "test"] = temp
        assert df.iloc[0, 2] == "-----"

    def test_multitype_list_index_access(self):
        # GH 10610
        df = DataFrame(np.random.random((10, 5)), columns=["a"] + [20, 21, 22, 23])

        with pytest.raises(KeyError, match=re.escape("'[-8, 26] not in index'")):
            df[[22, 26, -8]]
        assert df[21].shape[0] == df.shape[0]

    def test_set_index_nan(self):

        # GH 3586
        df = DataFrame(
            {
                "PRuid": {
                    17: "nonQC",
                    18: "nonQC",
                    19: "nonQC",
                    20: "10",
                    21: "11",
                    22: "12",
                    23: "13",
                    24: "24",
                    25: "35",
                    26: "46",
                    27: "47",
                    28: "48",
                    29: "59",
                    30: "10",
                },
                "QC": {
                    17: 0.0,
                    18: 0.0,
                    19: 0.0,
                    20: np.nan,
                    21: np.nan,
                    22: np.nan,
                    23: np.nan,
                    24: 1.0,
                    25: np.nan,
                    26: np.nan,
                    27: np.nan,
                    28: np.nan,
                    29: np.nan,
                    30: np.nan,
                },
                "data": {
                    17: 7.9544899999999998,
                    18: 8.0142609999999994,
                    19: 7.8591520000000008,
                    20: 0.86140349999999999,
                    21: 0.87853110000000001,
                    22: 0.8427041999999999,
                    23: 0.78587700000000005,
                    24: 0.73062459999999996,
                    25: 0.81668560000000001,
                    26: 0.81927080000000008,
                    27: 0.80705009999999999,
                    28: 0.81440240000000008,
                    29: 0.80140849999999997,
                    30: 0.81307740000000006,
                },
                "year": {
                    17: 2006,
                    18: 2007,
                    19: 2008,
                    20: 1985,
                    21: 1985,
                    22: 1985,
                    23: 1985,
                    24: 1985,
                    25: 1985,
                    26: 1985,
                    27: 1985,
                    28: 1985,
                    29: 1985,
                    30: 1986,
                },
            }
        ).reset_index()

        result = (
            df.set_index(["year", "PRuid", "QC"])
            .reset_index()
            .reindex(columns=df.columns)
        )
        tm.assert_frame_equal(result, df)

    def test_multi_assign(self):

        # GH 3626, an assignment of a sub-df to a df
        df = DataFrame(
            {
                "FC": ["a", "b", "a", "b", "a", "b"],
                "PF": [0, 0, 0, 0, 1, 1],
                "col1": list(range(6)),
                "col2": list(range(6, 12)),
            }
        )
        df.iloc[1, 0] = np.nan
        df2 = df.copy()

        mask = ~df2.FC.isna()
        cols = ["col1", "col2"]

        dft = df2 * 2
        dft.iloc[3, 3] = np.nan

        expected = DataFrame(
            {
                "FC": ["a", np.nan, "a", "b", "a", "b"],
                "PF": [0, 0, 0, 0, 1, 1],
                "col1": Series([0, 1, 4, 6, 8, 10]),
                "col2": [12, 7, 16, np.nan, 20, 22],
            }
        )

        # frame on rhs
        df2.loc[mask, cols] = dft.loc[mask, cols]
        tm.assert_frame_equal(df2, expected)

        df2.loc[mask, cols] = dft.loc[mask, cols]
        tm.assert_frame_equal(df2, expected)

        # with an ndarray on rhs
        # coerces to float64 because values has float64 dtype
        # GH 14001
        expected = DataFrame(
            {
                "FC": ["a", np.nan, "a", "b", "a", "b"],
                "PF": [0, 0, 0, 0, 1, 1],
                "col1": [0.0, 1.0, 4.0, 6.0, 8.0, 10.0],
                "col2": [12, 7, 16, np.nan, 20, 22],
            }
        )
        df2 = df.copy()
        df2.loc[mask, cols] = dft.loc[mask, cols].values
        tm.assert_frame_equal(df2, expected)
        df2.loc[mask, cols] = dft.loc[mask, cols].values
        tm.assert_frame_equal(df2, expected)

        # broadcasting on the rhs is required
        df = DataFrame(
            dict(
                A=[1, 2, 0, 0, 0],
                B=[0, 0, 0, 10, 11],
                C=[0, 0, 0, 10, 11],
                D=[3, 4, 5, 6, 7],
            )
        )

        expected = df.copy()
        mask = expected["A"] == 0
        for col in ["A", "B"]:
            expected.loc[mask, col] = df["D"]

        df.loc[df["A"] == 0, ["A", "B"]] = df["D"]
        tm.assert_frame_equal(df, expected)

    def test_setitem_list(self):

        # GH 6043
        # iloc with a list
        df = DataFrame(index=[0, 1], columns=[0])
        df.iloc[1, 0] = [1, 2, 3]
        df.iloc[1, 0] = [1, 2]

        result = DataFrame(index=[0, 1], columns=[0])
        result.iloc[1, 0] = [1, 2]

        tm.assert_frame_equal(result, df)

        # iloc with an object
        class TO:
            def __init__(self, value):
                self.value = value

            def __str__(self) -> str:
                return "[{0}]".format(self.value)

            __repr__ = __str__

            def __eq__(self, other) -> bool:
                return self.value == other.value

            def view(self):
                return self

        df = DataFrame(index=[0, 1], columns=[0])
        df.iloc[1, 0] = TO(1)
        df.iloc[1, 0] = TO(2)

        result = DataFrame(index=[0, 1], columns=[0])
        result.iloc[1, 0] = TO(2)

        tm.assert_frame_equal(result, df)

        # remains object dtype even after setting it back
        df = DataFrame(index=[0, 1], columns=[0])
        df.iloc[1, 0] = TO(1)
        df.iloc[1, 0] = np.nan
        result = DataFrame(index=[0, 1], columns=[0])

        tm.assert_frame_equal(result, df)

    def test_string_slice(self):
        # GH 14424
        # string indexing against datetimelike with object
        # dtype should properly raises KeyError
        df = DataFrame([1], Index([pd.Timestamp("2011-01-01")], dtype=object))
        assert df.index.is_all_dates
        with pytest.raises(KeyError, match="'2011'"):
            df["2011"]

        with pytest.raises(KeyError, match="'2011'"):
            df.loc["2011", 0]

        df = DataFrame()
        assert not df.index.is_all_dates
        with pytest.raises(KeyError, match="'2011'"):
            df["2011"]

        with pytest.raises(KeyError, match="'2011'"):
            df.loc["2011", 0]

    def test_astype_assignment(self):

        # GH4312 (iloc)
        df_orig = DataFrame(
            [["1", "2", "3", ".4", 5, 6.0, "foo"]], columns=list("ABCDEFG")
        )

        df = df_orig.copy()
        df.iloc[:, 0:2] = df.iloc[:, 0:2].astype(np.int64)
        expected = DataFrame(
            [[1, 2, "3", ".4", 5, 6.0, "foo"]], columns=list("ABCDEFG")
        )
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.iloc[:, 0:2] = df.iloc[:, 0:2]._convert(datetime=True, numeric=True)
        expected = DataFrame(
            [[1, 2, "3", ".4", 5, 6.0, "foo"]], columns=list("ABCDEFG")
        )
        tm.assert_frame_equal(df, expected)

        # GH5702 (loc)
        df = df_orig.copy()
        df.loc[:, "A"] = df.loc[:, "A"].astype(np.int64)
        expected = DataFrame(
            [[1, "2", "3", ".4", 5, 6.0, "foo"]], columns=list("ABCDEFG")
        )
        tm.assert_frame_equal(df, expected)

        df = df_orig.copy()
        df.loc[:, ["B", "C"]] = df.loc[:, ["B", "C"]].astype(np.int64)
        expected = DataFrame(
            [["1", 2, 3, ".4", 5, 6.0, "foo"]], columns=list("ABCDEFG")
        )
        tm.assert_frame_equal(df, expected)

        # full replacements / no nans
        df = DataFrame({"A": [1.0, 2.0, 3.0, 4.0]})
        df.iloc[:, 0] = df["A"].astype(np.int64)
        expected = DataFrame({"A": [1, 2, 3, 4]})
        tm.assert_frame_equal(df, expected)

        df = DataFrame({"A": [1.0, 2.0, 3.0, 4.0]})
        df.loc[:, "A"] = df["A"].astype(np.int64)
        expected = DataFrame({"A": [1, 2, 3, 4]})
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "index,val",
        [
            (Index([0, 1, 2]), 2),
            (Index([0, 1, "2"]), "2"),
            (Index([0, 1, 2, np.inf, 4]), 4),
            (Index([0, 1, 2, np.nan, 4]), 4),
            (Index([0, 1, 2, np.inf]), np.inf),
            (Index([0, 1, 2, np.nan]), np.nan),
        ],
    )
    def test_index_contains(self, index, val):
        assert val in index

    @pytest.mark.parametrize(
        "index,val",
        [
            (Index([0, 1, 2]), "2"),
            (Index([0, 1, "2"]), 2),
            (Index([0, 1, 2, np.inf]), 4),
            (Index([0, 1, 2, np.nan]), 4),
            (Index([0, 1, 2, np.inf]), np.nan),
            (Index([0, 1, 2, np.nan]), np.inf),
            # Checking if np.inf in Int64Index should not cause an OverflowError
            # Related to GH 16957
            (pd.Int64Index([0, 1, 2]), np.inf),
            (pd.Int64Index([0, 1, 2]), np.nan),
            (pd.UInt64Index([0, 1, 2]), np.inf),
            (pd.UInt64Index([0, 1, 2]), np.nan),
        ],
    )
    def test_index_not_contains(self, index, val):
        assert val not in index

    @pytest.mark.parametrize(
        "index,val", [(Index([0, 1, "2"]), 0), (Index([0, 1, "2"]), "2")]
    )
    def test_mixed_index_contains(self, index, val):
        # GH 19860
        assert val in index

    @pytest.mark.parametrize(
        "index,val", [(Index([0, 1, "2"]), "1"), (Index([0, 1, "2"]), 2)]
    )
    def test_mixed_index_not_contains(self, index, val):
        # GH 19860
        assert val not in index

    def test_contains_with_float_index(self):
        # GH#22085
        integer_index = pd.Int64Index([0, 1, 2, 3])
        uinteger_index = pd.UInt64Index([0, 1, 2, 3])
        float_index = pd.Float64Index([0.1, 1.1, 2.2, 3.3])

        for index in (integer_index, uinteger_index):
            assert 1.1 not in index
            assert 1.0 in index
            assert 1 in index

        assert 1.1 in float_index
        assert 1.0 not in float_index
        assert 1 not in float_index

    def test_index_type_coercion(self):

        # GH 11836
        # if we have an index type and set it with something that looks
        # to numpy like the same, but is actually, not
        # (e.g. setting with a float or string '0')
        # then we need to coerce to object

        # integer indexes
        for s in [Series(range(5)), Series(range(5), index=range(1, 6))]:

            assert s.index.is_integer()

            for indexer in [lambda x: x.loc, lambda x: x]:
                s2 = s.copy()
                indexer(s2)[0.1] = 0
                assert s2.index.is_floating()
                assert indexer(s2)[0.1] == 0

                s2 = s.copy()
                indexer(s2)[0.0] = 0
                exp = s.index
                if 0 not in s:
                    exp = Index(s.index.tolist() + [0])
                tm.assert_index_equal(s2.index, exp)

                s2 = s.copy()
                indexer(s2)["0"] = 0
                assert s2.index.is_object()

        for s in [Series(range(5), index=np.arange(5.0))]:

            assert s.index.is_floating()

            for idxr in [lambda x: x.loc, lambda x: x]:

                s2 = s.copy()
                idxr(s2)[0.1] = 0
                assert s2.index.is_floating()
                assert idxr(s2)[0.1] == 0

                s2 = s.copy()
                idxr(s2)[0.0] = 0
                tm.assert_index_equal(s2.index, s.index)

                s2 = s.copy()
                idxr(s2)["0"] = 0
                assert s2.index.is_object()


class TestMisc(Base):
    def test_float_index_to_mixed(self):
        df = DataFrame({0.0: np.random.rand(10), 1.0: np.random.rand(10)})
        df["a"] = 10
        tm.assert_frame_equal(
            DataFrame({0.0: df[0.0], 1.0: df[1.0], "a": [10] * 10}), df
        )

    def test_float_index_non_scalar_assignment(self):
        df = DataFrame({"a": [1, 2, 3], "b": [3, 4, 5]}, index=[1.0, 2.0, 3.0])
        df.loc[df.index[:2]] = 1
        expected = DataFrame({"a": [1, 1, 3], "b": [1, 1, 5]}, index=df.index)
        tm.assert_frame_equal(expected, df)

        df = DataFrame({"a": [1, 2, 3], "b": [3, 4, 5]}, index=[1.0, 2.0, 3.0])
        df2 = df.copy()
        df.loc[df.index] = df.loc[df.index]
        tm.assert_frame_equal(df, df2)

    def test_float_index_at_iat(self):
        s = Series([1, 2, 3], index=[0.1, 0.2, 0.3])
        for el, item in s.items():
            assert s.at[el] == item
        for i in range(len(s)):
            assert s.iat[i] == i + 1

    def test_mixed_index_assignment(self):
        # GH 19860
        s = Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        s.at["a"] = 11
        assert s.iat[0] == 11
        s.at[1] = 22
        assert s.iat[3] == 22

    def test_mixed_index_no_fallback(self):
        # GH 19860
        s = Series([1, 2, 3, 4, 5], index=["a", "b", "c", 1, 2])
        with pytest.raises(KeyError, match="^0$"):
            s.at[0]
        with pytest.raises(KeyError, match="^4$"):
            s.at[4]

    def test_rhs_alignment(self):
        # GH8258, tests that both rows & columns are aligned to what is
        # assigned to. covers both uniform data-type & multi-type cases
        def run_tests(df, rhs, right):
            # label, index, slice
            lbl_one, idx_one, slice_one = list("bcd"), [1, 2, 3], slice(1, 4)
            lbl_two, idx_two, slice_two = ["joe", "jolie"], [1, 2], slice(1, 3)

            left = df.copy()
            left.loc[lbl_one, lbl_two] = rhs
            tm.assert_frame_equal(left, right)

            left = df.copy()
            left.iloc[idx_one, idx_two] = rhs
            tm.assert_frame_equal(left, right)

            left = df.copy()
            left.iloc[slice_one, slice_two] = rhs
            tm.assert_frame_equal(left, right)

        xs = np.arange(20).reshape(5, 4)
        cols = ["jim", "joe", "jolie", "joline"]
        df = DataFrame(xs, columns=cols, index=list("abcde"))

        # right hand side; permute the indices and multiplpy by -2
        rhs = -2 * df.iloc[3:0:-1, 2:0:-1]

        # expected `right` result; just multiply by -2
        right = df.copy()
        right.iloc[1:4, 1:3] *= -2

        # run tests with uniform dtypes
        run_tests(df, rhs, right)

        # make frames multi-type & re-run tests
        for frame in [df, rhs, right]:
            frame["joe"] = frame["joe"].astype("float64")
            frame["jolie"] = frame["jolie"].map("@{0}".format)

        run_tests(df, rhs, right)

    def test_str_label_slicing_with_negative_step(self):
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            tm.assert_series_equal(s.loc[l_slc], s.iloc[i_slc])

            if not idx.is_integer:
                # For integer indices, .loc and plain getitem are position-based.
                tm.assert_series_equal(s[l_slc], s.iloc[i_slc])
                tm.assert_series_equal(s.loc[l_slc], s.iloc[i_slc])

        for idx in [_mklbl("A", 20), np.arange(20) + 100, np.linspace(100, 150, 20)]:
            idx = Index(idx)
            s = Series(np.arange(20), index=idx)
            assert_slices_equivalent(SLC[idx[9] :: -1], SLC[9::-1])
            assert_slices_equivalent(SLC[: idx[9] : -1], SLC[:8:-1])
            assert_slices_equivalent(SLC[idx[13] : idx[9] : -1], SLC[13:8:-1])
            assert_slices_equivalent(SLC[idx[9] : idx[13] : -1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        s = Series(np.arange(20), index=_mklbl("A", 20))
        with pytest.raises(ValueError, match="slice step cannot be zero"):
            s[::0]
        with pytest.raises(ValueError, match="slice step cannot be zero"):
            s.loc[::0]

    def test_indexing_assignment_dict_already_exists(self):
        df = DataFrame({"x": [1, 2, 6], "y": [2, 2, 8], "z": [-5, 0, 5]}).set_index("z")
        expected = df.copy()
        rhs = dict(x=9, y=99)
        df.loc[5] = rhs
        expected.loc[5] = [9, 99]
        tm.assert_frame_equal(df, expected)

    def test_indexing_dtypes_on_empty(self):
        # Check that .iloc returns correct dtypes GH9983
        df = DataFrame({"a": [1, 2, 3], "b": ["b", "b2", "b3"]})
        df2 = df.iloc[[], :]

        assert df2.loc[:, "a"].dtype == np.int64
        tm.assert_series_equal(df2.loc[:, "a"], df2.iloc[:, 0])

    def test_range_in_series_indexing(self):
        # range can cause an indexing error
        # GH 11652
        for x in [5, 999999, 1000000]:
            s = Series(index=range(x), dtype=np.float64)
            s.loc[range(1)] = 42
            tm.assert_series_equal(s.loc[range(1)], Series(42.0, index=[0]))

            s.loc[range(2)] = 43
            tm.assert_series_equal(s.loc[range(2)], Series(43.0, index=[0, 1]))

    def test_non_reducing_slice(self):
        df = DataFrame([[0, 1], [2, 3]])

        slices = [
            # pd.IndexSlice[:, :],
            pd.IndexSlice[:, 1],
            pd.IndexSlice[1, :],
            pd.IndexSlice[[1], [1]],
            pd.IndexSlice[1, [1]],
            pd.IndexSlice[[1], 1],
            pd.IndexSlice[1],
            pd.IndexSlice[1, 1],
            slice(None, None, None),
            [0, 1],
            np.array([0, 1]),
            Series([0, 1]),
        ]
        for slice_ in slices:
            tslice_ = _non_reducing_slice(slice_)
            assert isinstance(df.loc[tslice_], DataFrame)

    def test_list_slice(self):
        # like dataframe getitem
        slices = [["A"], Series(["A"]), np.array(["A"])]
        df = DataFrame({"A": [1, 2], "B": [3, 4]}, index=["A", "B"])
        expected = pd.IndexSlice[:, ["A"]]
        for subset in slices:
            result = _non_reducing_slice(subset)
            tm.assert_frame_equal(df.loc[result], df.loc[expected])

    def test_maybe_numeric_slice(self):
        df = DataFrame({"A": [1, 2], "B": ["c", "d"], "C": [True, False]})
        result = _maybe_numeric_slice(df, slice_=None)
        expected = pd.IndexSlice[:, ["A"]]
        assert result == expected

        result = _maybe_numeric_slice(df, None, include_bool=True)
        expected = pd.IndexSlice[:, ["A", "C"]]
        result = _maybe_numeric_slice(df, [1])
        expected = [1]
        assert result == expected

    def test_partial_boolean_frame_indexing(self):
        # GH 17170
        df = DataFrame(
            np.arange(9.0).reshape(3, 3), index=list("abc"), columns=list("ABC")
        )
        index_df = DataFrame(1, index=list("ab"), columns=list("AB"))
        result = df[index_df.notnull()]
        expected = DataFrame(
            np.array([[0.0, 1.0, np.nan], [3.0, 4.0, np.nan], [np.nan] * 3]),
            index=list("abc"),
            columns=list("ABC"),
        )
        tm.assert_frame_equal(result, expected)

    def test_no_reference_cycle(self):
        df = DataFrame({"a": [0, 1], "b": [2, 3]})
        for name in ("loc", "iloc", "at", "iat"):
            getattr(df, name)
        wr = weakref.ref(df)
        del df
        assert wr() is None


class TestSeriesNoneCoercion:
    EXPECTED_RESULTS = [
        # For numeric series, we should coerce to NaN.
        ([1, 2, 3], [np.nan, 2, 3]),
        ([1.0, 2.0, 3.0], [np.nan, 2.0, 3.0]),
        # For datetime series, we should coerce to NaT.
        (
            [datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
            [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)],
        ),
        # For objects, we should preserve the None value.
        (["foo", "bar", "baz"], [None, "bar", "baz"]),
    ]

    def test_coercion_with_setitem(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series[0] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)

    def test_coercion_with_loc_setitem(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series.loc[0] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)

    def test_coercion_with_setitem_and_series(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series[start_series == start_series[0]] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)

    def test_coercion_with_loc_and_series(self):
        for start_data, expected_result in self.EXPECTED_RESULTS:
            start_series = Series(start_data)
            start_series.loc[start_series == start_series[0]] = None

            expected_series = Series(expected_result)
            tm.assert_series_equal(start_series, expected_series)


class TestDataframeNoneCoercion:
    EXPECTED_SINGLE_ROW_RESULTS = [
        # For numeric series, we should coerce to NaN.
        ([1, 2, 3], [np.nan, 2, 3]),
        ([1.0, 2.0, 3.0], [np.nan, 2.0, 3.0]),
        # For datetime series, we should coerce to NaT.
        (
            [datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
            [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)],
        ),
        # For objects, we should preserve the None value.
        (["foo", "bar", "baz"], [None, "bar", "baz"]),
    ]

    def test_coercion_with_loc(self):
        for start_data, expected_result in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({"foo": start_data})
            start_dataframe.loc[0, ["foo"]] = None

            expected_dataframe = DataFrame({"foo": expected_result})
            tm.assert_frame_equal(start_dataframe, expected_dataframe)

    def test_coercion_with_setitem_and_dataframe(self):
        for start_data, expected_result in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({"foo": start_data})
            start_dataframe[start_dataframe["foo"] == start_dataframe["foo"][0]] = None

            expected_dataframe = DataFrame({"foo": expected_result})
            tm.assert_frame_equal(start_dataframe, expected_dataframe)

    def test_none_coercion_loc_and_dataframe(self):
        for start_data, expected_result in self.EXPECTED_SINGLE_ROW_RESULTS:
            start_dataframe = DataFrame({"foo": start_data})
            start_dataframe.loc[
                start_dataframe["foo"] == start_dataframe["foo"][0]
            ] = None

            expected_dataframe = DataFrame({"foo": expected_result})
            tm.assert_frame_equal(start_dataframe, expected_dataframe)

    def test_none_coercion_mixed_dtypes(self):
        start_dataframe = DataFrame(
            {
                "a": [1, 2, 3],
                "b": [1.0, 2.0, 3.0],
                "c": [datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)],
                "d": ["a", "b", "c"],
            }
        )
        start_dataframe.iloc[0] = None

        exp = DataFrame(
            {
                "a": [np.nan, 2, 3],
                "b": [np.nan, 2.0, 3.0],
                "c": [NaT, datetime(2000, 1, 2), datetime(2000, 1, 3)],
                "d": [None, "b", "c"],
            }
        )
        tm.assert_frame_equal(start_dataframe, exp)


def test_validate_indices_ok():
    indices = np.asarray([0, 1])
    validate_indices(indices, 2)
    validate_indices(indices[:0], 0)
    validate_indices(np.array([-1, -1]), 0)


def test_validate_indices_low():
    indices = np.asarray([0, -2])
    with pytest.raises(ValueError, match="'indices' contains"):
        validate_indices(indices, 2)


def test_validate_indices_high():
    indices = np.asarray([0, 1, 2])
    with pytest.raises(IndexError, match="indices are out"):
        validate_indices(indices, 2)


def test_validate_indices_empty():
    with pytest.raises(IndexError, match="indices are out"):
        validate_indices(np.array([0, 1]), 0)


def test_extension_array_cross_section():
    # A cross-section of a homogeneous EA should be an EA
    df = pd.DataFrame(
        {
            "A": pd.core.arrays.integer_array([1, 2]),
            "B": pd.core.arrays.integer_array([3, 4]),
        },
        index=["a", "b"],
    )
    expected = pd.Series(
        pd.core.arrays.integer_array([1, 3]), index=["A", "B"], name="a"
    )
    result = df.loc["a"]
    tm.assert_series_equal(result, expected)

    result = df.iloc[0]
    tm.assert_series_equal(result, expected)


def test_extension_array_cross_section_converts():
    df = pd.DataFrame(
        {"A": pd.core.arrays.integer_array([1, 2]), "B": np.array([1, 2])},
        index=["a", "b"],
    )
    result = df.loc["a"]
    expected = pd.Series([1, 1], dtype=object, index=["A", "B"], name="a")
    tm.assert_series_equal(result, expected)

    result = df.iloc[0]
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "idxr, error, error_message",
    [
        (lambda x: x, AbstractMethodError, None),
        (
            lambda x: x.loc,
            AttributeError,
            "type object 'NDFrame' has no attribute '_AXIS_NAMES'",
        ),
        (
            lambda x: x.iloc,
            AttributeError,
            "type object 'NDFrame' has no attribute '_AXIS_NAMES'",
        ),
    ],
)
def test_ndframe_indexing_raises(idxr, error, error_message):
    # GH 25567
    frame = NDFrame(np.random.randint(5, size=(2, 2, 2)))
    with pytest.raises(error, match=error_message):
        idxr(frame)[0]


def test_readonly_indices():
    # GH#17192 iloc with read-only array raising TypeError
    df = pd.DataFrame({"data": np.ones(100, dtype="float64")})
    indices = np.array([1, 3, 6])
    indices.flags.writeable = False

    result = df.iloc[indices]
    expected = df.loc[[1, 3, 6]]
    tm.assert_frame_equal(result, expected)

    result = df["data"].iloc[indices]
    expected = df["data"].loc[[1, 3, 6]]
    tm.assert_series_equal(result, expected)


def test_1tuple_without_multiindex():
    ser = pd.Series(range(5))
    key = (slice(3),)

    result = ser[key]
    expected = ser[key[0]]
    tm.assert_series_equal(result, expected)


def test_duplicate_index_mistyped_key_raises_keyerror():
    # GH#29189 float_index.get_loc(None) should raise KeyError, not TypeError
    ser = pd.Series([2, 5, 6, 8], index=[2.0, 4.0, 4.0, 5.0])
    with pytest.raises(KeyError):
        ser[None]

    with pytest.raises(KeyError):
        ser.index.get_loc(None)

    with pytest.raises(KeyError):
        ser.index._engine.get_loc(None)
