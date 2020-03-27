import re

import numpy as np
import pytest

from pandas.errors import PerformanceWarning

import pandas as pd
from pandas import DataFrame, Index, MultiIndex
import pandas._testing as tm


@pytest.mark.parametrize(
    "msg,labels,level",
    [
        (r"labels \[4\] not found in level", 4, "a"),
        (r"labels \[7\] not found in level", 7, "b"),
    ],
)
def test_drop_raise_exception_if_labels_not_in_level(msg, labels, level):
    # GH 8594
    mi = pd.MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
    s = pd.Series([10, 20, 30], index=mi)
    df = pd.DataFrame([10, 20, 30], index=mi)

    with pytest.raises(KeyError, match=msg):
        s.drop(labels, level=level)
    with pytest.raises(KeyError, match=msg):
        df.drop(labels, level=level)


@pytest.mark.parametrize("labels,level", [(4, "a"), (7, "b")])
def test_drop_errors_ignore(labels, level):
    # GH 8594
    mi = pd.MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]], names=["a", "b"])
    s = pd.Series([10, 20, 30], index=mi)
    df = pd.DataFrame([10, 20, 30], index=mi)

    expected_s = s.drop(labels, level=level, errors="ignore")
    tm.assert_series_equal(s, expected_s)

    expected_df = df.drop(labels, level=level, errors="ignore")
    tm.assert_frame_equal(df, expected_df)


def test_drop_with_non_unique_datetime_index_and_invalid_keys():
    # GH 30399

    # define dataframe with unique datetime index
    df = pd.DataFrame(
        np.random.randn(5, 3),
        columns=["a", "b", "c"],
        index=pd.date_range("2012", freq="H", periods=5),
    )
    # create dataframe with non-unique datetime index
    df = df.iloc[[0, 2, 2, 3]].copy()

    with pytest.raises(KeyError, match="not found in axis"):
        df.drop(["a", "b"])  # Dropping with labels not exist in the index


class TestDataFrameDrop:
    def test_drop_names(self):
        df = DataFrame(
            [[1, 2, 3], [3, 4, 5], [5, 6, 7]],
            index=["a", "b", "c"],
            columns=["d", "e", "f"],
        )
        df.index.name, df.columns.name = "first", "second"
        df_dropped_b = df.drop("b")
        df_dropped_e = df.drop("e", axis=1)
        df_inplace_b, df_inplace_e = df.copy(), df.copy()
        df_inplace_b.drop("b", inplace=True)
        df_inplace_e.drop("e", axis=1, inplace=True)
        for obj in (df_dropped_b, df_dropped_e, df_inplace_b, df_inplace_e):
            assert obj.index.name == "first"
            assert obj.columns.name == "second"
        assert list(df.columns) == ["d", "e", "f"]

        msg = r"\['g'\] not found in axis"
        with pytest.raises(KeyError, match=msg):
            df.drop(["g"])
        with pytest.raises(KeyError, match=msg):
            df.drop(["g"], 1)

        # errors = 'ignore'
        dropped = df.drop(["g"], errors="ignore")
        expected = Index(["a", "b", "c"], name="first")
        tm.assert_index_equal(dropped.index, expected)

        dropped = df.drop(["b", "g"], errors="ignore")
        expected = Index(["a", "c"], name="first")
        tm.assert_index_equal(dropped.index, expected)

        dropped = df.drop(["g"], axis=1, errors="ignore")
        expected = Index(["d", "e", "f"], name="second")
        tm.assert_index_equal(dropped.columns, expected)

        dropped = df.drop(["d", "g"], axis=1, errors="ignore")
        expected = Index(["e", "f"], name="second")
        tm.assert_index_equal(dropped.columns, expected)

        # GH 16398
        dropped = df.drop([], errors="ignore")
        expected = Index(["a", "b", "c"], name="first")
        tm.assert_index_equal(dropped.index, expected)

    def test_drop(self):
        simple = DataFrame({"A": [1, 2, 3, 4], "B": [0, 1, 2, 3]})
        tm.assert_frame_equal(simple.drop("A", axis=1), simple[["B"]])
        tm.assert_frame_equal(simple.drop(["A", "B"], axis="columns"), simple[[]])
        tm.assert_frame_equal(simple.drop([0, 1, 3], axis=0), simple.loc[[2], :])
        tm.assert_frame_equal(simple.drop([0, 3], axis="index"), simple.loc[[1, 2], :])

        with pytest.raises(KeyError, match=r"\[5\] not found in axis"):
            simple.drop(5)
        with pytest.raises(KeyError, match=r"\['C'\] not found in axis"):
            simple.drop("C", 1)
        with pytest.raises(KeyError, match=r"\[5\] not found in axis"):
            simple.drop([1, 5])
        with pytest.raises(KeyError, match=r"\['C'\] not found in axis"):
            simple.drop(["A", "C"], 1)

        # errors = 'ignore'
        tm.assert_frame_equal(simple.drop(5, errors="ignore"), simple)
        tm.assert_frame_equal(
            simple.drop([0, 5], errors="ignore"), simple.loc[[1, 2, 3], :]
        )
        tm.assert_frame_equal(simple.drop("C", axis=1, errors="ignore"), simple)
        tm.assert_frame_equal(
            simple.drop(["A", "C"], axis=1, errors="ignore"), simple[["B"]]
        )

        # non-unique - wheee!
        nu_df = DataFrame(
            list(zip(range(3), range(-3, 1), list("abc"))), columns=["a", "a", "b"]
        )
        tm.assert_frame_equal(nu_df.drop("a", axis=1), nu_df[["b"]])
        tm.assert_frame_equal(nu_df.drop("b", axis="columns"), nu_df["a"])
        tm.assert_frame_equal(nu_df.drop([]), nu_df)  # GH 16398

        nu_df = nu_df.set_index(pd.Index(["X", "Y", "X"]))
        nu_df.columns = list("abc")
        tm.assert_frame_equal(nu_df.drop("X", axis="rows"), nu_df.loc[["Y"], :])
        tm.assert_frame_equal(nu_df.drop(["X", "Y"], axis=0), nu_df.loc[[], :])

        # inplace cache issue
        # GH#5628
        df = pd.DataFrame(np.random.randn(10, 3), columns=list("abc"))
        expected = df[~(df.b > 0)]
        df.drop(labels=df[df.b > 0].index, inplace=True)
        tm.assert_frame_equal(df, expected)

    def test_drop_multiindex_not_lexsorted(self):
        # GH#11640

        # define the lexsorted version
        lexsorted_mi = MultiIndex.from_tuples(
            [("a", ""), ("b1", "c1"), ("b2", "c2")], names=["b", "c"]
        )
        lexsorted_df = DataFrame([[1, 3, 4]], columns=lexsorted_mi)
        assert lexsorted_df.columns.is_lexsorted()

        # define the non-lexsorted version
        not_lexsorted_df = DataFrame(
            columns=["a", "b", "c", "d"], data=[[1, "b1", "c1", 3], [1, "b2", "c2", 4]]
        )
        not_lexsorted_df = not_lexsorted_df.pivot_table(
            index="a", columns=["b", "c"], values="d"
        )
        not_lexsorted_df = not_lexsorted_df.reset_index()
        assert not not_lexsorted_df.columns.is_lexsorted()

        # compare the results
        tm.assert_frame_equal(lexsorted_df, not_lexsorted_df)

        expected = lexsorted_df.drop("a", axis=1)
        with tm.assert_produces_warning(PerformanceWarning):
            result = not_lexsorted_df.drop("a", axis=1)

        tm.assert_frame_equal(result, expected)

    def test_drop_api_equivalence(self):
        # equivalence of the labels/axis and index/columns API's (GH#12392)
        df = DataFrame(
            [[1, 2, 3], [3, 4, 5], [5, 6, 7]],
            index=["a", "b", "c"],
            columns=["d", "e", "f"],
        )

        res1 = df.drop("a")
        res2 = df.drop(index="a")
        tm.assert_frame_equal(res1, res2)

        res1 = df.drop("d", 1)
        res2 = df.drop(columns="d")
        tm.assert_frame_equal(res1, res2)

        res1 = df.drop(labels="e", axis=1)
        res2 = df.drop(columns="e")
        tm.assert_frame_equal(res1, res2)

        res1 = df.drop(["a"], axis=0)
        res2 = df.drop(index=["a"])
        tm.assert_frame_equal(res1, res2)

        res1 = df.drop(["a"], axis=0).drop(["d"], axis=1)
        res2 = df.drop(index=["a"], columns=["d"])
        tm.assert_frame_equal(res1, res2)

        msg = "Cannot specify both 'labels' and 'index'/'columns'"
        with pytest.raises(ValueError, match=msg):
            df.drop(labels="a", index="b")

        with pytest.raises(ValueError, match=msg):
            df.drop(labels="a", columns="b")

        msg = "Need to specify at least one of 'labels', 'index' or 'columns'"
        with pytest.raises(ValueError, match=msg):
            df.drop(axis=1)

    data = [[1, 2, 3], [1, 2, 3]]

    @pytest.mark.parametrize(
        "actual",
        [
            DataFrame(data=data, index=["a", "a"]),
            DataFrame(data=data, index=["a", "b"]),
            DataFrame(data=data, index=["a", "b"]).set_index([0, 1]),
            DataFrame(data=data, index=["a", "a"]).set_index([0, 1]),
        ],
    )
    def test_raise_on_drop_duplicate_index(self, actual):

        # GH#19186
        level = 0 if isinstance(actual.index, MultiIndex) else None
        msg = re.escape("\"['c'] not found in axis\"")
        with pytest.raises(KeyError, match=msg):
            actual.drop("c", level=level, axis=0)
        with pytest.raises(KeyError, match=msg):
            actual.T.drop("c", level=level, axis=1)
        expected_no_err = actual.drop("c", axis=0, level=level, errors="ignore")
        tm.assert_frame_equal(expected_no_err, actual)
        expected_no_err = actual.T.drop("c", axis=1, level=level, errors="ignore")
        tm.assert_frame_equal(expected_no_err.T, actual)

    @pytest.mark.parametrize("index", [[1, 2, 3], [1, 1, 2]])
    @pytest.mark.parametrize("drop_labels", [[], [1], [2]])
    def test_drop_empty_list(self, index, drop_labels):
        # GH#21494
        expected_index = [i for i in index if i not in drop_labels]
        frame = pd.DataFrame(index=index).drop(drop_labels)
        tm.assert_frame_equal(frame, pd.DataFrame(index=expected_index))

    @pytest.mark.parametrize("index", [[1, 2, 3], [1, 2, 2]])
    @pytest.mark.parametrize("drop_labels", [[1, 4], [4, 5]])
    def test_drop_non_empty_list(self, index, drop_labels):
        # GH# 21494
        with pytest.raises(KeyError, match="not found in axis"):
            pd.DataFrame(index=index).drop(drop_labels)
