import numpy as np
import pytest

import pandas as pd
from pandas import CategoricalIndex, Index
import pandas._testing as tm


class TestTake:
    def test_take_fill_value(self):
        # GH 12631

        # numeric category
        idx = pd.CategoricalIndex([1, 2, 3], name="xxx")
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.CategoricalIndex([2, 1, 3], name="xxx")
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.CategoricalIndex([2, 1, np.nan], categories=[1, 2, 3], name="xxx")
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = pd.CategoricalIndex([2, 1, 3], name="xxx")
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # object category
        idx = pd.CategoricalIndex(
            list("CBA"), categories=list("ABC"), ordered=True, name="xxx"
        )
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.CategoricalIndex(
            list("BCA"), categories=list("ABC"), ordered=True, name="xxx"
        )
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.CategoricalIndex(
            ["B", "C", np.nan], categories=list("ABC"), ordered=True, name="xxx"
        )
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = pd.CategoricalIndex(
            list("BCA"), categories=list("ABC"), ordered=True, name="xxx"
        )
        tm.assert_index_equal(result, expected)
        tm.assert_categorical_equal(result.values, expected.values)

        msg = (
            "When allow_fill=True and fill_value is not None, "
            "all indices must be >= -1"
        )
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        msg = "index -5 is out of bounds for (axis 0 with )?size 3"
        with pytest.raises(IndexError, match=msg):
            idx.take(np.array([1, -5]))

    def test_take_fill_value_datetime(self):

        # datetime category
        idx = pd.DatetimeIndex(["2011-01-01", "2011-02-01", "2011-03-01"], name="xxx")
        idx = pd.CategoricalIndex(idx)
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.DatetimeIndex(
            ["2011-02-01", "2011-01-01", "2011-03-01"], name="xxx"
        )
        expected = pd.CategoricalIndex(expected)
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.DatetimeIndex(["2011-02-01", "2011-01-01", "NaT"], name="xxx")
        exp_cats = pd.DatetimeIndex(["2011-01-01", "2011-02-01", "2011-03-01"])
        expected = pd.CategoricalIndex(expected, categories=exp_cats)
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False, fill_value=True)
        expected = pd.DatetimeIndex(
            ["2011-02-01", "2011-01-01", "2011-03-01"], name="xxx"
        )
        expected = pd.CategoricalIndex(expected)
        tm.assert_index_equal(result, expected)

        msg = (
            "When allow_fill=True and fill_value is not None, "
            "all indices must be >= -1"
        )
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        msg = "index -5 is out of bounds for (axis 0 with )?size 3"
        with pytest.raises(IndexError, match=msg):
            idx.take(np.array([1, -5]))

    def test_take_invalid_kwargs(self):
        idx = pd.CategoricalIndex([1, 2, 3], name="foo")
        indices = [1, 0, -1]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        with pytest.raises(TypeError, match=msg):
            idx.take(indices, foo=2)

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            idx.take(indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            idx.take(indices, mode="clip")


class TestGetLoc:
    def test_get_loc(self):
        # GH 12531
        cidx1 = CategoricalIndex(list("abcde"), categories=list("edabc"))
        idx1 = Index(list("abcde"))
        assert cidx1.get_loc("a") == idx1.get_loc("a")
        assert cidx1.get_loc("e") == idx1.get_loc("e")

        for i in [cidx1, idx1]:
            with pytest.raises(KeyError, match="'NOT-EXIST'"):
                i.get_loc("NOT-EXIST")

        # non-unique
        cidx2 = CategoricalIndex(list("aacded"), categories=list("edabc"))
        idx2 = Index(list("aacded"))

        # results in bool array
        res = cidx2.get_loc("d")
        tm.assert_numpy_array_equal(res, idx2.get_loc("d"))
        tm.assert_numpy_array_equal(
            res, np.array([False, False, False, True, False, True])
        )
        # unique element results in scalar
        res = cidx2.get_loc("e")
        assert res == idx2.get_loc("e")
        assert res == 4

        for i in [cidx2, idx2]:
            with pytest.raises(KeyError, match="'NOT-EXIST'"):
                i.get_loc("NOT-EXIST")

        # non-unique, sliceable
        cidx3 = CategoricalIndex(list("aabbb"), categories=list("abc"))
        idx3 = Index(list("aabbb"))

        # results in slice
        res = cidx3.get_loc("a")
        assert res == idx3.get_loc("a")
        assert res == slice(0, 2, None)

        res = cidx3.get_loc("b")
        assert res == idx3.get_loc("b")
        assert res == slice(2, 5, None)

        for i in [cidx3, idx3]:
            with pytest.raises(KeyError, match="'c'"):
                i.get_loc("c")

    def test_get_loc_unique(self):
        cidx = pd.CategoricalIndex(list("abc"))
        result = cidx.get_loc("b")
        assert result == 1

    def test_get_loc_monotonic_nonunique(self):
        cidx = pd.CategoricalIndex(list("abbc"))
        result = cidx.get_loc("b")
        expected = slice(1, 3, None)
        assert result == expected

    def test_get_loc_nonmonotonic_nonunique(self):
        cidx = pd.CategoricalIndex(list("abcb"))
        result = cidx.get_loc("b")
        expected = np.array([False, True, False, True], dtype=bool)
        tm.assert_numpy_array_equal(result, expected)


class TestGetIndexer:
    def test_get_indexer_base(self):
        # Determined by cat ordering.
        idx = CategoricalIndex(list("cab"), categories=list("cab"))
        expected = np.arange(len(idx), dtype=np.intp)

        actual = idx.get_indexer(idx)
        tm.assert_numpy_array_equal(expected, actual)

        with pytest.raises(ValueError, match="Invalid fill method"):
            idx.get_indexer(idx, method="invalid")

    def test_get_indexer_non_unique(self):
        np.random.seed(123456789)

        ci = CategoricalIndex(list("aabbca"), categories=list("cab"), ordered=False)
        oidx = Index(np.array(ci))

        for n in [1, 2, 5, len(ci)]:
            finder = oidx[np.random.randint(0, len(ci), size=n)]
            expected = oidx.get_indexer_non_unique(finder)[0]

            actual = ci.get_indexer(finder)
            tm.assert_numpy_array_equal(expected, actual)

        # see gh-17323
        #
        # Even when indexer is equal to the
        # members in the index, we should
        # respect duplicates instead of taking
        # the fast-track path.
        for finder in [list("aabbca"), list("aababca")]:
            expected = oidx.get_indexer_non_unique(finder)[0]

            actual = ci.get_indexer(finder)
            tm.assert_numpy_array_equal(expected, actual)

    def test_get_indexer(self):

        idx1 = CategoricalIndex(list("aabcde"), categories=list("edabc"))
        idx2 = CategoricalIndex(list("abf"))

        for indexer in [idx2, list("abf"), Index(list("abf"))]:
            r1 = idx1.get_indexer(idx2)
            tm.assert_almost_equal(r1, np.array([0, 1, 2, -1], dtype=np.intp))

        msg = (
            "method='pad' and method='backfill' not implemented yet for "
            "CategoricalIndex"
        )
        with pytest.raises(NotImplementedError, match=msg):
            idx2.get_indexer(idx1, method="pad")
        with pytest.raises(NotImplementedError, match=msg):
            idx2.get_indexer(idx1, method="backfill")

        msg = "method='nearest' not implemented yet for CategoricalIndex"
        with pytest.raises(NotImplementedError, match=msg):
            idx2.get_indexer(idx1, method="nearest")
