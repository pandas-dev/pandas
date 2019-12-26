import numpy as np
import pytest

import pandas._config.config as cf

from pandas._libs import index as libindex

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import Categorical, IntervalIndex
from pandas.core.indexes.api import CategoricalIndex, Index
import pandas.util.testing as tm

from ..common import Base


class TestCategoricalIndex(Base):
    _holder = CategoricalIndex

    @pytest.fixture
    def indices(self, request):
        return tm.makeCategoricalIndex(100)

    def create_index(self, categories=None, ordered=False):
        if categories is None:
            categories = list("cab")
        return CategoricalIndex(list("aabbca"), categories=categories, ordered=ordered)

    def test_can_hold_identifiers(self):
        idx = self.create_index(categories=list("abcd"))
        key = idx[0]
        assert idx._can_hold_identifiers_and_holds_name(key) is True

    @pytest.mark.parametrize(
        "func,op_name",
        [
            (lambda idx: idx - idx, "__sub__"),
            (lambda idx: idx + idx, "__add__"),
            (lambda idx: idx - ["a", "b"], "__sub__"),
            (lambda idx: idx + ["a", "b"], "__add__"),
            (lambda idx: ["a", "b"] - idx, "__rsub__"),
            (lambda idx: ["a", "b"] + idx, "__radd__"),
        ],
    )
    def test_disallow_set_ops(self, func, op_name):
        # GH 10039
        # set ops (+/-) raise TypeError
        idx = pd.Index(pd.Categorical(["a", "b"]))
        msg = f"cannot perform {op_name} with this index type: CategoricalIndex"
        with pytest.raises(TypeError, match=msg):
            func(idx)

    def test_method_delegation(self):

        ci = CategoricalIndex(list("aabbca"), categories=list("cabdef"))
        result = ci.set_categories(list("cab"))
        tm.assert_index_equal(
            result, CategoricalIndex(list("aabbca"), categories=list("cab"))
        )

        ci = CategoricalIndex(list("aabbca"), categories=list("cab"))
        result = ci.rename_categories(list("efg"))
        tm.assert_index_equal(
            result, CategoricalIndex(list("ffggef"), categories=list("efg"))
        )

        # GH18862 (let rename_categories take callables)
        result = ci.rename_categories(lambda x: x.upper())
        tm.assert_index_equal(
            result, CategoricalIndex(list("AABBCA"), categories=list("CAB"))
        )

        ci = CategoricalIndex(list("aabbca"), categories=list("cab"))
        result = ci.add_categories(["d"])
        tm.assert_index_equal(
            result, CategoricalIndex(list("aabbca"), categories=list("cabd"))
        )

        ci = CategoricalIndex(list("aabbca"), categories=list("cab"))
        result = ci.remove_categories(["c"])
        tm.assert_index_equal(
            result,
            CategoricalIndex(list("aabb") + [np.nan] + ["a"], categories=list("ab")),
        )

        ci = CategoricalIndex(list("aabbca"), categories=list("cabdef"))
        result = ci.as_unordered()
        tm.assert_index_equal(result, ci)

        ci = CategoricalIndex(list("aabbca"), categories=list("cabdef"))
        result = ci.as_ordered()
        tm.assert_index_equal(
            result,
            CategoricalIndex(list("aabbca"), categories=list("cabdef"), ordered=True),
        )

        # invalid
        msg = "cannot use inplace with CategoricalIndex"
        with pytest.raises(ValueError, match=msg):
            ci.set_categories(list("cab"), inplace=True)

    def test_contains(self):

        ci = self.create_index(categories=list("cabdef"))

        assert "a" in ci
        assert "z" not in ci
        assert "e" not in ci
        assert np.nan not in ci

        # assert codes NOT in index
        assert 0 not in ci
        assert 1 not in ci

        ci = CategoricalIndex(list("aabbca") + [np.nan], categories=list("cabdef"))
        assert np.nan in ci

    @pytest.mark.parametrize(
        "item, expected",
        [
            (pd.Interval(0, 1), True),
            (1.5, True),
            (pd.Interval(0.5, 1.5), False),
            ("a", False),
            (pd.Timestamp(1), False),
            (pd.Timedelta(1), False),
        ],
        ids=str,
    )
    def test_contains_interval(self, item, expected):
        # GH 23705
        ci = CategoricalIndex(IntervalIndex.from_breaks(range(3)))
        result = item in ci
        assert result is expected

    def test_contains_list(self):
        # GH#21729
        idx = pd.CategoricalIndex([1, 2, 3])

        assert "a" not in idx

        with pytest.raises(TypeError, match="unhashable type"):
            ["a"] in idx

        with pytest.raises(TypeError, match="unhashable type"):
            ["a", "b"] in idx

    def test_map(self):
        ci = pd.CategoricalIndex(list("ABABC"), categories=list("CBA"), ordered=True)
        result = ci.map(lambda x: x.lower())
        exp = pd.CategoricalIndex(list("ababc"), categories=list("cba"), ordered=True)
        tm.assert_index_equal(result, exp)

        ci = pd.CategoricalIndex(
            list("ABABC"), categories=list("BAC"), ordered=False, name="XXX"
        )
        result = ci.map(lambda x: x.lower())
        exp = pd.CategoricalIndex(
            list("ababc"), categories=list("bac"), ordered=False, name="XXX"
        )
        tm.assert_index_equal(result, exp)

        # GH 12766: Return an index not an array
        tm.assert_index_equal(
            ci.map(lambda x: 1), Index(np.array([1] * 5, dtype=np.int64), name="XXX")
        )

        # change categories dtype
        ci = pd.CategoricalIndex(list("ABABC"), categories=list("BAC"), ordered=False)

        def f(x):
            return {"A": 10, "B": 20, "C": 30}.get(x)

        result = ci.map(f)
        exp = pd.CategoricalIndex(
            [10, 20, 10, 20, 30], categories=[20, 10, 30], ordered=False
        )
        tm.assert_index_equal(result, exp)

        result = ci.map(pd.Series([10, 20, 30], index=["A", "B", "C"]))
        tm.assert_index_equal(result, exp)

        result = ci.map({"A": 10, "B": 20, "C": 30})
        tm.assert_index_equal(result, exp)

    def test_map_with_categorical_series(self):
        # GH 12756
        a = pd.Index([1, 2, 3, 4])
        b = pd.Series(["even", "odd", "even", "odd"], dtype="category")
        c = pd.Series(["even", "odd", "even", "odd"])

        exp = CategoricalIndex(["odd", "even", "odd", np.nan])
        tm.assert_index_equal(a.map(b), exp)
        exp = pd.Index(["odd", "even", "odd", np.nan])
        tm.assert_index_equal(a.map(c), exp)

    @pytest.mark.parametrize(
        ("data", "f"),
        (
            ([1, 1, np.nan], pd.isna),
            ([1, 2, np.nan], pd.isna),
            ([1, 1, np.nan], {1: False}),
            ([1, 2, np.nan], {1: False, 2: False}),
            ([1, 1, np.nan], pd.Series([False, False])),
            ([1, 2, np.nan], pd.Series([False, False, False])),
        ),
    )
    def test_map_with_nan(self, data, f):  # GH 24241
        values = pd.Categorical(data)
        result = values.map(f)
        if data[1] == 1:
            expected = pd.Categorical([False, False, np.nan])
            tm.assert_categorical_equal(result, expected)
        else:
            expected = pd.Index([False, False, np.nan])
            tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("klass", [list, tuple, np.array, pd.Series])
    def test_where(self, klass):
        i = self.create_index()
        cond = [True] * len(i)
        expected = i
        result = i.where(klass(cond))
        tm.assert_index_equal(result, expected)

        cond = [False] + [True] * (len(i) - 1)
        expected = CategoricalIndex([np.nan] + i[1:].tolist(), categories=i.categories)
        result = i.where(klass(cond))
        tm.assert_index_equal(result, expected)

    def test_append(self):

        ci = self.create_index()
        categories = ci.categories

        # append cats with the same categories
        result = ci[:3].append(ci[3:])
        tm.assert_index_equal(result, ci, exact=True)

        foos = [ci[:1], ci[1:3], ci[3:]]
        result = foos[0].append(foos[1:])
        tm.assert_index_equal(result, ci, exact=True)

        # empty
        result = ci.append([])
        tm.assert_index_equal(result, ci, exact=True)

        # appending with different categories or reordered is not ok
        msg = "all inputs must be Index"
        with pytest.raises(TypeError, match=msg):
            ci.append(ci.values.set_categories(list("abcd")))
        with pytest.raises(TypeError, match=msg):
            ci.append(ci.values.reorder_categories(list("abc")))

        # with objects
        result = ci.append(Index(["c", "a"]))
        expected = CategoricalIndex(list("aabbcaca"), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # invalid objects
        msg = "cannot append a non-category item to a CategoricalIndex"
        with pytest.raises(TypeError, match=msg):
            ci.append(Index(["a", "d"]))

        # GH14298 - if base object is not categorical -> coerce to object
        result = Index(["c", "a"]).append(ci)
        expected = Index(list("caaabbca"))
        tm.assert_index_equal(result, expected, exact=True)

    def test_append_to_another(self):
        # hits Index._concat_same_dtype
        fst = Index(["a", "b"])
        snd = CategoricalIndex(["d", "e"])
        result = fst.append(snd)
        expected = Index(["a", "b", "d", "e"])
        tm.assert_index_equal(result, expected)

    def test_insert(self):

        ci = self.create_index()
        categories = ci.categories

        # test 0th element
        result = ci.insert(0, "a")
        expected = CategoricalIndex(list("aaabbca"), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # test Nth element that follows Python list behavior
        result = ci.insert(-1, "a")
        expected = CategoricalIndex(list("aabbcaa"), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # test empty
        result = CategoricalIndex(categories=categories).insert(0, "a")
        expected = CategoricalIndex(["a"], categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        # invalid
        msg = (
            "cannot insert an item into a CategoricalIndex that is not"
            " already an existing category"
        )
        with pytest.raises(TypeError, match=msg):
            ci.insert(0, "d")

        # GH 18295 (test missing)
        expected = CategoricalIndex(["a", np.nan, "a", "b", "c", "b"])
        for na in (np.nan, pd.NaT, None):
            result = CategoricalIndex(list("aabcb")).insert(1, na)
            tm.assert_index_equal(result, expected)

    def test_delete(self):

        ci = self.create_index()
        categories = ci.categories

        result = ci.delete(0)
        expected = CategoricalIndex(list("abbca"), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        result = ci.delete(-1)
        expected = CategoricalIndex(list("aabbc"), categories=categories)
        tm.assert_index_equal(result, expected, exact=True)

        with pytest.raises((IndexError, ValueError)):
            # Either depending on NumPy version
            ci.delete(10)

    def test_astype(self):

        ci = self.create_index()
        result = ci.astype(object)
        tm.assert_index_equal(result, Index(np.array(ci)))

        # this IS equal, but not the same class
        assert result.equals(ci)
        assert isinstance(result, Index)
        assert not isinstance(result, CategoricalIndex)

        # interval
        ii = IntervalIndex.from_arrays(left=[-0.001, 2.0], right=[2, 4], closed="right")

        ci = CategoricalIndex(
            Categorical.from_codes([0, 1, -1], categories=ii, ordered=True)
        )

        result = ci.astype("interval")
        expected = ii.take([0, 1, -1])
        tm.assert_index_equal(result, expected)

        result = IntervalIndex(result.values)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("name", [None, "foo"])
    @pytest.mark.parametrize("dtype_ordered", [True, False])
    @pytest.mark.parametrize("index_ordered", [True, False])
    def test_astype_category(self, name, dtype_ordered, index_ordered):
        # GH 18630
        index = self.create_index(ordered=index_ordered)
        if name:
            index = index.rename(name)

        # standard categories
        dtype = CategoricalDtype(ordered=dtype_ordered)
        result = index.astype(dtype)
        expected = CategoricalIndex(
            index.tolist(),
            name=name,
            categories=index.categories,
            ordered=dtype_ordered,
        )
        tm.assert_index_equal(result, expected)

        # non-standard categories
        dtype = CategoricalDtype(index.unique().tolist()[:-1], dtype_ordered)
        result = index.astype(dtype)
        expected = CategoricalIndex(index.tolist(), name=name, dtype=dtype)
        tm.assert_index_equal(result, expected)

        if dtype_ordered is False:
            # dtype='category' can't specify ordered, so only test once
            result = index.astype("category")
            expected = index
            tm.assert_index_equal(result, expected)

    def test_reindex_base(self):
        # Determined by cat ordering.
        idx = CategoricalIndex(list("cab"), categories=list("cab"))
        expected = np.arange(len(idx), dtype=np.intp)

        actual = idx.get_indexer(idx)
        tm.assert_numpy_array_equal(expected, actual)

        with pytest.raises(ValueError, match="Invalid fill method"):
            idx.get_indexer(idx, method="invalid")

    def test_reindexing(self):
        np.random.seed(123456789)

        ci = self.create_index()
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

    def test_reindex_dtype(self):
        c = CategoricalIndex(["a", "b", "c", "a"])
        res, indexer = c.reindex(["a", "c"])
        tm.assert_index_equal(res, Index(["a", "a", "c"]), exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

        c = CategoricalIndex(["a", "b", "c", "a"])
        res, indexer = c.reindex(Categorical(["a", "c"]))

        exp = CategoricalIndex(["a", "a", "c"], categories=["a", "c"])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

        c = CategoricalIndex(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        res, indexer = c.reindex(["a", "c"])
        exp = Index(["a", "a", "c"], dtype="object")
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

        c = CategoricalIndex(["a", "b", "c", "a"], categories=["a", "b", "c", "d"])
        res, indexer = c.reindex(Categorical(["a", "c"]))
        exp = CategoricalIndex(["a", "a", "c"], categories=["a", "c"])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 3, 2], dtype=np.intp))

    def test_reindex_duplicate_target(self):
        # See GH25459
        cat = CategoricalIndex(["a", "b", "c"], categories=["a", "b", "c", "d"])
        res, indexer = cat.reindex(["a", "c", "c"])
        exp = Index(["a", "c", "c"], dtype="object")
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 2, 2], dtype=np.intp))

        res, indexer = cat.reindex(
            CategoricalIndex(["a", "c", "c"], categories=["a", "b", "c", "d"])
        )
        exp = CategoricalIndex(["a", "c", "c"], categories=["a", "b", "c", "d"])
        tm.assert_index_equal(res, exp, exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([0, 2, 2], dtype=np.intp))

    def test_reindex_empty_index(self):
        # See GH16770
        c = CategoricalIndex([])
        res, indexer = c.reindex(["a", "b"])
        tm.assert_index_equal(res, Index(["a", "b"]), exact=True)
        tm.assert_numpy_array_equal(indexer, np.array([-1, -1], dtype=np.intp))

    @pytest.mark.parametrize(
        "data, non_lexsorted_data",
        [[[1, 2, 3], [9, 0, 1, 2, 3]], [list("abc"), list("fabcd")]],
    )
    def test_is_monotonic(self, data, non_lexsorted_data):
        c = CategoricalIndex(data)
        assert c.is_monotonic_increasing is True
        assert c.is_monotonic_decreasing is False

        c = CategoricalIndex(data, ordered=True)
        assert c.is_monotonic_increasing is True
        assert c.is_monotonic_decreasing is False

        c = CategoricalIndex(data, categories=reversed(data))
        assert c.is_monotonic_increasing is False
        assert c.is_monotonic_decreasing is True

        c = CategoricalIndex(data, categories=reversed(data), ordered=True)
        assert c.is_monotonic_increasing is False
        assert c.is_monotonic_decreasing is True

        # test when data is neither monotonic increasing nor decreasing
        reordered_data = [data[0], data[2], data[1]]
        c = CategoricalIndex(reordered_data, categories=reversed(data))
        assert c.is_monotonic_increasing is False
        assert c.is_monotonic_decreasing is False

        # non lexsorted categories
        categories = non_lexsorted_data

        c = CategoricalIndex(categories[:2], categories=categories)
        assert c.is_monotonic_increasing is True
        assert c.is_monotonic_decreasing is False

        c = CategoricalIndex(categories[1:3], categories=categories)
        assert c.is_monotonic_increasing is True
        assert c.is_monotonic_decreasing is False

    def test_has_duplicates(self):

        idx = CategoricalIndex([0, 0, 0], name="foo")
        assert idx.is_unique is False
        assert idx.has_duplicates is True

    def test_drop_duplicates(self):

        idx = CategoricalIndex([0, 0, 0], name="foo")
        expected = CategoricalIndex([0], name="foo")
        tm.assert_index_equal(idx.drop_duplicates(), expected)
        tm.assert_index_equal(idx.unique(), expected)

    def test_get_indexer(self):

        idx1 = CategoricalIndex(list("aabcde"), categories=list("edabc"))
        idx2 = CategoricalIndex(list("abf"))

        for indexer in [idx2, list("abf"), Index(list("abf"))]:
            r1 = idx1.get_indexer(idx2)
            tm.assert_almost_equal(r1, np.array([0, 1, 2, -1], dtype=np.intp))

        msg = (
            "method='pad' and method='backfill' not implemented yet for"
            " CategoricalIndex"
        )
        with pytest.raises(NotImplementedError, match=msg):
            idx2.get_indexer(idx1, method="pad")
        with pytest.raises(NotImplementedError, match=msg):
            idx2.get_indexer(idx1, method="backfill")

        msg = "method='nearest' not implemented yet for CategoricalIndex"
        with pytest.raises(NotImplementedError, match=msg):
            idx2.get_indexer(idx1, method="nearest")

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

    def test_repr_roundtrip(self):

        ci = CategoricalIndex(["a", "b"], categories=["a", "b"], ordered=True)
        str(ci)
        tm.assert_index_equal(eval(repr(ci)), ci, exact=True)

        # formatting
        str(ci)

        # long format
        # this is not reprable
        ci = CategoricalIndex(np.random.randint(0, 5, size=100))
        str(ci)

    def test_isin(self):

        ci = CategoricalIndex(list("aabca") + [np.nan], categories=["c", "a", "b"])
        tm.assert_numpy_array_equal(
            ci.isin(["c"]), np.array([False, False, False, True, False, False])
        )
        tm.assert_numpy_array_equal(
            ci.isin(["c", "a", "b"]), np.array([True] * 5 + [False])
        )
        tm.assert_numpy_array_equal(
            ci.isin(["c", "a", "b", np.nan]), np.array([True] * 6)
        )

        # mismatched categorical -> coerced to ndarray so doesn't matter
        result = ci.isin(ci.set_categories(list("abcdefghi")))
        expected = np.array([True] * 6)
        tm.assert_numpy_array_equal(result, expected)

        result = ci.isin(ci.set_categories(list("defghi")))
        expected = np.array([False] * 5 + [True])
        tm.assert_numpy_array_equal(result, expected)

    def test_identical(self):

        ci1 = CategoricalIndex(["a", "b"], categories=["a", "b"], ordered=True)
        ci2 = CategoricalIndex(["a", "b"], categories=["a", "b", "c"], ordered=True)
        assert ci1.identical(ci1)
        assert ci1.identical(ci1.copy())
        assert not ci1.identical(ci2)

    def test_ensure_copied_data(self, indices):
        # gh-12309: Check the "copy" argument of each
        # Index.__new__ is honored.
        #
        # Must be tested separately from other indexes because
        # self.values is not an ndarray.
        # GH#29918 Index.base has been removed
        # FIXME: is this test still meaningful?
        _base = lambda ar: ar if getattr(ar, "base", None) is None else ar.base

        result = CategoricalIndex(indices.values, copy=True)
        tm.assert_index_equal(indices, result)
        assert _base(indices.values) is not _base(result.values)

        result = CategoricalIndex(indices.values, copy=False)
        assert _base(indices.values) is _base(result.values)

    def test_equals_categorical(self):
        ci1 = CategoricalIndex(["a", "b"], categories=["a", "b"], ordered=True)
        ci2 = CategoricalIndex(["a", "b"], categories=["a", "b", "c"], ordered=True)

        assert ci1.equals(ci1)
        assert not ci1.equals(ci2)
        assert ci1.equals(ci1.astype(object))
        assert ci1.astype(object).equals(ci1)

        assert (ci1 == ci1).all()
        assert not (ci1 != ci1).all()
        assert not (ci1 > ci1).all()
        assert not (ci1 < ci1).all()
        assert (ci1 <= ci1).all()
        assert (ci1 >= ci1).all()

        assert not (ci1 == 1).all()
        assert (ci1 == Index(["a", "b"])).all()
        assert (ci1 == ci1.values).all()

        # invalid comparisons
        with pytest.raises(ValueError, match="Lengths must match"):
            ci1 == Index(["a", "b", "c"])

        msg = (
            "categorical index comparisons must have the same categories"
            " and ordered attributes"
            "|"
            "Categoricals can only be compared if 'categories' are the same. "
            "Categories are different lengths"
            "|"
            "Categoricals can only be compared if 'ordered' is the same"
        )
        with pytest.raises(TypeError, match=msg):
            ci1 == ci2
        with pytest.raises(TypeError, match=msg):
            ci1 == Categorical(ci1.values, ordered=False)
        with pytest.raises(TypeError, match=msg):
            ci1 == Categorical(ci1.values, categories=list("abc"))

        # tests
        # make sure that we are testing for category inclusion properly
        ci = CategoricalIndex(list("aabca"), categories=["c", "a", "b"])
        assert not ci.equals(list("aabca"))
        # Same categories, but different order
        # Unordered
        assert ci.equals(CategoricalIndex(list("aabca")))
        # Ordered
        assert not ci.equals(CategoricalIndex(list("aabca"), ordered=True))
        assert ci.equals(ci.copy())

        ci = CategoricalIndex(list("aabca") + [np.nan], categories=["c", "a", "b"])
        assert not ci.equals(list("aabca"))
        assert not ci.equals(CategoricalIndex(list("aabca")))
        assert ci.equals(ci.copy())

        ci = CategoricalIndex(list("aabca") + [np.nan], categories=["c", "a", "b"])
        assert not ci.equals(list("aabca") + [np.nan])
        assert ci.equals(CategoricalIndex(list("aabca") + [np.nan]))
        assert not ci.equals(CategoricalIndex(list("aabca") + [np.nan], ordered=True))
        assert ci.equals(ci.copy())

    def test_equals_categoridcal_unordered(self):
        # https://github.com/pandas-dev/pandas/issues/16603
        a = pd.CategoricalIndex(["A"], categories=["A", "B"])
        b = pd.CategoricalIndex(["A"], categories=["B", "A"])
        c = pd.CategoricalIndex(["C"], categories=["B", "A"])
        assert a.equals(b)
        assert not a.equals(c)
        assert not b.equals(c)

    def test_frame_repr(self):
        df = pd.DataFrame({"A": [1, 2, 3]}, index=pd.CategoricalIndex(["a", "b", "c"]))
        result = repr(df)
        expected = "   A\na  1\nb  2\nc  3"
        assert result == expected

    def test_string_categorical_index_repr(self):
        # short
        idx = pd.CategoricalIndex(["a", "bb", "ccc"])
        expected = """CategoricalIndex(['a', 'bb', 'ccc'], categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""  # noqa
        assert repr(idx) == expected

        # multiple lines
        idx = pd.CategoricalIndex(["a", "bb", "ccc"] * 10)
        expected = """CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb',
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""  # noqa

        assert repr(idx) == expected

        # truncated
        idx = pd.CategoricalIndex(["a", "bb", "ccc"] * 100)
        expected = """CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  ...
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category', length=300)"""  # noqa

        assert repr(idx) == expected

        # larger categories
        idx = pd.CategoricalIndex(list("abcdefghijklmmo"))
        expected = """CategoricalIndex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                  'm', 'm', 'o'],
                 categories=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', ...], ordered=False, dtype='category')"""  # noqa

        assert repr(idx) == expected

        # short
        idx = pd.CategoricalIndex(["あ", "いい", "ううう"])
        expected = """CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa
        assert repr(idx) == expected

        # multiple lines
        idx = pd.CategoricalIndex(["あ", "いい", "ううう"] * 10)
        expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa

        assert repr(idx) == expected

        # truncated
        idx = pd.CategoricalIndex(["あ", "いい", "ううう"] * 100)
        expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""  # noqa

        assert repr(idx) == expected

        # larger categories
        idx = pd.CategoricalIndex(list("あいうえおかきくけこさしすせそ"))
        expected = """CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ', 'さ', 'し',
                  'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', ...], ordered=False, dtype='category')"""  # noqa

        assert repr(idx) == expected

        # Emable Unicode option -----------------------------------------
        with cf.option_context("display.unicode.east_asian_width", True):

            # short
            idx = pd.CategoricalIndex(["あ", "いい", "ううう"])
            expected = """CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa
            assert repr(idx) == expected

            # multiple lines
            idx = pd.CategoricalIndex(["あ", "いい", "ううう"] * 10)
            expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa

            assert repr(idx) == expected

            # truncated
            idx = pd.CategoricalIndex(["あ", "いい", "ううう"] * 100)
            expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""  # noqa

            assert repr(idx) == expected

            # larger categories
            idx = pd.CategoricalIndex(list("あいうえおかきくけこさしすせそ"))
            expected = """CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ',
                  'さ', 'し', 'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', ...], ordered=False, dtype='category')"""  # noqa

            assert repr(idx) == expected

    def test_fillna_categorical(self):
        # GH 11343
        idx = CategoricalIndex([1.0, np.nan, 3.0, 1.0], name="x")
        # fill by value in categories
        exp = CategoricalIndex([1.0, 1.0, 3.0, 1.0], name="x")
        tm.assert_index_equal(idx.fillna(1.0), exp)

        # fill by value not in categories raises ValueError
        msg = "fill value must be in categories"
        with pytest.raises(ValueError, match=msg):
            idx.fillna(2.0)

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

        with pytest.raises(IndexError):
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

        with pytest.raises(IndexError):
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

    @pytest.mark.parametrize(
        "dtype, engine_type",
        [
            (np.int8, libindex.Int8Engine),
            (np.int16, libindex.Int16Engine),
            (np.int32, libindex.Int32Engine),
            (np.int64, libindex.Int64Engine),
        ],
    )
    def test_engine_type(self, dtype, engine_type):
        if dtype != np.int64:
            # num. of uniques required to push CategoricalIndex.codes to a
            # dtype (128 categories required for .codes dtype to be int16 etc.)
            num_uniques = {np.int8: 1, np.int16: 128, np.int32: 32768}[dtype]
            ci = pd.CategoricalIndex(range(num_uniques))
        else:
            # having 2**32 - 2**31 categories would be very memory-intensive,
            # so we cheat a bit with the dtype
            ci = pd.CategoricalIndex(range(32768))  # == 2**16 - 2**(16 - 1)
            ci.values._codes = ci.values._codes.astype("int64")
        assert np.issubdtype(ci.codes.dtype, dtype)
        assert isinstance(ci._engine, engine_type)
