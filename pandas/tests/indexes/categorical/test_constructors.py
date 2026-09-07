import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestCategoricalIndexConstructors:
    def test_construction_disallows_scalar(self):
        msg = "must be called with a collection of some kind"
        with pytest.raises(TypeError, match=msg):
            pd.CategoricalIndex(data=1, categories=list("abcd"), ordered=False)
        with pytest.raises(TypeError, match=msg):
            pd.CategoricalIndex(categories=list("abcd"), ordered=False)

    def test_construction(self):
        ci = pd.CategoricalIndex(list("aabbca"), categories=list("abcd"), ordered=False)
        categories = ci.categories

        result = pd.Index(ci)
        tm.assert_index_equal(result, ci, exact=True)
        assert not result.ordered

        result = pd.Index(ci.values)
        tm.assert_index_equal(result, ci, exact=True)
        assert not result.ordered

        # empty
        result = pd.CategoricalIndex([], categories=categories)
        tm.assert_index_equal(result.categories, pd.Index(categories))
        tm.assert_numpy_array_equal(result.codes, np.array([], dtype="int8"))
        assert not result.ordered

        # passing categories
        result = pd.CategoricalIndex(list("aabbca"), categories=categories)
        tm.assert_index_equal(result.categories, pd.Index(categories))
        tm.assert_numpy_array_equal(
            result.codes, np.array([0, 0, 1, 1, 2, 0], dtype="int8")
        )

        c = pd.Categorical(list("aabbca"))
        result = pd.CategoricalIndex(c)
        tm.assert_index_equal(result.categories, pd.Index(list("abc")))
        tm.assert_numpy_array_equal(
            result.codes, np.array([0, 0, 1, 1, 2, 0], dtype="int8")
        )
        assert not result.ordered

        result = pd.CategoricalIndex(c, categories=categories)
        tm.assert_index_equal(result.categories, pd.Index(categories))
        tm.assert_numpy_array_equal(
            result.codes, np.array([0, 0, 1, 1, 2, 0], dtype="int8")
        )
        assert not result.ordered

        ci = pd.CategoricalIndex(c, categories=list("abcd"))
        result = pd.CategoricalIndex(ci)
        tm.assert_index_equal(result.categories, pd.Index(categories))
        tm.assert_numpy_array_equal(
            result.codes, np.array([0, 0, 1, 1, 2, 0], dtype="int8")
        )
        assert not result.ordered

        msg = "Constructing a Categorical with a dtype and values containing"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.CategoricalIndex(ci, categories=list("ab"))
        tm.assert_index_equal(result.categories, pd.Index(list("ab")))
        tm.assert_numpy_array_equal(
            result.codes, np.array([0, 0, 1, 1, -1, 0], dtype="int8")
        )
        assert not result.ordered

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.CategoricalIndex(ci, categories=list("ab"), ordered=True)
        tm.assert_index_equal(result.categories, pd.Index(list("ab")))
        tm.assert_numpy_array_equal(
            result.codes, np.array([0, 0, 1, 1, -1, 0], dtype="int8")
        )
        assert result.ordered

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = pd.CategoricalIndex(ci, categories=list("ab"), ordered=True)
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            expected = pd.CategoricalIndex(
                ci, categories=list("ab"), ordered=True, dtype="category"
            )
        tm.assert_index_equal(result, expected, exact=True)

        # turn me to an Index
        result = pd.Index(np.array(ci))
        assert isinstance(result, pd.Index)
        assert not isinstance(result, pd.CategoricalIndex)

    def test_construction_with_dtype(self):
        # specify dtype
        ci = pd.CategoricalIndex(list("aabbca"), categories=list("abc"), ordered=False)

        result = pd.Index(np.array(ci), dtype="category")
        tm.assert_index_equal(result, ci, exact=True)

        result = pd.Index(np.array(ci).tolist(), dtype="category")
        tm.assert_index_equal(result, ci, exact=True)

        # these are generally only equal when the categories are reordered
        ci = pd.CategoricalIndex(list("aabbca"), categories=list("cab"), ordered=False)

        result = pd.Index(np.array(ci), dtype="category").reorder_categories(
            ci.categories
        )
        tm.assert_index_equal(result, ci, exact=True)

        # make sure indexes are handled
        idx = pd.Index(range(3))
        expected = pd.CategoricalIndex([0, 1, 2], categories=idx, ordered=True)
        result = pd.CategoricalIndex(idx, categories=idx, ordered=True)
        tm.assert_index_equal(result, expected, exact=True)

    def test_construction_empty_with_bool_categories(self):
        # see GH#22702
        cat = pd.CategoricalIndex([], categories=[True, False])
        categories = sorted(cat.categories.tolist())
        assert categories == [False, True]

    def test_construction_with_categorical_dtype(self):
        # construction with CategoricalDtype
        # GH#18109
        data, cats, ordered = "a a b b".split(), "c b a".split(), True
        dtype = pd.CategoricalDtype(categories=cats, ordered=ordered)

        result = pd.CategoricalIndex(data, dtype=dtype)
        expected = pd.CategoricalIndex(data, categories=cats, ordered=ordered)
        tm.assert_index_equal(result, expected, exact=True)

        # GH#19032
        result = pd.Index(data, dtype=dtype)
        tm.assert_index_equal(result, expected, exact=True)

        # error when combining categories/ordered and dtype kwargs
        msg = "Cannot specify `categories` or `ordered` together with `dtype`."
        with pytest.raises(ValueError, match=msg):
            pd.CategoricalIndex(data, categories=cats, dtype=dtype)

        with pytest.raises(ValueError, match=msg):
            pd.CategoricalIndex(data, ordered=ordered, dtype=dtype)
