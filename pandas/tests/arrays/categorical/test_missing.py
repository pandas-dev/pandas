import collections

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas import (
    Categorical,
    Index,
    Series,
    isna,
)
import pandas._testing as tm


class TestCategoricalMissing:
    def test_isna(self):
        exp = np.array([False, False, True])
        cat = Categorical(["a", "b", np.nan])
        res = cat.isna()

        tm.assert_numpy_array_equal(res, exp)

    def test_na_flags_int_categories(self):
        # #1457

        categories = list(range(10))
        labels = np.random.default_rng(2).integers(0, 10, 20)
        labels[::5] = -1
        msg = "Constructing a Categorical with a dtype and values containing"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            cat = Categorical(labels, categories)
        repr(cat)

        tm.assert_numpy_array_equal(isna(cat), labels == -1)

    def test_nan_handling(self):
        # Nans are represented as -1 in codes
        c = Categorical(["a", "b", np.nan, "a"])
        tm.assert_index_equal(c.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(c._codes, np.array([0, 1, -1, 0], dtype=np.int8))
        c[1] = np.nan
        tm.assert_index_equal(c.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(c._codes, np.array([0, -1, -1, 0], dtype=np.int8))

        # Adding nan to categories should make assigned nan point to the
        # category!
        c = Categorical(["a", "b", np.nan, "a"])
        tm.assert_index_equal(c.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(c._codes, np.array([0, 1, -1, 0], dtype=np.int8))

    def test_set_dtype_nans(self):
        c = Categorical(["a", "b", np.nan])
        result = c._set_dtype(CategoricalDtype(["a", "c"]), copy=True)
        tm.assert_numpy_array_equal(result.codes, np.array([0, -1, -1], dtype="int8"))

    def test_set_item_nan(self):
        cat = Categorical([1, 2, 3])
        cat[1] = np.nan

        exp = Categorical([1, np.nan, 3], categories=[1, 2, 3])
        tm.assert_categorical_equal(cat, exp)

    @pytest.mark.parametrize("named", [True, False])
    def test_fillna_iterable_category(self, named):
        # https://github.com/pandas-dev/pandas/issues/21097
        if named:
            Point = collections.namedtuple("Point", "x y")
        else:
            Point = lambda *args: args  # tuple
        cat = Categorical(np.array([Point(0, 0), Point(0, 1), None], dtype=object))
        result = cat.fillna(Point(0, 0))
        expected = Categorical([Point(0, 0), Point(0, 1), Point(0, 0)])

        tm.assert_categorical_equal(result, expected)

        # Case where the Point is not among our categories; we want ValueError,
        #  not NotImplementedError GH#41914
        cat = Categorical(np.array([Point(1, 0), Point(0, 1), None], dtype=object))
        msg = "Cannot setitem on a Categorical with a new category"
        with pytest.raises(TypeError, match=msg):
            cat.fillna(Point(0, 0))

    def test_fillna_array(self):
        # accept Categorical or ndarray value if it holds appropriate values
        cat = Categorical(["A", "B", "C", None, None])

        other = cat.fillna("C")
        result = cat.fillna(other)
        tm.assert_categorical_equal(result, other)
        assert isna(cat[-1])  # didn't modify original inplace

        other = np.array(["A", "B", "C", "B", "A"])
        result = cat.fillna(other)
        expected = Categorical(["A", "B", "C", "B", "A"], dtype=cat.dtype)
        tm.assert_categorical_equal(result, expected)
        assert isna(cat[-1])  # didn't modify original inplace

    @pytest.mark.parametrize(
        "a1, a2, categories",
        [
            (["a", "b", "c"], [np.nan, "a", "b"], ["a", "b", "c"]),
            ([1, 2, 3], [np.nan, 1, 2], [1, 2, 3]),
        ],
    )
    def test_compare_categorical_with_missing(self, a1, a2, categories):
        # GH 28384
        cat_type = CategoricalDtype(categories)

        # !=
        result = Series(a1, dtype=cat_type) != Series(a2, dtype=cat_type)
        expected = Series(a1) != Series(a2)
        tm.assert_series_equal(result, expected)

        # ==
        result = Series(a1, dtype=cat_type) == Series(a2, dtype=cat_type)
        expected = Series(a1) == Series(a2)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "na_value, dtype",
        [
            (pd.NaT, "datetime64[s]"),
            (None, "float64"),
            (np.nan, "float64"),
            (pd.NA, "float64"),
        ],
    )
    def test_categorical_only_missing_values_no_cast(self, na_value, dtype):
        # GH#44900
        result = Categorical([na_value, na_value])
        tm.assert_index_equal(result.categories, Index([], dtype=dtype))
