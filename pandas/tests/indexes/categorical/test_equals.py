import numpy as np
import pytest

from pandas import (
    Categorical,
    CategoricalIndex,
    Index,
    MultiIndex,
)


class TestEquals:
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

        msg = "Categoricals can only be compared if 'categories' are the same"
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

        ci = CategoricalIndex([*list("aabca"), np.nan], categories=["c", "a", "b"])
        assert not ci.equals(list("aabca"))
        assert not ci.equals(CategoricalIndex(list("aabca")))
        assert ci.equals(ci.copy())

        ci = CategoricalIndex([*list("aabca"), np.nan], categories=["c", "a", "b"])
        assert not ci.equals([*list("aabca"), np.nan])
        assert ci.equals(CategoricalIndex([*list("aabca"), np.nan]))
        assert not ci.equals(CategoricalIndex([*list("aabca"), np.nan], ordered=True))
        assert ci.equals(ci.copy())

    def test_equals_categorical_unordered(self):
        # https://github.com/pandas-dev/pandas/issues/16603
        a = CategoricalIndex(["A"], categories=["A", "B"])
        b = CategoricalIndex(["A"], categories=["B", "A"])
        c = CategoricalIndex([None], categories=["B", "A"])
        assert a.equals(b)
        assert not a.equals(c)
        assert not b.equals(c)

    def test_equals_non_category(self):
        # GH#37667 Case where other contains a value not among ci's
        #  categories ("D") and also contains np.nan
        ci = CategoricalIndex(["A", "B", np.nan, np.nan])
        other = Index(["A", "B", "D", np.nan])

        assert not ci.equals(other)

    def test_equals_multiindex(self):
        # dont raise NotImplementedError when calling is_dtype_compat

        mi = MultiIndex.from_arrays([["A", "B", "C", "D"], range(4)])
        ci = mi.to_flat_index().astype("category")

        assert not ci.equals(mi)

    def test_equals_string_dtype(self, any_string_dtype):
        # GH#55364
        idx = CategoricalIndex(list("abc"), name="B")
        other = Index(["a", "b", "c"], name="B", dtype=any_string_dtype)
        assert idx.equals(other)


class TestCategoricalIndexEqualsCheckDtype:
    """Tests for CategoricalIndex.equals() with check_dtype parameter."""

    def test_equals_categorical_ordered_check_dtype(self):
        """Test categorical indexes with ordered=True/False."""
        idx1 = pd.CategoricalIndex(["a", "b", "c"], ordered=True)
        idx2 = pd.CategoricalIndex(["a", "b", "c"], ordered=False)

        # Different ordered flags mean different dtypes
        assert idx1.equals(idx2, check_dtype=True) is False
        # But values are the same
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_categorical_different_categories_order(self):
        """Test categorical indexes with same categories but different order."""
        idx1 = pd.CategoricalIndex(
            ["a", "b"], categories=["a", "b", "c"]
        )
        idx2 = pd.CategoricalIndex(
            ["a", "b"], categories=["c", "b", "a"]
        )

        # Different category ordering is a dtype difference
        assert idx1.equals(idx2, check_dtype=True) is False
        # Values are the same
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_categorical_vs_object(self):
        """Test categorical index vs object index with same values."""
        idx1 = pd.CategoricalIndex(["a", "b", "c"])
        idx2 = pd.Index(["a", "b", "c"], dtype=object)

        # Different types
        assert idx1.equals(idx2, check_dtype=True) is False
        # Same values
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_categorical_same_dtype(self):
        """Test categorical indexes with identical dtype."""
        idx1 = pd.CategoricalIndex(["a", "b", "c"], categories=["a", "b", "c"])
        idx2 = pd.CategoricalIndex(["a", "b", "c"], categories=["a", "b", "c"])

        assert idx1.equals(idx2, check_dtype=True) is True
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_categorical_different_values(self):
        """Even with check_dtype=False, different values should not be equal."""
        idx1 = pd.CategoricalIndex(["a", "b", "c"])
        idx2 = pd.CategoricalIndex(["a", "b", "d"])

        assert idx1.equals(idx2, check_dtype=False) is False