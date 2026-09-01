import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("c", [None, [1, 2, 3, 4, 5]])
def test_categorical_equal(c):
    c = pd.Categorical([1, 2, 3, 4], categories=c)
    tm.assert_categorical_equal(c, c)


@pytest.mark.parametrize("check_category_order", [True, False])
def test_categorical_equal_order_mismatch(check_category_order):
    c1 = pd.Categorical([1, 2, 3, 4], categories=[1, 2, 3, 4])
    c2 = pd.Categorical([1, 2, 3, 4], categories=[4, 3, 2, 1])
    kwargs = {"check_category_order": check_category_order}

    if check_category_order:
        msg = """Categorical\\.categories are different

Categorical\\.categories values are different \\(100\\.0 %\\)
\\[left\\]:  Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: Index\\(\\[4, 3, 2, 1\\], dtype='int64'\\)"""
        with pytest.raises(AssertionError, match=msg):
            tm.assert_categorical_equal(c1, c2, **kwargs)
    else:
        tm.assert_categorical_equal(c1, c2, **kwargs)


def test_categorical_equal_categories_mismatch():
    msg = """Categorical\\.categories are different

Categorical\\.categories values are different \\(25\\.0 %\\)
\\[left\\]:  Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: Index\\(\\[1, 2, 3, 5\\], dtype='int64'\\)"""

    c1 = pd.Categorical([1, 2, 3, 4])
    c2 = pd.Categorical([1, 2, 3, 5])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c1, c2)


def test_categorical_equal_codes_mismatch():
    categories = [1, 2, 3, 4]
    msg = """Categorical\\.codes are different

Categorical\\.codes values are different \\(50\\.0 %\\)
\\[left\\]:  \\[0, 1, 3, 2\\]
\\[right\\]: \\[0, 1, 2, 3\\]"""

    c1 = pd.Categorical([1, 2, 4, 3], categories=categories)
    c2 = pd.Categorical([1, 2, 3, 4], categories=categories)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c1, c2)


def test_categorical_equal_ordered_mismatch():
    data = [1, 2, 3, 4]
    msg = """Categorical are different

Attribute "ordered" are different
\\[left\\]:  False
\\[right\\]: True"""

    c1 = pd.Categorical(data, ordered=False)
    c2 = pd.Categorical(data, ordered=True)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c1, c2)


@pytest.mark.parametrize("obj", ["index", "foo", "pandas"])
def test_categorical_equal_object_override(obj):
    data = [1, 2, 3, 4]
    msg = f"""{obj} are different

Attribute "ordered" are different
\\[left\\]:  False
\\[right\\]: True"""

    c1 = pd.Categorical(data, ordered=False)
    c2 = pd.Categorical(data, ordered=True)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c1, c2, obj=obj)


def test_categorical_equal_order_mismatch_with_nan():
    # GH#62008 the -1 code for NA must not be treated as a positional indexer
    c1 = Categorical(["B", None, "D"], categories=["B", "D"])
    c2 = Categorical(["B", None, "D"], categories=["D", "B"])

    tm.assert_categorical_equal(c1, c2, check_category_order=False)


def test_categorical_equal_nan_vs_value():
    # GH#62008 an NA on one side must not compare equal to a value on the other
    msg = """Categorical\\.codes are different

Categorical\\.codes values are different \\(50\\.0 %\\)
\\[left\\]:  \\[0, -1\\]
\\[right\\]: \\[0, 1\\]"""

    c1 = Categorical(["B", None], categories=["B", "D"])
    c2 = Categorical(["B", "D"], categories=["D", "B"])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c1, c2, check_category_order=False)


def test_categorical_equal_large_int_no_float_cast():
    # GH#62008 guards against the rejected fix of comparing values via
    #  take(codes, allow_fill=True, fill_value=np.nan), which casts int64
    #  categories to float64 and collapses 2**53 and 2**53 + 1 into one value
    msg = """Categorical\\.codes are different

Categorical\\.codes values are different \\(66\\.66667 %\\)
\\[left\\]:  \\[0, -1, 1\\]
\\[right\\]: \\[1, -1, 0\\]"""

    categories = [2**53, 2**53 + 1]
    c1 = Categorical([2**53, None, 2**53 + 1], categories=categories)
    c2 = Categorical([2**53 + 1, None, 2**53], categories=categories)

    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c1, c2, check_category_order=False)


def test_categorical_equal_unsortable_categories():
    # GH#62008 categories that cannot be sorted fall back to comparing in the
    #  original order; NA must still not be read as a positional indexer
    categories = [1, "a"]
    c1 = Categorical([1, None, "a"], categories=categories)
    c2 = Categorical([1, None, "a"], categories=categories)

    tm.assert_categorical_equal(c1, c2, check_category_order=False)

    c3 = Categorical([1, None, "a"], categories=categories)
    c4 = Categorical(["a", None, 1], categories=categories)
    msg = """Categorical\\.codes are different

Categorical\\.codes values are different \\(66\\.66667 %\\)
\\[left\\]:  \\[0, -1, 1\\]
\\[right\\]: \\[1, -1, 0\\]"""
    with pytest.raises(AssertionError, match=msg):
        tm.assert_categorical_equal(c3, c4, check_category_order=False)
