import numpy as np
import pytest

from pandas import Categorical, Index, MultiIndex, NaT
import pandas.util.testing as tm


def test_index_equal_levels_mismatch():
    msg = """Index are different

Index levels are different
\\[left\\]:  1, Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: 2, MultiIndex\\(\\[\\('A', 1\\),
            \\('A', 2\\),
            \\('B', 3\\),
            \\('B', 4\\)\\],
           \\)"""

    idx1 = Index([1, 2, 3])
    idx2 = MultiIndex.from_tuples([("A", 1), ("A", 2), ("B", 3), ("B", 4)])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_index_equal(idx1, idx2, exact=False)


def test_index_equal_values_mismatch(check_exact):
    msg = """MultiIndex level \\[1\\] are different

MultiIndex level \\[1\\] values are different \\(25\\.0 %\\)
\\[left\\]:  Int64Index\\(\\[2, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""

    idx1 = MultiIndex.from_tuples([("A", 2), ("A", 2), ("B", 3), ("B", 4)])
    idx2 = MultiIndex.from_tuples([("A", 1), ("A", 2), ("B", 3), ("B", 4)])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_index_equal(idx1, idx2, check_exact=check_exact)


def test_index_equal_length_mismatch(check_exact):
    msg = """Index are different

Index length are different
\\[left\\]:  3, Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: 4, Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""

    idx1 = Index([1, 2, 3])
    idx2 = Index([1, 2, 3, 4])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_index_equal(idx1, idx2, check_exact=check_exact)


def test_index_equal_class_mismatch(check_exact):
    msg = """Index are different

Index classes are different
\\[left\\]:  Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: Float64Index\\(\\[1\\.0, 2\\.0, 3\\.0\\], dtype='float64'\\)"""

    idx1 = Index([1, 2, 3])
    idx2 = Index([1, 2, 3.0])

    with pytest.raises(AssertionError, match=msg):
        tm.assert_index_equal(idx1, idx2, exact=True, check_exact=check_exact)


def test_index_equal_values_close(check_exact):
    idx1 = Index([1, 2, 3.0])
    idx2 = Index([1, 2, 3.0000000001])

    if check_exact:
        msg = """Index are different

Index values are different \\(33\\.33333 %\\)
\\[left\\]:  Float64Index\\(\\[1.0, 2.0, 3.0], dtype='float64'\\)
\\[right\\]: Float64Index\\(\\[1.0, 2.0, 3.0000000001\\], dtype='float64'\\)"""

        with pytest.raises(AssertionError, match=msg):
            tm.assert_index_equal(idx1, idx2, check_exact=check_exact)
    else:
        tm.assert_index_equal(idx1, idx2, check_exact=check_exact)


def test_index_equal_values_less_close(check_exact, check_less_precise):
    idx1 = Index([1, 2, 3.0])
    idx2 = Index([1, 2, 3.0001])
    kwargs = dict(check_exact=check_exact, check_less_precise=check_less_precise)

    if check_exact or not check_less_precise:
        msg = """Index are different

Index values are different \\(33\\.33333 %\\)
\\[left\\]:  Float64Index\\(\\[1.0, 2.0, 3.0], dtype='float64'\\)
\\[right\\]: Float64Index\\(\\[1.0, 2.0, 3.0001\\], dtype='float64'\\)"""

        with pytest.raises(AssertionError, match=msg):
            tm.assert_index_equal(idx1, idx2, **kwargs)
    else:
        tm.assert_index_equal(idx1, idx2, **kwargs)


def test_index_equal_values_too_far(check_exact, check_less_precise):
    idx1 = Index([1, 2, 3])
    idx2 = Index([1, 2, 4])
    kwargs = dict(check_exact=check_exact, check_less_precise=check_less_precise)

    msg = """Index are different

Index values are different \\(33\\.33333 %\\)
\\[left\\]:  Int64Index\\(\\[1, 2, 3\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 4\\], dtype='int64'\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_index_equal(idx1, idx2, **kwargs)


def test_index_equal_level_values_mismatch(check_exact, check_less_precise):
    idx1 = MultiIndex.from_tuples([("A", 2), ("A", 2), ("B", 3), ("B", 4)])
    idx2 = MultiIndex.from_tuples([("A", 1), ("A", 2), ("B", 3), ("B", 4)])
    kwargs = dict(check_exact=check_exact, check_less_precise=check_less_precise)

    msg = """MultiIndex level \\[1\\] are different

MultiIndex level \\[1\\] values are different \\(25\\.0 %\\)
\\[left\\]:  Int64Index\\(\\[2, 2, 3, 4\\], dtype='int64'\\)
\\[right\\]: Int64Index\\(\\[1, 2, 3, 4\\], dtype='int64'\\)"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_index_equal(idx1, idx2, **kwargs)


@pytest.mark.parametrize(
    "name1,name2",
    [(None, "x"), ("x", "x"), (np.nan, np.nan), (NaT, NaT), (np.nan, NaT)],
)
def test_index_equal_names(name1, name2):
    msg = """Index are different

Attribute "names" are different
\\[left\\]:  \\[{name1}\\]
\\[right\\]: \\[{name2}\\]"""

    idx1 = Index([1, 2, 3], name=name1)
    idx2 = Index([1, 2, 3], name=name2)

    if name1 == name2 or name1 is name2:
        tm.assert_index_equal(idx1, idx2)
    else:
        name1 = "'x'" if name1 == "x" else name1
        name2 = "'x'" if name2 == "x" else name2
        msg = msg.format(name1=name1, name2=name2)

        with pytest.raises(AssertionError, match=msg):
            tm.assert_index_equal(idx1, idx2)


def test_index_equal_category_mismatch(check_categorical):
    msg = """Index are different

Attribute "dtype" are different
\\[left\\]:  CategoricalDtype\\(categories=\\['a', 'b'\\], ordered=False\\)
\\[right\\]: CategoricalDtype\\(categories=\\['a', 'b', 'c'\\], \
ordered=False\\)"""

    idx1 = Index(Categorical(["a", "b"]))
    idx2 = Index(Categorical(["a", "b"], categories=["a", "b", "c"]))

    if check_categorical:
        with pytest.raises(AssertionError, match=msg):
            tm.assert_index_equal(idx1, idx2, check_categorical=check_categorical)
    else:
        tm.assert_index_equal(idx1, idx2, check_categorical=check_categorical)
