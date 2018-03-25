import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm


@pytest.mark.parametrize('ordered', [True, False])
@pytest.mark.parametrize('categories', [
    ['b', 'a', 'c'],
    ['a', 'b', 'c', 'd'],
])
def test_factorize(categories, ordered):
    cat = pd.Categorical(['b', 'b', 'a', 'c', None],
                         categories=categories,
                         ordered=ordered)
    labels, uniques = pd.factorize(cat)
    expected_labels = np.array([0, 0, 1, 2, -1], dtype=np.intp)
    expected_uniques = pd.Categorical(['b', 'a', 'c'],
                                      categories=categories,
                                      ordered=ordered)

    tm.assert_numpy_array_equal(labels, expected_labels)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_factorized_sort():
    cat = pd.Categorical(['b', 'b', None, 'a'])
    labels, uniques = pd.factorize(cat, sort=True)
    expected_labels = np.array([1, 1, -1, 0], dtype=np.intp)
    expected_uniques = pd.Categorical(['a', 'b'])

    tm.assert_numpy_array_equal(labels, expected_labels)
    tm.assert_categorical_equal(uniques, expected_uniques)


def test_factorized_sort_ordered():
    cat = pd.Categorical(['b', 'b', None, 'a'],
                         categories=['c', 'b', 'a'],
                         ordered=True)

    labels, uniques = pd.factorize(cat, sort=True)
    expected_labels = np.array([0, 0, -1, 1], dtype=np.intp)
    expected_uniques = pd.Categorical(['b', 'a'],
                                      categories=['c', 'b', 'a'],
                                      ordered=True)

    tm.assert_numpy_array_equal(labels, expected_labels)
    tm.assert_categorical_equal(uniques, expected_uniques)
