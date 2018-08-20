
import numpy as np
from pandas import date_range, Index
import pandas.util.testing as tm
from pandas.core.reshape.util import cartesian_product

import string
from datetime import date
from dateutil import relativedelta

from pandas.util import _hypothesis as hp

NO_OF_EXAMPLES_PER_TEST_CASE = 20


class TestCartesianProduct(object):

    @hp.settings(max_examples=20)
    @hp.given(hp.st.lists(hp.st.text(string.ascii_letters,
                                     min_size=1, max_size=1),
                          min_size=1, max_size=3),
              hp.get_seq((int,), False, 1, 2))
    def test_simple(self, x, y):
        result1, result2 = cartesian_product([x, y])
        expected1 = np.array([item1 for item1 in x for item2 in y])
        expected2 = np.array([item2 for item1 in x for item2 in y])

        tm.assert_numpy_array_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    @hp.settings(max_examples=20)
    @hp.given(hp.st.dates(min_value=date(1900, 1, 1),
                          max_value=date(2100, 1, 1)))
    def test_datetimeindex(self, d):
        # regression test for GitHub issue #6439
        # make sure that the ordering on datetimeindex is consistent
        n = d + relativedelta.relativedelta(days=1)
        x = date_range(d, periods=2)
        result1, result2 = [Index(y).day for y in cartesian_product([x, x])]
        expected1 = Index([d.day, d.day, n.day, n.day])
        expected2 = Index([d.day, n.day, d.day, n.day])

        tm.assert_index_equal(result1, expected1)
        tm.assert_index_equal(result2, expected2)

    @hp.settings(max_examples=20)
    @hp.given(hp.st.lists(hp.st.nothing()),
              hp.get_seq((int,), False, min_size=1, max_size=10),
              hp.get_seq((str,), False, min_size=1, max_size=10))
    def test_empty(self, empty_list, list_of_int, list_of_str):
        # product of empty factors
        X = [empty_list, list_of_int, empty_list]
        Y = [empty_list, empty_list, list_of_str]

        for x, y in zip(X, Y):
            expected1 = np.array([], dtype=np.asarray(x).dtype)
            expected2 = np.array([], dtype=np.asarray(y).dtype)
            result1, result2 = cartesian_product([x, y])
            tm.assert_numpy_array_equal(result1, expected1)
            tm.assert_numpy_array_equal(result2, expected2)

        # empty product (empty input):
        result = cartesian_product(empty_list)
        expected = []
        assert result == expected

    @hp.settings(max_examples=20)
    @hp.given(hp.st.integers(),
              hp.st.text(string.ascii_letters, min_size=1),
              hp.get_seq((int, str), True, min_size=1),
              hp.st.builds(lambda *x: list(x), hp.st.integers(),
                           hp.st.text(string.ascii_letters, min_size=1),
                           hp.st.lists(hp.st.integers(), min_size=1)))
    def test_invalid_input(self, number, text, seq, mixed_seq):

        invalid_inputs = [number,
                          text,
                          seq,
                          mixed_seq]

        msg = "Input must be a list-like of list-likes"
        for X in invalid_inputs:
            tm.assert_raises_regex(TypeError, msg, cartesian_product, X=X)
