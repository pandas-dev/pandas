
import numpy as np
from pandas import date_range, Index
import pandas.util.testing as tm
from pandas.core.reshape.util import cartesian_product

from hypothesis import strategies as st
from hypothesis import given, settings, assume
from datetime import date
from dateutil import relativedelta
import string


NO_OF_EXAMPLES_PER_TEST_CASE = 20


def get_elements(elem_type):
    strategy = st.nothing()
    if elem_type == bool:
        strategy = st.booleans()
    elif elem_type == int:
        strategy = st.integers()
    elif elem_type == float:
        strategy = st.floats()
    elif elem_type == str:
        strategy = st.text(string.ascii_letters, max_size=10)
    return strategy


@st.composite
def get_seq(draw, types, mixed=False, min_size=None, max_size=None, transform_func=None):
    """helper function to generate strategy for creating lists. parameters define the nature of to be generated list.
    :param types: what type of elements constitute the list
    :param mixed: if True, list will contains elements from all types listed in arg, oterwise it will have elements only from types[0].
    :param min_size: minimum size of the list.
    :param max_size: maximum size of the list.
    :param transform_func: a callable which can be applied to whole list after it has been generated.
    """
    strategy = st.nothing()
    if min_size is None:
        min_size = draw(st.integers(min_value=0, max_value=100))

    if max_size is None:
        max_size = draw(st.integers(min_value=min_size, max_value=100))

    assert min_size <= max_size, 'max_size must be greater than equal to min_size'

    elem_strategies = []
    for elem_type in types:
        elem_strategies.append(get_elements(elem_type))
        if not mixed:
            break

    if transform_func:
        strategy = draw(st.lists(st.one_of(elem_strategies),
                                 min_size=min_size, max_size=max_size).map(transform_func))
    else:
        strategy = draw(st.lists(st.one_of(elem_strategies),
                                 min_size=min_size, max_size=max_size))
    return strategy


class TestCartesianProduct(object):

    @settings(max_examples=NO_OF_EXAMPLES_PER_TEST_CASE)
    @given(get_seq((str,), False, 1, 1),
           get_seq((int,), False, 1, 2))
    def test_simple(self, x, y):
        x = list(x[0])
        # non-empty test case is handled in test_empty, therefore ignore it here
        assume(len(x) != 0)
        result1, result2 = cartesian_product([x, y])
        expected1 = np.array([item1 for item1 in x for item2 in y])
        expected2 = np.array([item2 for item1 in x for item2 in y])

        tm.assert_numpy_array_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    @settings(max_examples=NO_OF_EXAMPLES_PER_TEST_CASE)
    def test_datetimeindex(self):
        # regression test for GitHub issue #6439
        # make sure that the ordering on datetimeindex is consistent
        d = st.dates(min_value=date(1900, 1, 1), max_value=date(2100, 1, 1)).example()
        n = d + relativedelta.relativedelta(days=1)
        x = date_range(d, periods=2)
        result1, result2 = [Index(y).day for y in cartesian_product([x, x])]
        expected1 = Index([d.day, d.day, n.day, n.day])
        expected2 = Index([d.day, n.day, d.day, n.day])

        tm.assert_index_equal(result1, expected1)
        tm.assert_index_equal(result2, expected2)

    @settings(max_examples=NO_OF_EXAMPLES_PER_TEST_CASE)
    @given(st.lists(st.nothing()),
           get_seq((int,), False),
           get_seq((str,), False))
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

    @settings(max_examples=NO_OF_EXAMPLES_PER_TEST_CASE)
    def test_invalid_input(self):
        invalid_inputs = [st.integers().example(),
                          st.tuples(st.integers()).example(),
                          st.tuples(st.integers(), st.integers()).example(),
                          st.text(string.ascii_letters, min_size=1, max_size=1).example(),
                          st.tuples(st.text(string.ascii_letters, min_size=1, max_size=1)).example(),
                          st.tuples(st.text(string.ascii_letters, min_size=1, max_size=1),
                                    st.text(string.ascii_letters, min_size=1, max_size=1)).example(),
                          st.tuples(st.tuples(st.text(string.ascii_letters, min_size=1, max_size=1)),
                                    st.text(string.ascii_letters, min_size=1, max_size=1)).example()]

        msg = "Input must be a list-like of list-likes"
        for X in invalid_inputs:
            tm.assert_raises_regex(TypeError, msg, cartesian_product, X=X)
