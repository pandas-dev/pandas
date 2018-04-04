
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
def get_seq(draw, types, mixed=False, min_size=None, max_size=None,
            transform_func=None):
    """
    Helper function to generate strategy for creating lists.
    What constitute in the generated list is driven by the different
    parameters.

    Parameters
    ----------
    types: iterable sequence like tuple or list
        types which can be in the generated list.
    mixed: bool
        if True, list will contains elements from all types listed in arg,
        otherwise it will have elements only from types[0].
    min_size: int
        minimum size of the list.
    max_size: int
        maximum size of the list.
    transform_func: callable
        a callable which can be applied to whole list after it has been
         generated. It can think of as providing functionality of filter
         and map function.

    Returns
    -------
    hypothesis lists strategy.

    Examples
    --------
    seq_strategy = get_seq((int, str, bool),
                            mixed=True, min_size=1, max_size=5)
    seq_strategy.example()
    Out[12]: ['lkYMSn', -2501, 35, 'J']
    seq_strategy.example()
    Out[13]: [True]
    seq_strategy.example()
    Out[14]: ['dRWgQYrBrW', True, False, 'gmsujJVDBM', 'Z']

    seq_strategy = get_seq((int, bool),
                            mixed=False,
                            min_size=1,
                            max_size=5,
                            transform_func=lambda seq: [str(x) for x in seq])
    seq_strategy.example()
    Out[19]: ['-1892']
    seq_strategy.example()
    Out[20]: ['22', '66', '14785', '-26312', '32']
    seq_strategy.example()
    Out[21]: ['22890', '-15537', '96']
    """
    strategy = st.nothing()
    if min_size is None:
        min_size = draw(st.integers(min_value=0, max_value=100))

    if max_size is None:
        max_size = draw(st.integers(min_value=min_size, max_value=100))

    assert min_size <= max_size, \
        'max_size must be greater than equal to min_size'

    elem_strategies = []
    for elem_type in types:
        elem_strategies.append(get_elements(elem_type))
        if not mixed:
            break

    if transform_func:
        strategy = draw(st.lists(st.one_of(elem_strategies),
                                 min_size=min_size,
                                 max_size=max_size).map(transform_func))
    else:
        strategy = draw(st.lists(st.one_of(elem_strategies),
                                 min_size=min_size,
                                 max_size=max_size))
    return strategy


class TestCartesianProduct(object):

    @settings(max_examples=NO_OF_EXAMPLES_PER_TEST_CASE)
    @given(get_seq((str,), False, 1, 1),
           get_seq((int,), False, 1, 2))
    def test_simple(self, x, y):
        x = list(x[0])
        # non-empty test case is handled in test_empty,
        # therefore ignore it here.
        assume(len(x) != 0)
        result1, result2 = cartesian_product([x, y])
        expected1 = np.array([item1 for item1 in x for item2 in y])
        expected2 = np.array([item2 for item1 in x for item2 in y])

        tm.assert_numpy_array_equal(result1, expected1)
        tm.assert_numpy_array_equal(result2, expected2)

    @settings(max_examples=NO_OF_EXAMPLES_PER_TEST_CASE)
    @given(st.dates(min_value=date(1900, 1, 1), max_value=date(2100, 1, 1)))
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
    @given(st.integers(),
           st.text(string.ascii_letters, min_size=1),
           get_seq((int, str), True, min_size=1),
           st.builds(lambda *x: list(x), st.integers(),
                     st.text(string.ascii_letters, min_size=1),
                     st.lists(st.integers(), min_size=1)))
    def test_invalid_input(self, number, text, seq, mixed_seq):

        invalid_inputs = [number,
                          text,
                          seq,
                          mixed_seq]

        msg = "Input must be a list-like of list-likes"
        for X in invalid_inputs:
            tm.assert_raises_regex(TypeError, msg, cartesian_product, X=X)
