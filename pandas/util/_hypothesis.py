"""
This module houses utility functions to generate hypothesis strategies which
 can be used to generate random input test data for various test cases.
It is for internal use by different test case files like pandas/test/test*.py
 files only and should not be used beyond this purpose.
For more information on hypothesis, check
(http://hypothesis.readthedocs.io/en/latest/).
"""
import string
from hypothesis import (given,  # noqa:F401
                        settings,   # noqa:F401
                        assume, # noqa:F401
                        strategies as st,
                        )


def get_elements(elem_type):
    """
    Helper function to return hypothesis strategy whose elements depends on
    the input data-type.
    Currently only four types are supported namely, bool, int, float and str.

    Parameters
    ----------
    elem_type: type
        type of the elements for the strategy.

    Returns
    -------
    hypothesis strategy.

    Examples
    --------
    >>> strat = get_elements(str)
    >>> strat.example()
    'KWAo'

    >>> strat.example()
    'OfAlBH'

    >>> strat = get_elements(int)
    >>> strat.example()
    31911

    >>> strat.example()
    25288

    >>> strat = get_elements(float)
    >>> strat.example()
    nan

    >>> strat.example()
    inf

    >>> strat.example()
    -2.2250738585072014e-308

    >>> strat.example()
    0.5

    >>> strat.example()
    1.7976931348623157e+308

    >>> strat = get_elements(bool)
    >>> strat.example()
    True

    >>> strat.example()
    True

    >>> strat.example()
    False
    """
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
    >>> seq_strategy = get_seq((int, str, bool), mixed=True, min_size=1,
...    max_size=5)

    >>> seq_strategy.example()
    ['lkYMSn', -2501, 35, 'J']

    >>> seq_strategy.example()
    [True]

    >>> seq_strategy.example()
    ['dRWgQYrBrW', True, False, 'gmsujJVDBM', 'Z']

    >>> seq_strategy = get_seq((int, bool),
...                             mixed=False,
...                             min_size=1,
...                             max_size=5,
...                             transform_func=lambda seq:
...                             [str(x) for x in seq])

    >>> seq_strategy.example()
    ['9552', '124', '-24024']

    >>> seq_strategy.example()
    ['-1892']

    >>> seq_strategy.example()
    ['22', '66', '14785', '-26312', '32']
    """
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
