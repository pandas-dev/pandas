import string
from hypothesis import (given,
                        settings,
                        assume,
                        strategies as st,
                        )


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
