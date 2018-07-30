# -*- coding: utf-8 -*-

import collections
from functools import partial

import numpy as np
import pytest

from pandas import Series, Timestamp
from pandas.core import (
    _maybe_match_name,
    common as com,
)


def test_get_callable_name():
    getname = com.get_callable_name

    def fn(x):
        return x

    lambda_ = lambda x: x  # noqa: E731
    part1 = partial(fn)
    part2 = partial(part1)

    class somecall(object):

        def __call__(self):
            return x  # noqa

    assert getname(fn) == 'fn'
    assert getname(lambda_)
    assert getname(part1) == 'fn'
    assert getname(part2) == 'fn'
    assert getname(somecall()) == 'somecall'
    assert getname(1) is None


def test_any_none():
    assert (com._any_none(1, 2, 3, None))
    assert (not com._any_none(1, 2, 3, 4))


def test_all_not_none():
    assert (com._all_not_none(1, 2, 3, 4))
    assert (not com._all_not_none(1, 2, 3, None))
    assert (not com._all_not_none(None, None, None, None))


def test_random_state():
    import numpy.random as npr
    # Check with seed
    state = com.random_state(5)
    assert state.uniform() == npr.RandomState(5).uniform()

    # Check with random state object
    state2 = npr.RandomState(10)
    assert com.random_state(state2).uniform() == npr.RandomState(10).uniform()

    # check with no arg random state
    assert com.random_state() is np.random

    # Error for floats or strings
    with pytest.raises(ValueError):
        com.random_state('test')

    with pytest.raises(ValueError):
        com.random_state(5.5)


@pytest.mark.parametrize('left, right, expected', [
    (Series([1], name='x'), Series([2], name='x'), 'x'),
    (Series([1], name='x'), Series([2], name='y'), None),
    (Series([1]), Series([2], name='x'), None),
    (Series([1], name='x'), Series([2]), None),
    (Series([1], name='x'), [2], 'x'),
    ([1], Series([2], name='y'), 'y')])
def test_maybe_match_name(left, right, expected):
    assert _maybe_match_name(left, right) == expected


def test_dict_compat():
    data_datetime64 = {np.datetime64('1990-03-15'): 1,
                       np.datetime64('2015-03-15'): 2}
    data_unchanged = {1: 2, 3: 4, 5: 6}
    expected = {Timestamp('1990-3-15'): 1, Timestamp('2015-03-15'): 2}
    assert (com.dict_compat(data_datetime64) == expected)
    assert (com.dict_compat(expected) == expected)
    assert (com.dict_compat(data_unchanged) == data_unchanged)


def test_standardize_mapping():
    # No uninitialized defaultdicts
    with pytest.raises(TypeError):
        com.standardize_mapping(collections.defaultdict)

    # No non-mapping subtypes, instance
    with pytest.raises(TypeError):
        com.standardize_mapping([])

    # No non-mapping subtypes, class
    with pytest.raises(TypeError):
        com.standardize_mapping(list)

    fill = {'bad': 'data'}
    assert (com.standardize_mapping(fill) == dict)

    # Convert instance to type
    assert (com.standardize_mapping({}) == dict)

    dd = collections.defaultdict(list)
    assert isinstance(com.standardize_mapping(dd), partial)
