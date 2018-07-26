# -*- coding: utf-8 -*-

import pytest
import os
import collections
from functools import partial

import numpy as np

import pandas
from pandas import Series, DataFrame, Timestamp
from pandas.compat import range, lmap
import pandas.core.common as com
from pandas.core import ops
from pandas.io.common import (
    _compression_to_extension,
    _get_handle,
)
import pandas.util.testing as tm


def test_mut_exclusive():
    msg = "mutually exclusive arguments: '[ab]' and '[ab]'"
    with tm.assert_raises_regex(TypeError, msg):
        com._mut_exclusive(a=1, b=2)
    assert com._mut_exclusive(a=1, b=None) == 1
    assert com._mut_exclusive(major=None, major_axis=None) is None
    assert com._mut_exclusive(a=None, b=2) == 2


def test_get_callable_name():
    getname = com._get_callable_name

    def fn(x):
        return x

    lambda_ = lambda x: x
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


def test_iterpairs():
    data = [1, 2, 3, 4]
    expected = [(1, 2), (2, 3), (3, 4)]

    result = list(com.iterpairs(data))

    assert (result == expected)


def test_split_ranges():
    def _bin(x, width):
        "return int(x) as a base2 string of given width"
        return ''.join(str((x >> i) & 1) for i in range(width - 1, -1, -1))

    def test_locs(mask):
        nfalse = sum(np.array(mask) == 0)

        remaining = 0
        for s, e in com.split_ranges(mask):
            remaining += e - s

            assert 0 not in mask[s:e]

        # make sure the total items covered by the ranges are a complete cover
        assert remaining + nfalse == len(mask)

    # exhaustively test all possible mask sequences of length 8
    ncols = 8
    for i in range(2 ** ncols):
        cols = lmap(int, list(_bin(i, ncols)))  # count up in base2
        mask = [cols[i] == 1 for i in range(len(cols))]
        test_locs(mask)

    # base cases
    test_locs([])
    test_locs([0])
    test_locs([1])


def test_map_indices_py():
    data = [4, 3, 2, 1]
    expected = {4: 0, 3: 1, 2: 2, 1: 3}

    result = com.map_indices_py(data)

    assert (result == expected)


def test_union():
    a = [1, 2, 3]
    b = [4, 5, 6]

    union = sorted(com.union(a, b))

    assert ((a + b) == union)


def test_difference():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(com.difference(b, a))

    assert ([4, 5, 6] == inter)


def test_intersection():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(com.intersection(a, b))

    assert (a == inter)


def test_groupby():
    values = ['foo', 'bar', 'baz', 'baz2', 'qux', 'foo3']
    expected = {'f': ['foo', 'foo3'],
                'b': ['bar', 'baz', 'baz2'],
                'q': ['qux']}

    grouped = com.groupby(values, lambda x: x[0])

    for k, v in grouped:
        assert v == expected[k]


def test_random_state():
    import numpy.random as npr
    # Check with seed
    state = com._random_state(5)
    assert state.uniform() == npr.RandomState(5).uniform()

    # Check with random state object
    state2 = npr.RandomState(10)
    assert com._random_state(state2).uniform() == npr.RandomState(10).uniform()

    # check with no arg random state
    assert com._random_state() is np.random

    # Error for floats or strings
    with pytest.raises(ValueError):
        com._random_state('test')

    with pytest.raises(ValueError):
        com._random_state(5.5)


@pytest.mark.parametrize('left, right, expected', [
    (Series([1], name='x'), Series([2], name='x'), 'x'),
    (Series([1], name='x'), Series([2], name='y'), None),
    (Series([1]), Series([2], name='x'), None),
    (Series([1], name='x'), Series([2]), None),
    (Series([1], name='x'), [2], 'x'),
    ([1], Series([2], name='y'), 'y')])
def test_maybe_match_name(left, right, expected):
    assert ops._maybe_match_name(left, right) == expected


def test_dict_compat():
    data_datetime64 = {np.datetime64('1990-03-15'): 1,
                       np.datetime64('2015-03-15'): 2}
    data_unchanged = {1: 2, 3: 4, 5: 6}
    expected = {Timestamp('1990-3-15'): 1, Timestamp('2015-03-15'): 2}
    assert (com._dict_compat(data_datetime64) == expected)
    assert (com._dict_compat(expected) == expected)
    assert (com._dict_compat(data_unchanged) == data_unchanged)


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


@pytest.mark.parametrize('obj', [
    DataFrame(100 * [[0.123456, 0.234567, 0.567567],
                     [12.32112, 123123.2, 321321.2]],
              columns=['X', 'Y', 'Z']),
    Series(100 * [0.123456, 0.234567, 0.567567], name='X')])
@pytest.mark.parametrize('method', ['to_pickle', 'to_json', 'to_csv'])
def test_compression_size(obj, method, compression_only):
    # Tests that compression is occurring by comparing to the bytes on disk of
    # the uncompressed file.
    with tm.ensure_clean() as filename:
        getattr(obj, method)(filename, compression=compression_only)
        compressed = os.path.getsize(filename)
        getattr(obj, method)(filename, compression=None)
        uncompressed = os.path.getsize(filename)
        assert uncompressed > compressed


@pytest.mark.parametrize('obj', [
    DataFrame(100 * [[0.123456, 0.234567, 0.567567],
                     [12.32112, 123123.2, 321321.2]],
              columns=['X', 'Y', 'Z']),
    Series(100 * [0.123456, 0.234567, 0.567567], name='X')])
@pytest.mark.parametrize('method', ['to_csv', 'to_json'])
def test_compression_size_fh(obj, method, compression_only):

    with tm.ensure_clean() as filename:
        f, _handles = _get_handle(filename, 'w', compression=compression_only)
        with f:
            getattr(obj, method)(f)
            assert not f.closed
        assert f.closed
        compressed = os.path.getsize(filename)
    with tm.ensure_clean() as filename:
        f, _handles = _get_handle(filename, 'w', compression=None)
        with f:
            getattr(obj, method)(f)
            assert not f.closed
        assert f.closed
        uncompressed = os.path.getsize(filename)
        assert uncompressed > compressed


@pytest.mark.parametrize('write_method, write_kwargs, read_method', [
    ('to_csv', {'index': False}, pandas.read_csv),
    ('to_json', {}, pandas.read_json),
    ('to_pickle', {}, pandas.read_pickle),
])
def test_dataframe_compression_defaults_to_infer(
        write_method, write_kwargs, read_method, compression_only):
    # Test that DataFrame.to_* methods default to inferring compression from
    # paths. https://github.com/pandas-dev/pandas/pull/22011
    input = DataFrame([[1.0, 0, -4.4], [3.4, 5, 2.4]], columns=['X', 'Y', 'Z'])
    extension = _compression_to_extension[compression_only]
    with tm.ensure_clean('compressed' + extension) as path:
        # assumes that compression='infer' is the default
        getattr(input, write_method)(path, **write_kwargs)
        output = read_method(path, compression=compression_only)
    tm.assert_frame_equal(output, input)


@pytest.mark.parametrize('write_method,write_kwargs,read_method,read_kwargs', [
    ('to_csv', {'index': False, 'header': True}, pandas.read_csv, {'squeeze': True}),
    ('to_json', {}, pandas.read_json, {'typ': 'series'}),
    ('to_pickle', {}, pandas.read_pickle, {}),
])
def test_series_compression_defaults_to_infer(
        write_method, write_kwargs, read_method, read_kwargs, compression_only):
    # Test that Series.to_* methods default to inferring compression from
    # paths. https://github.com/pandas-dev/pandas/pull/22011
    input = Series([0, 5, -2, 10], name='X')
    extension = _compression_to_extension[compression_only]
    with tm.ensure_clean('compressed' + extension) as path:
        # assumes that compression='infer' is the default
        getattr(input, write_method)(path, **write_kwargs)
        output = read_method(path, compression=compression_only, **read_kwargs)
    tm.assert_series_equal(output, input)


def test_compression_warning(compression_only):
    # Assert that passing a file object to to_csv while explicitly specifying a
    # compression protocol triggers a RuntimeWarning, as per
    # https://github.com/pandas-dev/pandas/issues/21227.
    # Note that pytest has an issue that causes assert_produces_warning to fail
    # in Python 2 if the warning has occurred in previous tests
    # (see https://git.io/fNEBm & https://git.io/fNEBC). Hence, should this
    # test fail in just Python 2 builds, it likely indicates that other tests
    # are producing RuntimeWarnings, thereby triggering the pytest bug.
    df = DataFrame(100 * [[0.123456, 0.234567, 0.567567],
                          [12.32112, 123123.2, 321321.2]],
                   columns=['X', 'Y', 'Z'])
    with tm.ensure_clean() as filename:
        f, _handles = _get_handle(filename, 'w', compression=compression_only)
        with tm.assert_produces_warning(RuntimeWarning,
                                        check_stacklevel=False):
            with f:
                df.to_csv(f, compression=compression_only)
