# -*- coding: utf-8 -*-

import numpy as np
import pytest

import pandas as pd
from pandas import MultiIndex, Series
import pandas.util.testing as tm


@pytest.mark.parametrize("case", [0.5, "xxx"])
@pytest.mark.parametrize("sort", [None, False])
@pytest.mark.parametrize("method", ["intersection", "union",
                                    "difference", "symmetric_difference"])
def test_set_ops_error_cases(idx, case, sort, method):
    # non-iterable input
    msg = "Input must be Index or array-like"
    with pytest.raises(TypeError, match=msg):
        getattr(idx, method)(case, sort=sort)


@pytest.mark.parametrize("sort", [None, False])
def test_intersection_base(idx, sort):
    first = idx[:5]
    second = idx[:3]
    intersect = first.intersection(second, sort=sort)

    if sort is None:
        tm.assert_index_equal(intersect, second.sort_values())
    assert tm.equalContents(intersect, second)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.intersection(case, sort=sort)
        if sort is None:
            tm.assert_index_equal(result, second.sort_values())
        assert tm.equalContents(result, second)

    msg = "other must be a MultiIndex or a list of tuples"
    with pytest.raises(TypeError, match=msg):
        first.intersection([1, 2, 3], sort=sort)


@pytest.mark.parametrize("sort", [None, False])
def test_union_base(idx, sort):
    first = idx[3:]
    second = idx[:5]
    everything = idx
    union = first.union(second, sort=sort)
    if sort is None:
        tm.assert_index_equal(union, everything.sort_values())
    assert tm.equalContents(union, everything)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.union(case, sort=sort)
        if sort is None:
            tm.assert_index_equal(result, everything.sort_values())
        assert tm.equalContents(result, everything)

    msg = "other must be a MultiIndex or a list of tuples"
    with pytest.raises(TypeError, match=msg):
        first.union([1, 2, 3], sort=sort)


@pytest.mark.parametrize("sort", [None, False])
def test_difference_base(idx, sort):
    second = idx[4:]
    answer = idx[:4]
    result = idx.difference(second, sort=sort)

    if sort is None:
        answer = answer.sort_values()

    assert result.equals(answer)
    tm.assert_index_equal(result, answer)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = idx.difference(case, sort=sort)
        tm.assert_index_equal(result, answer)

    msg = "other must be a MultiIndex or a list of tuples"
    with pytest.raises(TypeError, match=msg):
        idx.difference([1, 2, 3], sort=sort)


@pytest.mark.parametrize("sort", [None, False])
def test_symmetric_difference(idx, sort):
    first = idx[1:]
    second = idx[:-1]
    answer = idx[[-1, 0]]
    result = first.symmetric_difference(second, sort=sort)

    if sort is None:
        answer = answer.sort_values()

    tm.assert_index_equal(result, answer)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.symmetric_difference(case, sort=sort)
        tm.assert_index_equal(result, answer)

    msg = "other must be a MultiIndex or a list of tuples"
    with pytest.raises(TypeError, match=msg):
        first.symmetric_difference([1, 2, 3], sort=sort)


def test_empty(idx):
    # GH 15270
    assert not idx.empty
    assert idx[:0].empty


@pytest.mark.parametrize("sort", [None, False])
def test_difference(idx, sort):

    first = idx
    result = first.difference(idx[-3:], sort=sort)
    vals = idx[:-3].values

    if sort is None:
        vals = sorted(vals)

    expected = MultiIndex.from_tuples(vals,
                                      sortorder=0,
                                      names=idx.names)

    assert isinstance(result, MultiIndex)
    assert result.equals(expected)
    assert result.names == idx.names
    tm.assert_index_equal(result, expected)

    # empty difference: reflexive
    result = idx.difference(idx, sort=sort)
    expected = idx[:0]
    assert result.equals(expected)
    assert result.names == idx.names

    # empty difference: superset
    result = idx[-3:].difference(idx, sort=sort)
    expected = idx[:0]
    assert result.equals(expected)
    assert result.names == idx.names

    # empty difference: degenerate
    result = idx[:0].difference(idx, sort=sort)
    expected = idx[:0]
    assert result.equals(expected)
    assert result.names == idx.names

    # names not the same
    chunklet = idx[-3:]
    chunklet.names = ['foo', 'baz']
    result = first.difference(chunklet, sort=sort)
    assert result.names == (None, None)

    # empty, but non-equal
    result = idx.difference(idx.sortlevel(1)[0], sort=sort)
    assert len(result) == 0

    # raise Exception called with non-MultiIndex
    result = first.difference(first.values, sort=sort)
    assert result.equals(first[:0])

    # name from empty array
    result = first.difference([], sort=sort)
    assert first.equals(result)
    assert first.names == result.names

    # name from non-empty array
    result = first.difference([('foo', 'one')], sort=sort)
    expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'), (
        'foo', 'two'), ('qux', 'one'), ('qux', 'two')])
    expected.names = first.names
    assert first.names == result.names

    msg = "other must be a MultiIndex or a list of tuples"
    with pytest.raises(TypeError, match=msg):
        first.difference([1, 2, 3, 4, 5], sort=sort)


def test_difference_sort_special():
    # GH-24959
    idx = pd.MultiIndex.from_product([[1, 0], ['a', 'b']])
    # sort=None, the default
    result = idx.difference([])
    tm.assert_index_equal(result, idx)


@pytest.mark.xfail(reason="Not implemented.")
def test_difference_sort_special_true():
    # TODO decide on True behaviour
    idx = pd.MultiIndex.from_product([[1, 0], ['a', 'b']])
    result = idx.difference([], sort=True)
    expected = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
    tm.assert_index_equal(result, expected)


def test_difference_sort_incomparable():
    # GH-24959
    idx = pd.MultiIndex.from_product([[1, pd.Timestamp('2000'), 2],
                                      ['a', 'b']])

    other = pd.MultiIndex.from_product([[3, pd.Timestamp('2000'), 4],
                                        ['c', 'd']])
    # sort=None, the default
    # MultiIndex.difference deviates here from other difference
    # implementations in not catching the TypeError
    with pytest.raises(TypeError):
        result = idx.difference(other)

    # sort=False
    result = idx.difference(other, sort=False)
    tm.assert_index_equal(result, idx)


@pytest.mark.xfail(reason="Not implemented.")
def test_difference_sort_incomparable_true():
    # TODO decide on True behaviour
    # # sort=True, raises
    idx = pd.MultiIndex.from_product([[1, pd.Timestamp('2000'), 2],
                                      ['a', 'b']])
    other = pd.MultiIndex.from_product([[3, pd.Timestamp('2000'), 4],
                                        ['c', 'd']])

    with pytest.raises(TypeError):
        idx.difference(other, sort=True)


@pytest.mark.parametrize("sort", [None, False])
def test_union(idx, sort):
    piece1 = idx[:5][::-1]
    piece2 = idx[3:]

    the_union = piece1.union(piece2, sort=sort)

    if sort is None:
        tm.assert_index_equal(the_union, idx.sort_values())

    assert tm.equalContents(the_union, idx)

    # corner case, pass self or empty thing:
    the_union = idx.union(idx, sort=sort)
    assert the_union is idx

    the_union = idx.union(idx[:0], sort=sort)
    assert the_union is idx

    # won't work in python 3
    # tuples = _index.values
    # result = _index[:4] | tuples[4:]
    # assert result.equals(tuples)

    # not valid for python 3
    # def test_union_with_regular_index(self):
    #     other = Index(['A', 'B', 'C'])

    #     result = other.union(idx)
    #     assert ('foo', 'one') in result
    #     assert 'B' in result

    #     result2 = _index.union(other)
    #     assert result.equals(result2)


@pytest.mark.parametrize("sort", [None, False])
def test_intersection(idx, sort):
    piece1 = idx[:5][::-1]
    piece2 = idx[3:]

    the_int = piece1.intersection(piece2, sort=sort)

    if sort is None:
        tm.assert_index_equal(the_int, idx[3:5])
    assert tm.equalContents(the_int, idx[3:5])

    # corner case, pass self
    the_int = idx.intersection(idx, sort=sort)
    assert the_int is idx

    # empty intersection: disjoint
    empty = idx[:2].intersection(idx[2:], sort=sort)
    expected = idx[:0]
    assert empty.equals(expected)

    # can't do in python 3
    # tuples = _index.values
    # result = _index & tuples
    # assert result.equals(tuples)


def test_intersect_equal_sort():
    # GH-24959
    idx = pd.MultiIndex.from_product([[1, 0], ['a', 'b']])
    tm.assert_index_equal(idx.intersection(idx, sort=False), idx)
    tm.assert_index_equal(idx.intersection(idx, sort=None), idx)


@pytest.mark.xfail(reason="Not implemented.")
def test_intersect_equal_sort_true():
    # TODO decide on True behaviour
    idx = pd.MultiIndex.from_product([[1, 0], ['a', 'b']])
    sorted_ = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
    tm.assert_index_equal(idx.intersection(idx, sort=True), sorted_)


@pytest.mark.parametrize('slice_', [slice(None), slice(0)])
def test_union_sort_other_empty(slice_):
    # https://github.com/pandas-dev/pandas/issues/24959
    idx = pd.MultiIndex.from_product([[1, 0], ['a', 'b']])

    # default, sort=None
    other = idx[slice_]
    tm.assert_index_equal(idx.union(other), idx)
    # MultiIndex does not special case empty.union(idx)
    # tm.assert_index_equal(other.union(idx), idx)

    # sort=False
    tm.assert_index_equal(idx.union(other, sort=False), idx)


@pytest.mark.xfail(reason="Not implemented.")
def test_union_sort_other_empty_sort(slice_):
    # TODO decide on True behaviour
    # # sort=True
    idx = pd.MultiIndex.from_product([[1, 0], ['a', 'b']])
    other = idx[:0]
    result = idx.union(other, sort=True)
    expected = pd.MultiIndex.from_product([[0, 1], ['a', 'b']])
    tm.assert_index_equal(result, expected)


def test_union_sort_other_incomparable():
    # https://github.com/pandas-dev/pandas/issues/24959
    idx = pd.MultiIndex.from_product([[1, pd.Timestamp('2000')], ['a', 'b']])

    # default, sort=None
    result = idx.union(idx[:1])
    tm.assert_index_equal(result, idx)

    # sort=False
    result = idx.union(idx[:1], sort=False)
    tm.assert_index_equal(result, idx)


@pytest.mark.xfail(reason="Not implemented.")
def test_union_sort_other_incomparable_sort():
    # TODO decide on True behaviour
    # # sort=True
    idx = pd.MultiIndex.from_product([[1, pd.Timestamp('2000')], ['a', 'b']])
    with pytest.raises(TypeError, match='Cannot compare'):
        idx.union(idx[:1], sort=True)


@pytest.mark.parametrize("method", ['union', 'intersection', 'difference',
                                    'symmetric_difference'])
def test_setops_disallow_true(method):
    idx1 = pd.MultiIndex.from_product([['a', 'b'], [1, 2]])
    idx2 = pd.MultiIndex.from_product([['b', 'c'], [1, 2]])

    with pytest.raises(ValueError, match="The 'sort' keyword only takes"):
        getattr(idx1, method)(idx2, sort=True)
