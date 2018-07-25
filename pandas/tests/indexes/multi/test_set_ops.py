# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas import MultiIndex, Series


def test_setops_errorcases(idx):
    # # non-iterable input
    cases = [0.5, 'xxx']
    methods = [idx.intersection, idx.union, idx.difference,
               idx.symmetric_difference]

    for method in methods:
        for case in cases:
            tm.assert_raises_regex(TypeError,
                                   "Input must be Index "
                                   "or array-like",
                                   method, case)


def test_intersection_base(idx):
    first = idx[:5]
    second = idx[:3]
    intersect = first.intersection(second)

    assert tm.equalContents(intersect, second)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.intersection(case)
        assert tm.equalContents(result, second)

    msg = "other must be a MultiIndex or a list of tuples"
    with tm.assert_raises_regex(TypeError, msg):
        result = first.intersection([1, 2, 3])


def test_union_base(idx):
    first = idx[3:]
    second = idx[:5]
    everything = idx
    union = first.union(second)
    assert tm.equalContents(union, everything)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.union(case)
        assert tm.equalContents(result, everything)

    msg = "other must be a MultiIndex or a list of tuples"
    with tm.assert_raises_regex(TypeError, msg):
        result = first.union([1, 2, 3])


def test_difference_base(idx):
    first = idx[2:]
    second = idx[:4]
    answer = idx[4:]
    result = first.difference(second)

    assert tm.equalContents(result, answer)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.difference(case)
        assert tm.equalContents(result, answer)

    msg = "other must be a MultiIndex or a list of tuples"
    with tm.assert_raises_regex(TypeError, msg):
        result = first.difference([1, 2, 3])


def test_symmetric_difference(idx):
    first = idx[1:]
    second = idx[:-1]
    answer = idx[[0, -1]]
    result = first.symmetric_difference(second)
    assert tm.equalContents(result, answer)

    # GH 10149
    cases = [klass(second.values)
             for klass in [np.array, Series, list]]
    for case in cases:
        result = first.symmetric_difference(case)
        assert tm.equalContents(result, answer)

    msg = "other must be a MultiIndex or a list of tuples"
    with tm.assert_raises_regex(TypeError, msg):
        first.symmetric_difference([1, 2, 3])


def test_empty(idx):
    # GH 15270
    assert not idx.empty
    assert idx[:0].empty


def test_difference(idx):

    first = idx
    result = first.difference(idx[-3:])
    expected = MultiIndex.from_tuples(sorted(idx[:-3].values),
                                      sortorder=0,
                                      names=idx.names)

    assert isinstance(result, MultiIndex)
    assert result.equals(expected)
    assert result.names == idx.names

    # empty difference: reflexive
    result = idx.difference(idx)
    expected = idx[:0]
    assert result.equals(expected)
    assert result.names == idx.names

    # empty difference: superset
    result = idx[-3:].difference(idx)
    expected = idx[:0]
    assert result.equals(expected)
    assert result.names == idx.names

    # empty difference: degenerate
    result = idx[:0].difference(idx)
    expected = idx[:0]
    assert result.equals(expected)
    assert result.names == idx.names

    # names not the same
    chunklet = idx[-3:]
    chunklet.names = ['foo', 'baz']
    result = first.difference(chunklet)
    assert result.names == (None, None)

    # empty, but non-equal
    result = idx.difference(idx.sortlevel(1)[0])
    assert len(result) == 0

    # raise Exception called with non-MultiIndex
    result = first.difference(first.values)
    assert result.equals(first[:0])

    # name from empty array
    result = first.difference([])
    assert first.equals(result)
    assert first.names == result.names

    # name from non-empty array
    result = first.difference([('foo', 'one')])
    expected = pd.MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'), (
        'foo', 'two'), ('qux', 'one'), ('qux', 'two')])
    expected.names = first.names
    assert first.names == result.names
    tm.assert_raises_regex(TypeError, "other must be a MultiIndex "
                           "or a list of tuples",
                           first.difference, [1, 2, 3, 4, 5])


def test_union(idx):
    piece1 = idx[:5][::-1]
    piece2 = idx[3:]

    the_union = piece1 | piece2

    tups = sorted(idx.values)
    expected = MultiIndex.from_tuples(tups)

    assert the_union.equals(expected)

    # corner case, pass self or empty thing:
    the_union = idx.union(idx)
    assert the_union is idx

    the_union = idx.union(idx[:0])
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


def test_intersection(idx):
    piece1 = idx[:5][::-1]
    piece2 = idx[3:]

    the_int = piece1 & piece2
    tups = sorted(idx[3:5].values)
    expected = MultiIndex.from_tuples(tups)
    assert the_int.equals(expected)

    # corner case, pass self
    the_int = idx.intersection(idx)
    assert the_int is idx

    # empty intersection: disjoint
    empty = idx[:2] & idx[2:]
    expected = idx[:0]
    assert empty.equals(expected)

    # can't do in python 3
    # tuples = _index.values
    # result = _index & tuples
    # assert result.equals(tuples)
