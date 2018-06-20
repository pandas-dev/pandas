from pandas import (CategoricalIndex, DatetimeIndex, Float64Index, Index,
                    Int64Index, IntervalIndex, MultiIndex, PeriodIndex,
                    RangeIndex, Series, TimedeltaIndex, UInt64Index, compat,
                    isna)

import pandas.util.testing as tm
import numpy as np
import pandas as pd


def test_setops_errorcases(named_index):
    for name, idx in compat.iteritems(named_index):
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


def test_intersection_base(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[:5]
        second = idx[:3]
        intersect = first.intersection(second)

        if isinstance(idx, CategoricalIndex):
            pass
        else:
            assert tm.equalContents(intersect, second)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.intersection(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            else:
                result = first.intersection(case)
                assert tm.equalContents(result, second)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                result = first.intersection([1, 2, 3])


def test_union_base(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[3:]
        second = idx[:5]
        everything = idx
        union = first.union(second)
        assert tm.equalContents(union, everything)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.union(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            else:
                result = first.union(case)
                assert tm.equalContents(result, everything)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                result = first.union([1, 2, 3])


def test_difference_base(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[2:]
        second = idx[:4]
        answer = idx[4:]
        result = first.difference(second)

        if isinstance(idx, CategoricalIndex):
            pass
        else:
            assert tm.equalContents(result, answer)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.difference(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            elif isinstance(idx, (DatetimeIndex, TimedeltaIndex)):
                assert result.__class__ == answer.__class__
                tm.assert_numpy_array_equal(result.sort_values().asi8,
                                            answer.sort_values().asi8)
            else:
                result = first.difference(case)
                assert tm.equalContents(result, answer)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                result = first.difference([1, 2, 3])


def test_symmetric_difference(named_index):
    for name, idx in compat.iteritems(named_index):
        first = idx[1:]
        second = idx[:-1]
        if isinstance(idx, CategoricalIndex):
            pass
        else:
            answer = idx[[0, -1]]
            result = first.symmetric_difference(second)
            assert tm.equalContents(result, answer)

        # GH 10149
        cases = [klass(second.values)
                 for klass in [np.array, Series, list]]
        for case in cases:
            if isinstance(idx, PeriodIndex):
                msg = "can only call with other PeriodIndex-ed objects"
                with tm.assert_raises_regex(ValueError, msg):
                    result = first.symmetric_difference(case)
            elif isinstance(idx, CategoricalIndex):
                pass
            else:
                result = first.symmetric_difference(case)
                assert tm.equalContents(result, answer)

        if isinstance(idx, MultiIndex):
            msg = "other must be a MultiIndex or a list of tuples"
            with tm.assert_raises_regex(TypeError, msg):
                first.symmetric_difference([1, 2, 3])


def test_empty(_index):
    # GH 15270
    index = _index
    assert not index.empty
    assert index[:0].empty



def test_unique_na():
    idx = pd.Index([2, np.nan, 2, 1], name='my_index')
    expected = pd.Index([2, np.nan, 1], name='my_index')
    result = idx.unique()
    tm.assert_index_equal(result, expected)


def test_difference(_index):

    first = _index
    result = first.difference(_index[-3:])
    expected = MultiIndex.from_tuples(sorted(_index[:-3].values),
                                      sortorder=0,
                                      names=_index.names)

    assert isinstance(result, MultiIndex)
    assert result.equals(expected)
    assert result.names == _index.names

    # empty difference: reflexive
    result = _index.difference(_index)
    expected = _index[:0]
    assert result.equals(expected)
    assert result.names == _index.names

    # empty difference: superset
    result = _index[-3:].difference(_index)
    expected = _index[:0]
    assert result.equals(expected)
    assert result.names == _index.names

    # empty difference: degenerate
    result = _index[:0].difference(_index)
    expected = _index[:0]
    assert result.equals(expected)
    assert result.names == _index.names

    # names not the same
    chunklet = _index[-3:]
    chunklet.names = ['foo', 'baz']
    result = first.difference(chunklet)
    assert result.names == (None, None)

    # empty, but non-equal
    result = _index.difference(_index.sortlevel(1)[0])
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


def test_union(_index):
    piece1 = _index[:5][::-1]
    piece2 = _index[3:]

    the_union = piece1 | piece2

    tups = sorted(_index.values)
    expected = MultiIndex.from_tuples(tups)

    assert the_union.equals(expected)

    # corner case, pass self or empty thing:
    the_union = _index.union(_index)
    assert the_union is _index

    the_union = _index.union(_index[:0])
    assert the_union is _index

    # won't work in python 3
    # tuples = _index.values
    # result = _index[:4] | tuples[4:]
    # assert result.equals(tuples)

    # not valid for python 3
    # def test_union_with_regular_index(self):
    #     other = Index(['A', 'B', 'C'])

    #     result = other.union(_index)
    #     assert ('foo', 'one') in result
    #     assert 'B' in result

    #     result2 = _index.union(other)
    #     assert result.equals(result2)


def test_intersection(_index):
    piece1 = _index[:5][::-1]
    piece2 = _index[3:]

    the_int = piece1 & piece2
    tups = sorted(_index[3:5].values)
    expected = MultiIndex.from_tuples(tups)
    assert the_int.equals(expected)

    # corner case, pass self
    the_int = _index.intersection(_index)
    assert the_int is _index

    # empty intersection: disjoint
    empty = _index[:2] & _index[2:]
    expected = _index[:0]
    assert empty.equals(expected)

    # can't do in python 3
    # tuples = _index.values
    # result = _index & tuples
    # assert result.equals(tuples)
