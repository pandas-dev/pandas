# -*- coding: utf-8 -*-


from pandas import MultiIndex


def test_sortlevel(_index):
    import random

    tuples = list(_index)
    random.shuffle(tuples)

    index = MultiIndex.from_tuples(tuples)

    sorted_idx, _ = index.sortlevel(0)
    expected = MultiIndex.from_tuples(sorted(tuples))
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(0, ascending=False)
    assert sorted_idx.equals(expected[::-1])

    sorted_idx, _ = index.sortlevel(1)
    by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
    expected = MultiIndex.from_tuples(by1)
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(1, ascending=False)
    assert sorted_idx.equals(expected[::-1])


def test_sortlevel_not_sort_remaining():
    mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
    sorted_idx, _ = mi.sortlevel('A', sort_remaining=False)
    assert sorted_idx.equals(mi)


def test_sortlevel_deterministic():
    tuples = [('bar', 'one'), ('foo', 'two'), ('qux', 'two'),
              ('foo', 'one'), ('baz', 'two'), ('qux', 'one')]

    index = MultiIndex.from_tuples(tuples)

    sorted_idx, _ = index.sortlevel(0)
    expected = MultiIndex.from_tuples(sorted(tuples))
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(0, ascending=False)
    assert sorted_idx.equals(expected[::-1])

    sorted_idx, _ = index.sortlevel(1)
    by1 = sorted(tuples, key=lambda x: (x[1], x[0]))
    expected = MultiIndex.from_tuples(by1)
    assert sorted_idx.equals(expected)

    sorted_idx, _ = index.sortlevel(1, ascending=False)
    assert sorted_idx.equals(expected[::-1])
