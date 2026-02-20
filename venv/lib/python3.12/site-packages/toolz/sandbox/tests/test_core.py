from toolz import curry, unique, first, take
from toolz.sandbox.core import EqualityHashKey, unzip
from itertools import count, repeat

def test_EqualityHashKey_default_key():
    EqualityHashDefault = curry(EqualityHashKey, None)
    L1 = [1]
    L2 = [2]
    data1 = [L1, L1, L2, [], [], [1], [2], {}, ()]
    set1 = set(map(EqualityHashDefault, data1))
    set2 = set(map(EqualityHashDefault, [[], [1], [2], {}, ()]))
    assert set1 == set2
    assert len(set1) == 5

    # Test that ``EqualityHashDefault(item)`` is distinct from ``item``
    T0 = ()
    T1 = (1,)
    data2 = list(map(EqualityHashDefault, [T0, T0, T1, T1, (), (1,)]))
    data2.extend([T0, T1, (), (1,)])
    set3 = set(data2)
    assert set3 == {(), (1,), EqualityHashDefault(()),
                        EqualityHashDefault((1,))}
    assert len(set3) == 4
    assert EqualityHashDefault(()) in set3
    assert EqualityHashDefault((1,)) in set3

    # Miscellaneous
    E1 = EqualityHashDefault(L1)
    E2 = EqualityHashDefault(L2)
    assert str(E1) == '=[1]='
    assert repr(E1) == '=[1]='
    assert E1 != E2
    assert not (E1 == E2)
    assert E1 == EqualityHashDefault(L1)
    assert not (E1 != EqualityHashDefault(L1))
    assert E1 != L1
    assert not (E1 == L1)


def test_EqualityHashKey_callable_key():
    # Common simple hash key functions.
    EqualityHashLen = curry(EqualityHashKey, len)
    EqualityHashType = curry(EqualityHashKey, type)
    EqualityHashId = curry(EqualityHashKey, id)
    EqualityHashFirst = curry(EqualityHashKey, first)
    data1 = [[], [1], (), (1,), {}, {1: 2}]
    data2 = [[1, 2], (1, 2), (1, 3), [1, 3], [2, 1], {1: 2}]
    assert list(unique(data1*3, key=EqualityHashLen)) == data1
    assert list(unique(data2*3, key=EqualityHashLen)) == data2
    assert list(unique(data1*3, key=EqualityHashType)) == data1
    assert list(unique(data2*3, key=EqualityHashType)) == data2
    assert list(unique(data1*3, key=EqualityHashId)) == data1
    assert list(unique(data2*3, key=EqualityHashId)) == data2
    assert list(unique(data2*3, key=EqualityHashFirst)) == data2


def test_EqualityHashKey_index_key():
    d1 = {'firstname': 'Alice', 'age': 21, 'data': {}}
    d2 = {'firstname': 'Alice', 'age': 34, 'data': {}}
    d3a = {'firstname': 'Bob', 'age': 56, 'data': {}}
    d3b = {'firstname': 'Bob', 'age': 56, 'data': {}}
    EqualityHashFirstname = curry(EqualityHashKey, 'firstname')
    assert list(unique(3*[d1, d2, d3a, d3b],
                       key=EqualityHashFirstname)) == [d1, d2, d3a]
    EqualityHashFirstnameAge = curry(EqualityHashKey, ['firstname', 'age'])
    assert list(unique(3*[d1, d2, d3a, d3b],
                       key=EqualityHashFirstnameAge)) == [d1, d2, d3a]
    list1 = [0] * 10
    list2 = [0] * 100
    list3a = [1] * 10
    list3b = [1] * 10
    EqualityHash0 = curry(EqualityHashKey, 0)
    assert list(unique(3*[list1, list2, list3a, list3b],
                       key=EqualityHash0)) == [list1, list2, list3a]


def test_unzip():
    def _to_lists(seq, n=10):
        """iter of iters -> finite list of finite lists
        """
        def initial(s):
            return list(take(n, s))

        return initial(map(initial, seq))

    def _assert_initial_matches(a, b, n=10):
        assert list(take(n, a)) == list(take(n, b))

    # Unzips a simple list correctly
    assert _to_lists(unzip([('a', 1), ('b', 2), ('c', 3)])) \
        == [['a', 'b', 'c'], [1, 2, 3]]

    # Can handle a finite number of infinite iterators (the naive unzip
    # implementation `zip(*args)` implementation fails on this example).
    a, b, c = unzip(zip(count(1), repeat(0), repeat(1)))
    _assert_initial_matches(a, count(1))
    _assert_initial_matches(b, repeat(0))
    _assert_initial_matches(c, repeat(1))

    # Sensibly handles empty input
    assert list(unzip(zip([]))) == []
