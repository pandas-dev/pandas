from __future__ import annotations

import gc
import math
import os
import random
import warnings
import weakref
from bz2 import BZ2File
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from gzip import GzipFile
from itertools import repeat

import partd
import pytest
from tlz import groupby, identity, join, merge, pluck, unique, valmap

import dask
import dask.bag as db
from dask._task_spec import List, Task
from dask.bag.core import (
    Bag,
    collect,
    from_delayed,
    lazify,
    lazify_task,
    partition,
    reduceby,
    reify,
    total_mem_usage,
)
from dask.bag.utils import assert_eq
from dask.blockwise import Blockwise
from dask.delayed import Delayed
from dask.typing import Graph
from dask.utils import filetexts, tmpdir, tmpfile
from dask.utils_test import add, hlg_layer, inc

dsk: Graph = {("x", 0): (range, 5), ("x", 1): (range, 5), ("x", 2): (range, 5)}

L = list(range(5)) * 3
b = Bag(dsk, "x", 3)


def iseven(x):
    return x % 2 == 0


def isodd(x):
    return x % 2 == 1


def test_Bag():
    assert b.name == "x"
    assert b.npartitions == 3


def test_keys():
    assert b.__dask_keys__() == sorted(dsk.keys())


def test_bag_groupby_pure_hash():
    # https://github.com/dask/dask/issues/6640
    result = b.groupby(iseven).compute()
    assert result == [(True, [0, 2, 4] * 3), (False, [1, 3] * 3)]


def test_bag_groupby_normal_hash():
    # https://github.com/dask/dask/issues/6640
    result = b.groupby(lambda x: "even" if iseven(x) else "odd").compute()
    assert len(result) == 2
    assert ("odd", [1, 3] * 3) in result
    assert ("even", [0, 2, 4] * 3) in result


@pytest.mark.parametrize("shuffle", ["disk", "tasks"])
@pytest.mark.parametrize("scheduler", ["synchronous", "processes"])
def test_bag_groupby_none(shuffle, scheduler):
    with dask.config.set(scheduler=scheduler):
        seq = [(None, i) for i in range(50)]
        b = db.from_sequence(seq).groupby(lambda x: x[0], shuffle=shuffle)
        result = b.compute()
        assert len(result) == 1


@dataclass(frozen=True)
class Key:
    foo: int
    bar: int | None = None


@pytest.mark.parametrize(
    "key",
    # if a value for `bar` is not explicitly passed, Key.bar will default to `None`,
    # thereby introducing the risk of inter-process inconsistency for the value returned by
    # built-in `hash` (due to lack of deterministic hashing for `None` prior to python 3.12).
    # without https://github.com/dask/dask/pull/10734, this results in failures for this test.
    [Key(foo=1), Key(foo=1, bar=2)],
    ids=["none_field", "no_none_fields"],
)
@pytest.mark.parametrize("shuffle", ["disk", "tasks"])
@pytest.mark.parametrize("scheduler", ["synchronous", "processes"])
def test_bag_groupby_dataclass(key, shuffle, scheduler):
    seq = [(key, i) for i in range(50)]
    b = db.from_sequence(seq).groupby(lambda x: x[0], shuffle=shuffle)
    with dask.config.set(scheduler=scheduler):
        result = b.compute()
    assert len(result) == 1


def test_bag_map():
    b = db.from_sequence(range(100), npartitions=10)
    b2 = db.from_sequence(range(100, 200), npartitions=10)
    x = b.compute()
    x2 = b2.compute()

    def myadd(a=1, b=2, c=3):
        return a + b + c

    assert_eq(db.map(myadd, b), list(map(myadd, x)))
    assert_eq(db.map(myadd, a=b), list(map(myadd, x)))
    assert_eq(db.map(myadd, b, b2), list(map(myadd, x, x2)))
    assert_eq(db.map(myadd, b, 10), [myadd(i, 10) for i in x])
    assert_eq(db.map(myadd, 10, b=b), [myadd(10, b=i) for i in x])

    sol = [myadd(i, b=j, c=100) for (i, j) in zip(x, x2)]
    assert_eq(db.map(myadd, b, b=b2, c=100), sol)

    sol = [myadd(i, c=100) for (i, j) in zip(x, x2)]
    assert_eq(db.map(myadd, b, c=100), sol)

    x_sum = sum(x)
    sol = [myadd(x_sum, b=i, c=100) for i in x2]
    assert_eq(db.map(myadd, b.sum(), b=b2, c=100), sol)

    sol = [myadd(i, b=x_sum, c=100) for i in x2]
    assert_eq(db.map(myadd, b2, b.sum(), c=100), sol)

    sol = [myadd(a=100, b=x_sum, c=i) for i in x2]
    assert_eq(db.map(myadd, a=100, b=b.sum(), c=b2), sol)

    a = dask.delayed(10)
    assert_eq(db.map(myadd, b, a), [myadd(i, 10) for i in x])
    assert_eq(db.map(myadd, b, b=a), [myadd(i, b=10) for i in x])

    # Mispatched npartitions
    fewer_parts = db.from_sequence(range(100), npartitions=5)
    with pytest.raises(ValueError):
        db.map(myadd, b, fewer_parts)

    # No bags
    with pytest.raises(ValueError):
        db.map(myadd, b.sum(), 1, 2)

    # Unequal partitioning
    unequal = db.from_sequence(range(110), npartitions=10)
    with pytest.raises(ValueError):
        db.map(myadd, b, unequal, c=b2).compute()
    with pytest.raises(ValueError):
        db.map(myadd, b, b=unequal, c=b2).compute()


def test_map_method():
    b = db.from_sequence(range(100), npartitions=10)
    b2 = db.from_sequence(range(100, 200), npartitions=10)
    x = b.compute()
    x2 = b2.compute()

    def myadd(a, b=2, c=3):
        return a + b + c

    assert b.map(myadd).compute() == list(map(myadd, x))
    assert b.map(myadd, b2).compute() == list(map(myadd, x, x2))
    assert b.map(myadd, 10).compute() == [myadd(i, 10) for i in x]
    assert b.map(myadd, b=10).compute() == [myadd(i, b=10) for i in x]
    assert b.map(myadd, b2, c=10).compute() == [
        myadd(i, j, 10) for (i, j) in zip(x, x2)
    ]
    x_sum = sum(x)
    assert b.map(myadd, b.sum(), c=10).compute() == [myadd(i, x_sum, 10) for i in x]


def test_starmap():
    data = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10)]
    b = db.from_sequence(data, npartitions=2)

    def myadd(a, b, c=0):
        return a + b + c

    assert b.starmap(myadd).compute() == [myadd(*a) for a in data]
    assert b.starmap(myadd, c=10).compute() == [myadd(*a, c=10) for a in data]
    max_second = b.pluck(1).max()
    assert b.starmap(myadd, c=max_second).compute() == [
        myadd(*a, c=max_second.compute()) for a in data
    ]
    c = dask.delayed(10)
    assert b.starmap(myadd, c=c).compute() == [myadd(*a, c=10) for a in data]


def test_filter():
    c = b.filter(iseven)
    expected = merge(
        dsk,
        {
            (c.name, i): (reify, (filter, iseven, (b.name, i)))
            for i in range(b.npartitions)
        },
    )
    assert c.dask == expected
    assert c.name == b.filter(iseven).name


def test_remove():
    f = lambda x: x % 2 == 0
    c = b.remove(f)
    assert list(c) == [1, 3] * 3
    assert c.name == b.remove(f).name


def test_iter():
    assert sorted(b) == sorted(L)
    assert sorted(b.map(inc)) == sorted(list(range(1, 6)) * 3)


@pytest.mark.parametrize("func", [str, repr])
def test_repr(func):
    assert str(b.npartitions) in func(b)
    assert b.name[:5] in func(b)

    assert "from_sequence" in func(db.from_sequence(range(5)))


def test_pluck():
    d = {("x", 0): [(1, 10), (2, 20)], ("x", 1): [(3, 30), (4, 40)]}
    b = Bag(d, "x", 2)
    assert set(b.pluck(0)) == {1, 2, 3, 4}
    assert set(b.pluck(1)) == {10, 20, 30, 40}
    assert set(b.pluck([1, 0])) == {(10, 1), (20, 2), (30, 3), (40, 4)}
    assert b.pluck([1, 0]).name == b.pluck([1, 0]).name


def test_pluck_with_default():
    b = db.from_sequence(["Hello", "", "World"])
    pytest.raises(IndexError, lambda: list(b.pluck(0)))
    assert list(b.pluck(0, None)) == ["H", None, "W"]
    assert b.pluck(0, None).name == b.pluck(0, None).name
    assert b.pluck(0).name != b.pluck(0, None).name


def test_unzip():
    b = db.from_sequence(range(100)).map(lambda x: (x, x + 1, x + 2))
    one, two, three = b.unzip(3)
    assert list(one) == list(range(100))
    assert list(three) == [i + 2 for i in range(100)]
    assert one.name == b.unzip(3)[0].name
    assert one.name != two.name


def test_fold():
    c = b.fold(add)
    assert c.compute() == sum(L)
    assert c.key == b.fold(add).key

    c2 = b.fold(add, initial=10)
    assert c2.key != c.key
    assert c2.compute() == sum(L) + 10 * b.npartitions
    assert c2.key == b.fold(add, initial=10).key

    c = db.from_sequence(range(5), npartitions=3)

    def binop(acc, x):
        acc = acc.copy()
        acc.add(x)
        return acc

    d = c.fold(binop, set.union, initial=set())
    assert d.compute() == set(c)
    assert d.key == c.fold(binop, set.union, initial=set()).key

    d = db.from_sequence("hello")
    assert set(d.fold(lambda a, b: "".join([a, b]), initial="").compute()) == set(
        "hello"
    )

    e = db.from_sequence([[1], [2], [3]], npartitions=2)
    assert set(e.fold(add, initial=[]).compute(scheduler="sync")) == {1, 2, 3}


def test_fold_bag():
    def binop(tot, x):
        tot.add(x)
        return tot

    c = b.fold(binop, combine=set.union, initial=set(), out_type=Bag)
    assert isinstance(c, Bag)
    assert_eq(c, list(set(range(5))))


def test_distinct():
    assert sorted(b.distinct()) == [0, 1, 2, 3, 4]
    assert b.distinct().name == b.distinct().name
    assert "distinct" in b.distinct().name
    assert b.distinct().count().compute() == 5
    bag = db.from_sequence([0] * 50, npartitions=50)
    assert bag.filter(None).distinct().compute() == []


def test_distinct_with_key():
    seq = [{"a": i} for i in [0, 1, 2, 1, 2, 3, 2, 3, 4, 5]]
    bag = db.from_sequence(seq, npartitions=3)
    expected = list(unique(seq, key=lambda x: x["a"]))
    assert_eq(bag.distinct(key="a"), expected)
    assert_eq(bag.distinct(key=lambda x: x["a"]), expected)


def test_frequencies():
    c = b.frequencies()
    assert dict(c) == {0: 3, 1: 3, 2: 3, 3: 3, 4: 3}
    c2 = b.frequencies(split_every=2)
    assert dict(c2) == {0: 3, 1: 3, 2: 3, 3: 3, 4: 3}
    assert c.name == b.frequencies().name
    assert c.name != c2.name
    assert c2.name == b.frequencies(split_every=2).name
    # test bag with empty partitions
    b2 = db.from_sequence(range(20), partition_size=2)
    b2 = b2.filter(lambda x: x < 10)
    d = b2.frequencies()
    assert dict(d) == dict(zip(range(10), [1] * 10))
    bag = db.from_sequence([0, 0, 0, 0], npartitions=4)
    bag2 = bag.filter(None).frequencies(split_every=2)
    assert_eq(bag2, [])


def test_frequencies_sorted():
    b = db.from_sequence(["a", "b", "b", "b", "c", "c"])
    assert list(b.frequencies(sort=True).compute()) == [("b", 3), ("c", 2), ("a", 1)]


def test_topk():
    assert list(b.topk(4)) == [4, 4, 4, 3]
    c = b.topk(4, key=lambda x: -x)
    assert list(c) == [0, 0, 0, 1]
    c2 = b.topk(4, key=lambda x: -x, split_every=2)
    assert list(c2) == [0, 0, 0, 1]
    assert c.name != c2.name
    assert b.topk(4).name == b.topk(4).name


@pytest.mark.parametrize("npartitions", [1, 2])
def test_topk_with_non_callable_key(npartitions):
    b = db.from_sequence([(1, 10), (2, 9), (3, 8)], npartitions=npartitions)
    assert list(b.topk(2, key=1)) == [(1, 10), (2, 9)]
    assert list(b.topk(2, key=0)) == [(3, 8), (2, 9)]
    assert b.topk(2, key=1).name == b.topk(2, key=1).name


def test_topk_with_multiarg_lambda():
    b = db.from_sequence([(1, 10), (2, 9), (3, 8)], npartitions=2)
    assert list(b.topk(2, key=lambda a, b: b)) == [(1, 10), (2, 9)]


def test_lambdas():
    assert list(b.map(lambda x: x + 1)) == list(b.map(inc))


def test_reductions():
    assert int(b.count()) == 15
    assert int(b.sum()) == 30
    assert int(b.max()) == 4
    assert int(b.min()) == 0
    assert b.any().compute() is True
    assert b.all().compute() is False
    assert b.all().key == b.all().key
    assert b.all().key != b.any().key


def test_reduction_names():
    assert b.sum().name.startswith("sum")
    assert b.reduction(sum, sum).name.startswith("sum")
    assert any(
        isinstance(k, str) and k.startswith("max") for k in b.reduction(sum, max).dask
    )
    assert b.reduction(sum, sum, name="foo").name.startswith("foo")


def test_tree_reductions():
    b = db.from_sequence(range(12))
    c = b.reduction(sum, sum, split_every=2)
    d = b.reduction(sum, sum, split_every=6)
    e = b.reduction(sum, sum, split_every=5)

    assert c.compute() == d.compute() == e.compute()

    assert len(c.dask) > len(d.dask)

    c = b.sum(split_every=2)
    d = b.sum(split_every=5)

    assert c.compute() == d.compute()
    assert len(c.dask) > len(d.dask)

    assert c.key != d.key
    assert c.key == b.sum(split_every=2).key
    assert c.key != b.sum().key


@pytest.mark.parametrize("npartitions", [1, 3, 4])
def test_aggregation(npartitions):
    L = list(range(15))
    b = db.range(15, npartitions=npartitions)
    assert_eq(b.mean(), sum(L) / len(L))
    assert_eq(b.sum(), sum(L))
    assert_eq(b.count(), len(L))


@pytest.mark.parametrize("npartitions", [1, 10])
def test_non_splittable_reductions(npartitions):
    np = pytest.importorskip("numpy")
    data = list(range(100))
    c = db.from_sequence(data, npartitions=npartitions)

    assert_eq(c.mean(), np.mean(data))
    assert_eq(c.std(), np.std(data))


def test_std():
    assert_eq(b.std(), math.sqrt(2.0))
    assert float(b.std()) == math.sqrt(2.0)


def test_var():
    assert_eq(b.var(), 2.0)
    assert float(b.var()) == 2.0


@pytest.mark.parametrize(
    "transform", [identity, dask.delayed, lambda x: db.from_sequence(x, npartitions=1)]
)
def test_join(transform):
    other = transform([1, 2, 3])
    c = b.join(other, on_self=isodd, on_other=iseven)
    assert_eq(c, list(join(iseven, [1, 2, 3], isodd, list(b))))
    assert_eq(b.join(other, isodd), list(join(isodd, [1, 2, 3], isodd, list(b))))
    assert c.name == b.join(other, on_self=isodd, on_other=iseven).name


def test_foldby():
    c = b.foldby(iseven, add, 0, add, 0)
    assert (reduceby, iseven, add, (b.name, 0), 0) in list(c.dask.values())
    assert set(c) == set(reduceby(iseven, lambda acc, x: acc + x, L, 0).items())
    assert c.name == b.foldby(iseven, add, 0, add, 0).name

    c = b.foldby(iseven, lambda acc, x: acc + x)
    assert set(c) == set(reduceby(iseven, lambda acc, x: acc + x, L, 0).items())


def test_foldby_tree_reduction():
    dsk = list()
    for n in [1, 7, 32]:
        b = db.from_sequence(range(100), npartitions=n)
        c = b.foldby(iseven, add)
        dsk.extend([c])
        for m in [False, None, 2, 3]:
            d = b.foldby(iseven, add, split_every=m)
            e = b.foldby(iseven, add, 0, split_every=m)
            f = b.foldby(iseven, add, 0, add, split_every=m)
            g = b.foldby(iseven, add, 0, add, 0, split_every=m)
            dsk.extend([d, e, f, g])
    results = dask.compute(dsk)
    first = results[0]
    assert all([r == first for r in results])


def test_map_partitions():
    assert list(b.map_partitions(len)) == [5, 5, 5]
    assert b.map_partitions(len).name == b.map_partitions(len).name
    assert b.map_partitions(lambda a: len(a) + 1).name != b.map_partitions(len).name


def test_map_partitions_args_kwargs():
    x = [random.randint(-100, 100) for i in range(100)]
    y = [random.randint(-100, 100) for i in range(100)]

    dx = db.from_sequence(x, npartitions=10)
    dy = db.from_sequence(y, npartitions=10)

    def maximum(x, y=0):
        y = repeat(y) if isinstance(y, int) else y
        return [max(a, b) for (a, b) in zip(x, y)]

    sol = maximum(x, y=10)
    assert_eq(db.map_partitions(maximum, dx, y=10), sol)
    assert_eq(dx.map_partitions(maximum, y=10), sol)
    assert_eq(dx.map_partitions(maximum, 10), sol)

    sol = maximum(x, y)
    assert_eq(db.map_partitions(maximum, dx, dy), sol)
    assert_eq(dx.map_partitions(maximum, y=dy), sol)
    assert_eq(dx.map_partitions(maximum, dy), sol)

    dy_mean = dy.mean().apply(int)
    sol = maximum(x, int(sum(y) / len(y)))
    assert_eq(dx.map_partitions(maximum, y=dy_mean), sol)
    assert_eq(dx.map_partitions(maximum, dy_mean), sol)

    dy_mean = dask.delayed(dy_mean)
    assert_eq(dx.map_partitions(maximum, y=dy_mean), sol)
    assert_eq(dx.map_partitions(maximum, dy_mean), sol)


def test_map_partitions_blockwise():
    # Check that the `token` argument works,
    # and that `map_partitions`` is using `Blockwise`.
    layer = hlg_layer(
        b.map_partitions(lambda x: x, token="test-string").dask, "test-string"
    )
    assert layer
    assert isinstance(layer, Blockwise)


def test_random_sample_size():
    """
    Number of randomly sampled elements are in the expected range.
    """
    a = db.from_sequence(range(1000), npartitions=5)
    # we expect a size of approx. 100, but leave large margins to avoid
    # random failures
    assert 10 < len(list(a.random_sample(0.1, 42))) < 300


def test_random_sample_prob_range():
    """
    Specifying probabilities outside the range [0, 1] raises ValueError.
    """
    a = db.from_sequence(range(50), npartitions=5)
    with pytest.raises(ValueError):
        a.random_sample(-1)
    with pytest.raises(ValueError):
        a.random_sample(1.1)


def test_random_sample_repeated_computation():
    """
    Repeated computation of a defined random sampling operation
    generates identical results.
    """
    a = db.from_sequence(range(50), npartitions=5)
    b = a.random_sample(0.2)
    assert list(b) == list(b)  # computation happens here


def test_random_sample_different_definitions():
    """
    Repeatedly defining a random sampling operation yields different results
    upon computation if no random seed is specified.
    """
    a = db.from_sequence(range(50), npartitions=5)
    assert list(a.random_sample(0.5)) != list(a.random_sample(0.5))
    assert a.random_sample(0.5).name != a.random_sample(0.5).name


def test_random_sample_random_state():
    """
    Sampling with fixed random seed generates identical results.
    """
    a = db.from_sequence(range(50), npartitions=5)
    b = a.random_sample(0.5, 1234)
    c = a.random_sample(0.5, 1234)
    assert list(b) == list(c)


def test_lazify_task():
    task = Task(None, sum, Task(None, reify, Task(None, map, inc, [1, 2, 3])))
    assert lazify_task(task) == Task(None, sum, Task(None, map, inc, [1, 2, 3]))

    task = Task(None, reify, Task(None, map, inc, [1, 2, 3]))
    assert lazify_task(task) == task

    task = Task(
        None,
        reify,
        Task(None, map, inc, Task(None, reify, Task(None, filter, iseven, "y"))),
    )
    expected = Task(None, reify, Task(None, map, inc, Task(None, filter, iseven, "y")))
    assert lazify_task(task) == expected

    task = Task(
        None,
        sum,
        List(
            Task(
                None,
                reify,
                Task(None, map, inc, Task(None, reify, Task(None, filter, iseven, 1))),
            ),
            Task(
                None,
                reify,
                Task(None, map, inc, Task(None, reify, Task(None, filter, iseven, 2))),
            ),
        ),
    )
    assert lazify_task(task) == Task(
        None,
        sum,
        List(
            Task(None, map, inc, Task(None, filter, iseven, 1)),
            Task(None, map, inc, Task(None, filter, iseven, 2)),
        ),
    )


def test_lazify_task_legacy():
    task = (sum, (reify, (map, inc, [1, 2, 3])))
    assert lazify_task(task) == (sum, (map, inc, [1, 2, 3]))

    task = (reify, (map, inc, [1, 2, 3]))
    assert lazify_task(task) == task

    a = (reify, (map, inc, (reify, (filter, iseven, "y"))))
    b = (reify, (map, inc, (filter, iseven, "y")))
    assert lazify_task(a) == b


f = lambda x: x


def test_lazify():
    a = {
        "x": (reify, (map, inc, (reify, (filter, iseven, "y")))),
        "a": (f, "x"),
        "b": (f, "x"),
    }
    b = {"x": (reify, (map, inc, (filter, iseven, "y"))), "a": (f, "x"), "b": (f, "x")}
    assert lazify(a) == b


def test_take():
    assert list(b.take(2)) == [0, 1]
    assert b.take(2) == (0, 1)
    assert isinstance(b.take(2, compute=False), Bag)


def test_take_npartitions():
    assert list(b.take(6, npartitions=2)) == [0, 1, 2, 3, 4, 0]
    assert b.take(6, npartitions=-1) == (0, 1, 2, 3, 4, 0)
    assert b.take(3, npartitions=-1) == (0, 1, 2)
    with pytest.raises(ValueError):
        b.take(1, npartitions=5)


def test_take_npartitions_warn():
    # Use single-threaded scheduler so warnings are properly captured in the
    # same process
    with dask.config.set(scheduler="sync"):
        with pytest.warns(UserWarning):
            b.take(100)

        with pytest.warns(UserWarning):
            b.take(7)

        with warnings.catch_warnings(record=True) as record:
            b.take(7, npartitions=2)
            b.take(7, warn=False)
        assert not record


def test_map_is_lazy():
    assert isinstance(map(lambda x: x, [1, 2, 3]), Iterator)


def test_can_use_dict_to_make_concrete():
    assert isinstance(dict(b.frequencies()), dict)


@pytest.mark.slow
@pytest.mark.network
@pytest.mark.skip(reason="Hangs")
def test_from_url():
    a = db.from_url(["http://google.com", "http://github.com"])
    assert a.npartitions == 2

    b = db.from_url("http://raw.githubusercontent.com/dask/dask/main/README.rst")
    assert b.npartitions == 1
    assert b"Dask\n" in b.take(10)


def test_read_text():
    with filetexts({"a1.log": "A\nB", "a2.log": "C\nD"}) as fns:
        assert {line.strip() for line in db.read_text(fns)} == set("ABCD")
        assert {line.strip() for line in db.read_text("a*.log")} == set("ABCD")

    pytest.raises(ValueError, lambda: db.read_text("non-existent-*-path"))


def test_read_text_large():
    with tmpfile() as fn:
        with open(fn, "wb") as f:
            f.write(("Hello, world!" + os.linesep).encode() * 100)
        b = db.read_text(fn, blocksize=100)
        c = db.read_text(fn)
        assert len(b.dask) > 5
        assert list(map(str, b.str.strip())) == list(map(str, c.str.strip()))

        d = db.read_text([fn], blocksize=100)
        assert list(b) == list(d)


def test_read_text_encoding():
    with tmpfile() as fn:
        with open(fn, "wb") as f:
            f.write(("你好！" + os.linesep).encode("gb18030") * 100)
        b = db.read_text(fn, blocksize=100, encoding="gb18030")
        c = db.read_text(fn, encoding="gb18030")
        assert len(b.dask) > 5
        b_enc = b.str.strip().map(lambda x: x.encode("utf-8"))
        c_enc = c.str.strip().map(lambda x: x.encode("utf-8"))
        assert list(b_enc) == list(c_enc)

        d = db.read_text([fn], blocksize=100, encoding="gb18030")
        assert list(b) == list(d)


def test_read_text_large_gzip():
    with tmpfile("gz") as fn:
        data = b"Hello, world!\n" * 100
        f = GzipFile(fn, "wb")
        f.write(data)
        f.close()

        with pytest.raises(ValueError):
            # not allowed blocks when compressed
            db.read_text(fn, blocksize=50, linedelimiter="\n")

        c = db.read_text(fn, blocksize=None)
        assert c.npartitions == 1
        assert "".join(c.compute()) == data.decode()


@pytest.mark.xfail(reason="https://github.com/dask/dask/issues/6914")
@pytest.mark.slow
@pytest.mark.network
def test_from_s3():
    # note we don't test connection modes with aws_access_key and
    # aws_secret_key because these are not on travis-ci
    pytest.importorskip("s3fs")

    five_tips = (
        "total_bill,tip,sex,smoker,day,time,size\n",
        "16.99,1.01,Female,No,Sun,Dinner,2\n",
        "10.34,1.66,Male,No,Sun,Dinner,3\n",
        "21.01,3.5,Male,No,Sun,Dinner,3\n",
        "23.68,3.31,Male,No,Sun,Dinner,2\n",
    )

    # test compressed data
    e = db.read_text("s3://tip-data/t*.gz", storage_options=dict(anon=True))
    assert e.take(5) == five_tips

    # test multiple keys in bucket
    c = db.read_text(
        ["s3://tip-data/tips.gz", "s3://tip-data/tips.json", "s3://tip-data/tips.csv"],
        storage_options=dict(anon=True),
    )
    assert c.npartitions == 3


def test_from_sequence():
    b = db.from_sequence([1, 2, 3, 4, 5], npartitions=3)
    assert len(b.dask) == 3
    assert set(b) == {1, 2, 3, 4, 5}


def test_from_long_sequence():
    L = list(range(1001))
    b = db.from_sequence(L)
    assert set(b) == set(L)


def test_from_empty_sequence():
    pytest.importorskip("pandas")
    pytest.importorskip("dask.dataframe")
    b = db.from_sequence([])
    assert b.npartitions == 1
    df = b.to_dataframe(meta={"a": "int"}).compute()
    assert df.empty, "DataFrame is not empty"


def test_product():
    b2 = b.product(b)
    assert b2.npartitions == b.npartitions**2
    assert set(b2) == {(i, j) for i in L for j in L}

    x = db.from_sequence([1, 2, 3, 4])
    y = db.from_sequence([10, 20, 30])
    z = x.product(y)
    assert set(z) == {(i, j) for i in [1, 2, 3, 4] for j in [10, 20, 30]}

    assert z.name != b2.name
    assert z.name == x.product(y).name


def test_partition_collect():
    with partd.Pickle() as p:
        partition(identity, range(6), 3, p)
        for i in range(3):
            assert p.get(i)
            assert sorted(collect(identity, i, p, "")) == [(j, [j]) for j in p.get(i)]


def test_groupby():
    c = b.groupby(identity)
    result = dict(c)
    assert result == {
        0: [0, 0, 0],
        1: [1, 1, 1],
        2: [2, 2, 2],
        3: [3, 3, 3],
        4: [4, 4, 4],
    }
    assert c.npartitions == b.npartitions
    assert c.name == b.groupby(identity).name
    assert c.name != b.groupby(lambda x: x + 1).name


def test_groupby_with_indexer():
    b = db.from_sequence([[1, 2, 3], [1, 4, 9], [2, 3, 4]])
    result = dict(b.groupby(0))
    assert valmap(sorted, result) == {1: [[1, 2, 3], [1, 4, 9]], 2: [[2, 3, 4]]}


def test_groupby_with_npartitions_changed():
    result = b.groupby(lambda x: x, npartitions=1)
    result2 = dict(result)
    assert result2 == {
        0: [0, 0, 0],
        1: [1, 1, 1],
        2: [2, 2, 2],
        3: [3, 3, 3],
        4: [4, 4, 4],
    }

    assert result.npartitions == 1


def test_groupby_with_scheduler_func():
    from dask.threaded import get

    with dask.config.set(scheduler=get):
        b.groupby(lambda x: x, npartitions=1).compute()


def test_concat():
    a = db.from_sequence([1, 2, 3])
    b = db.from_sequence([4, 5, 6])
    c = db.concat([a, b])
    assert list(c) == [1, 2, 3, 4, 5, 6]
    assert c.name == db.concat([a, b]).name


def test_flatten():
    b = db.from_sequence([[1], [2, 3]])
    assert list(b.flatten()) == [1, 2, 3]
    assert b.flatten().name == b.flatten().name


def test_concat_after_map():
    a = db.from_sequence([1, 2])
    b = db.from_sequence([4, 5])
    result = db.concat([a.map(inc), b])
    assert list(result) == [2, 3, 4, 5]


def test_args():
    c = b.map(lambda x: x + 1)
    d = Bag(*c._args)

    assert list(c) == list(d)
    assert c.npartitions == d.npartitions


def test_to_dataframe():
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")

    def check_parts(df, sol):
        assert all(
            (p.dtypes == sol.dtypes).all() for p in dask.compute(*df.to_delayed())
        )

    dsk = {("test", 0): [(1, 2)], ("test", 1): [], ("test", 2): [(10, 20), (100, 200)]}
    b = Bag(dsk, "test", 3)
    sol = pd.DataFrame(b.compute(), columns=["a", "b"])

    # Elements are tuples
    df = b.to_dataframe()
    dd.utils.assert_eq(df, sol.rename(columns={"a": 0, "b": 1}), check_index=False)
    df = b.to_dataframe(columns=["a", "b"])
    dd.utils.assert_eq(df, sol, check_index=False)
    check_parts(df, sol)
    df = b.to_dataframe(meta=[("a", "i8"), ("b", "i8")])
    dd.utils.assert_eq(df, sol, check_index=False)
    check_parts(df, sol)

    # Elements are dictionaries
    b = b.map(lambda x: dict(zip(["a", "b"], x)))
    df = b.to_dataframe()
    dd.utils.assert_eq(df, sol, check_index=False)
    check_parts(df, sol)
    assert df._name == b.to_dataframe()._name

    # With metadata specified
    for meta in [sol, [("a", "i8"), ("b", "i8")]]:
        df = b.to_dataframe(meta=meta)
        dd.utils.assert_eq(df, sol, check_index=False)
        check_parts(df, sol)

    # Error to specify both columns and meta
    with pytest.raises(ValueError):
        b.to_dataframe(columns=["a", "b"], meta=sol)

    # Inference fails if empty first partition
    b2 = b.filter(lambda x: x["a"] > 200)
    with pytest.raises(ValueError):
        b2.to_dataframe()

    # Single column
    b = b.pluck("a")
    sol = sol[["a"]]
    df = b.to_dataframe(meta=sol)
    dd.utils.assert_eq(df, sol, check_index=False)
    check_parts(df, sol)

    # Works with iterators and tuples
    sol = pd.DataFrame({"a": range(100)})
    b = db.from_sequence(range(100), npartitions=5)
    for f in [iter, tuple]:
        df = b.map_partitions(f).to_dataframe(meta=sol)
        dd.utils.assert_eq(df, sol, check_index=False)
        check_parts(df, sol)


ext_open = [("gz", GzipFile), ("bz2", BZ2File), ("", open)]


@pytest.mark.parametrize("ext,myopen", ext_open)
def test_to_textfiles(ext, myopen):
    b = db.from_sequence(["abc", "123", "xyz"], npartitions=2)
    with tmpdir() as dir:
        c = b.to_textfiles(os.path.join(dir, "*." + ext), compute=False)
        dask.compute(*c, scheduler="sync")
        assert os.path.exists(os.path.join(dir, "1." + ext))

        f = myopen(os.path.join(dir, "1." + ext), "rb")
        text = f.read()
        if hasattr(text, "decode"):
            text = text.decode()
        assert "xyz" in text
        f.close()


def test_to_textfiles_name_function_preserves_order():
    seq = [
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
    ]
    b = db.from_sequence(seq, npartitions=16)
    with tmpdir() as dn:
        b.to_textfiles(dn)

        out = (
            db.read_text(os.path.join(dn, "*"), encoding="ascii")
            .map(str)
            .map(str.strip)
            .compute()
        )
        assert seq == out


def test_to_textfiles_name_function_warn():
    seq = [
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
    ]
    a = db.from_sequence(seq, npartitions=16)
    with tmpdir() as dn:
        a.to_textfiles(dn, name_function=str)


def test_to_textfiles_encoding():
    b = db.from_sequence(["汽车", "苹果", "天气"], npartitions=2)
    for ext, myopen in ext_open:
        with tmpdir() as dir:
            c = b.to_textfiles(
                os.path.join(dir, "*." + ext), encoding="gb18030", compute=False
            )
            dask.compute(*c)
            assert os.path.exists(os.path.join(dir, "1." + ext))

            f = myopen(os.path.join(dir, "1." + ext), "rb")
            text = f.read()
            if hasattr(text, "decode"):
                text = text.decode("gb18030")
            assert "天气" in text
            f.close()


def test_to_textfiles_inputs():
    B = db.from_sequence(["abc", "123", "xyz"], npartitions=2)
    with tmpfile() as a:
        with tmpfile() as b:
            B.to_textfiles([a, b])
            assert os.path.exists(a)
            assert os.path.exists(b)

    with tmpdir() as dirname:
        B.to_textfiles(dirname)
        assert os.path.exists(dirname)
        assert os.path.exists(os.path.join(dirname, "0.part"))

    with pytest.raises(TypeError):
        B.to_textfiles(5)


def test_to_textfiles_endlines():
    b = db.from_sequence(["a", "b", "c"], npartitions=1)
    with tmpfile() as fn:
        for last_endline in False, True:
            b.to_textfiles([fn], last_endline=last_endline)
            with open(fn) as f:
                result = f.readlines()
            assert result == ["a\n", "b\n", "c\n" if last_endline else "c"]


def test_string_namespace():
    b = db.from_sequence(["Alice Smith", "Bob Jones", "Charlie Smith"], npartitions=2)

    assert "split" in dir(b.str)
    assert "match" in dir(b.str)

    assert list(b.str.lower()) == ["alice smith", "bob jones", "charlie smith"]
    assert list(b.str.split(" ")) == [
        ["Alice", "Smith"],
        ["Bob", "Jones"],
        ["Charlie", "Smith"],
    ]
    assert list(b.str.match("*Smith")) == ["Alice Smith", "Charlie Smith"]

    pytest.raises(AttributeError, lambda: b.str.sfohsofhf)
    assert b.str.match("*Smith").name == b.str.match("*Smith").name
    assert b.str.match("*Smith").name != b.str.match("*John").name


def test_string_namespace_with_unicode():
    b = db.from_sequence(["Alice Smith", "Bob Jones", "Charlie Smith"], npartitions=2)
    assert list(b.str.lower()) == ["alice smith", "bob jones", "charlie smith"]


def test_str_empty_split():
    b = db.from_sequence(["Alice Smith", "Bob Jones", "Charlie Smith"], npartitions=2)
    assert list(b.str.split()) == [
        ["Alice", "Smith"],
        ["Bob", "Jones"],
        ["Charlie", "Smith"],
    ]


def test_map_with_iterator_function():
    b = db.from_sequence([[1, 2, 3], [4, 5, 6]], npartitions=2)

    def f(L):
        for x in L:
            yield x + 1

    c = b.map(f)

    assert list(c) == [[2, 3, 4], [5, 6, 7]]


def test_ensure_compute_output_is_concrete():
    b = db.from_sequence([1, 2, 3])
    result = b.map(lambda x: x + 1).compute()
    assert not isinstance(result, Iterator)


class BagOfDicts(db.Bag):
    def get(self, key, default=None):
        return self.map(lambda d: d.get(key, default))

    def set(self, key, value):
        def setter(d):
            d[key] = value
            return d

        return self.map(setter)


def test_bag_class_extend():
    dictbag = BagOfDicts(*db.from_sequence([{"a": {"b": "c"}}])._args)
    assert dictbag.get("a").get("b").compute()[0] == "c"
    assert dictbag.get("a").set("d", "EXTENSIBILITY!!!").compute()[0] == {
        "b": "c",
        "d": "EXTENSIBILITY!!!",
    }
    assert isinstance(dictbag.get("a").get("b"), BagOfDicts)


def test_gh715():
    bin_data = "\u20ac".encode()
    with tmpfile() as fn:
        with open(fn, "wb") as f:
            f.write(bin_data)
        a = db.read_text(fn)
        assert a.compute()[0] == bin_data.decode("utf-8")


def test_bag_compute_forward_kwargs():
    x = db.from_sequence([1, 2, 3]).map(lambda a: a + 1)
    x.compute(bogus_keyword=10)


def test_to_delayed():
    b = db.from_sequence([1, 2, 3, 4, 5, 6], npartitions=3)
    a, b, c = b.map(inc).to_delayed()
    assert all(isinstance(x, Delayed) for x in [a, b, c])
    assert b.compute() == [4, 5]

    b = db.from_sequence([1, 2, 3, 4, 5, 6], npartitions=3)
    t = b.sum().to_delayed()
    assert isinstance(t, Delayed)
    assert t.compute() == 21


def test_to_delayed_optimize_graph(tmpdir, monkeypatch):
    calls = 0
    from dask.bag.core import reify

    def count_reify(*args, **kwargs):
        nonlocal calls
        calls += 1
        return reify(*args, **kwargs)

    monkeypatch.setattr(dask.bag.core, "reify", count_reify)
    b = db.from_sequence([1, 2, 3, 4, 5, 6], npartitions=1)
    b2 = b.map(inc).map(inc).map(inc)

    [d] = b2.to_delayed()

    calls = 0
    d.compute()
    assert calls == 1

    calls = 0
    b = db.from_sequence([1, 2, 3, 4, 5, 6], npartitions=1)
    b2 = b.map(inc).map(inc).map(inc)
    [d] = b2.to_delayed()
    assert d.__dask_layers__() != b2.__dask_layers__()
    [d2] = b2.to_delayed(optimize_graph=False)
    assert dict(d2.dask) == dict(b2.dask)
    assert d2.__dask_layers__() == b2.__dask_layers__()
    calls = 0
    result = d.compute()
    assert calls == 1
    calls = 0
    result2 = d2.compute()
    assert calls == 3
    assert result == result2

    x = b2.sum()
    d = x.to_delayed()
    assert d.__dask_layers__() == x.__dask_layers__()
    d2 = x.to_delayed(optimize_graph=False)
    assert dict(d2.dask) == dict(x.dask)
    assert d2.__dask_layers__() == x.__dask_layers__()

    calls = 0
    assert d.compute() == d2.compute()
    assert calls == 3

    calls = 0
    [d] = b2.to_textfiles(str(tmpdir), compute=False)
    assert calls == 0


def test_from_delayed():
    from dask.delayed import delayed

    a, b, c = delayed([1, 2, 3]), delayed([4, 5, 6]), delayed([7, 8, 9])
    bb = from_delayed([a, b, c])
    assert bb.name == from_delayed([a, b, c]).name

    assert isinstance(bb, Bag)
    assert list(bb) == [1, 2, 3, 4, 5, 6, 7, 8, 9]

    asum_value = delayed(sum)(a)
    asum_item = db.Item.from_delayed(asum_value)
    assert asum_value.compute() == asum_item.compute() == 6


def test_from_delayed_iterator():
    from dask.delayed import delayed

    def lazy_records(n):
        return ({"operations": [1, 2]} for _ in range(n))

    delayed_records = delayed(lazy_records, pure=False)
    bag = db.from_delayed([delayed_records(5) for _ in range(5)])
    assert db.compute(
        bag.count(),
        bag.pluck("operations").count(),
        bag.pluck("operations").flatten().count(),
        scheduler="sync",
    ) == (25, 25, 50)


def test_range():
    for npartitions in [1, 7, 10, 28]:
        b = db.range(100, npartitions=npartitions)
        assert len(b.dask) == npartitions
        assert b.npartitions == npartitions
        assert list(b) == list(range(100))


@pytest.mark.parametrize("npartitions", [1, 7, 10, 28])
def test_zip(npartitions, hi=1000):
    evens = db.from_sequence(range(0, hi, 2), npartitions=npartitions)
    odds = db.from_sequence(range(1, hi, 2), npartitions=npartitions)
    pairs = db.zip(evens, odds)
    assert pairs.npartitions == evens.npartitions
    assert pairs.npartitions == odds.npartitions
    assert list(pairs) == list(zip(range(0, hi, 2), range(1, hi, 2)))


@pytest.mark.parametrize("nin", [1, 2, 7, 11, 23])
@pytest.mark.parametrize("nout", [1, 2, 5, 12, 23])
def test_repartition_npartitions(nin, nout):
    b = db.from_sequence(range(100), npartitions=nin)
    c = b.repartition(npartitions=nout)
    assert c.npartitions == nout
    assert_eq(b, c)
    results = dask.get(c.dask, c.__dask_keys__())
    assert all(results)


@pytest.mark.parametrize(
    "nin, nout",
    [
        (1, 1),
        (2, 1),
        (5, 1),
        (1, 2),
        (2, 2),
        (5, 2),
        (1, 5),
        (2, 5),
        (5, 5),
    ],
)
def test_repartition_partition_size(nin, nout):
    b = db.from_sequence(range(1, 100), npartitions=nin)
    total_mem = sum(b.map_partitions(total_mem_usage).compute())
    c = b.repartition(partition_size=(total_mem // nout))
    assert c.npartitions >= nout
    assert_eq(b, c)


def test_multiple_repartition_partition_size():
    b = db.from_sequence(range(1, 100), npartitions=1)
    total_mem = sum(b.map_partitions(total_mem_usage).compute())

    c = b.repartition(partition_size=(total_mem // 2))
    assert c.npartitions >= 2
    assert_eq(b, c)

    d = c.repartition(partition_size=(total_mem // 5))
    assert d.npartitions >= 5
    assert_eq(c, d)


def test_repartition_partition_size_complex_dtypes():
    np = pytest.importorskip("numpy")

    b = db.from_sequence([np.array(range(100)) for _ in range(4)], npartitions=1)
    total_mem = sum(b.map_partitions(total_mem_usage).compute())

    new_partition_size = total_mem // 4
    c = b.repartition(partition_size=new_partition_size)
    assert c.npartitions >= 4
    assert_eq(b, c)


def test_repartition_names():
    b = db.from_sequence(range(100), npartitions=5)
    c = b.repartition(2)
    assert b.name != c.name

    d = b.repartition(20)
    assert b.name != c.name
    assert c.name != d.name

    c = b.repartition(5)
    assert b is c


def test_repartition_input_errors():
    with pytest.raises(ValueError):
        bag = db.from_sequence(range(10))
        bag.repartition(npartitions=5, partition_size="5MiB")


def test_accumulate():
    parts = [[1, 2, 3], [4, 5], [], [6, 7]]
    dsk = {("test", i): p for (i, p) in enumerate(parts)}
    b = db.Bag(dsk, "test", len(parts))
    r = b.accumulate(add)
    assert r.name == b.accumulate(add).name
    assert r.name != b.accumulate(add, -1).name
    assert r.compute() == [1, 3, 6, 10, 15, 21, 28]
    assert b.accumulate(add, -1).compute() == [-1, 0, 2, 5, 9, 14, 20, 27]
    assert b.accumulate(add).map(inc).compute() == [2, 4, 7, 11, 16, 22, 29]

    b = db.from_sequence([1, 2, 3], npartitions=1)
    assert b.accumulate(add).compute() == [1, 3, 6]


def test_groupby_tasks():
    b = db.from_sequence(range(160), npartitions=4)
    out = b.groupby(lambda x: x % 10, max_branch=4, shuffle="tasks")
    partitions = dask.get(out.dask, out.__dask_keys__())

    for a in partitions:
        for b in partitions:
            if a is not b:
                assert not set(pluck(0, a)) & set(pluck(0, b))

    b = db.from_sequence(range(1000), npartitions=100)
    out = b.groupby(lambda x: x % 123, shuffle="tasks")
    assert len(out.dask) < 100**2
    partitions = dask.get(out.dask, out.__dask_keys__())

    for a in partitions:
        for b in partitions:
            if a is not b:
                assert not set(pluck(0, a)) & set(pluck(0, b))

    b = db.from_sequence(range(10000), npartitions=345)
    out = b.groupby(lambda x: x % 2834, max_branch=24, shuffle="tasks")
    partitions = dask.get(out.dask, out.__dask_keys__())

    for a in partitions:
        for b in partitions:
            if a is not b:
                assert not set(pluck(0, a)) & set(pluck(0, b))


def test_groupby_tasks_names():
    b = db.from_sequence(range(160), npartitions=4)
    func = lambda x: x % 10
    func2 = lambda x: x % 20
    assert set(b.groupby(func, max_branch=4, shuffle="tasks").dask) == set(
        b.groupby(func, max_branch=4, shuffle="tasks").dask
    )
    assert set(b.groupby(func, max_branch=4, shuffle="tasks").dask) != set(
        b.groupby(func, max_branch=2, shuffle="tasks").dask
    )
    assert set(b.groupby(func, max_branch=4, shuffle="tasks").dask) != set(
        b.groupby(func2, max_branch=4, shuffle="tasks").dask
    )


@pytest.mark.parametrize(
    "size,npartitions,groups", [(1000, 20, 100), (12345, 234, 1042), (100, 1, 50)]
)
def test_groupby_tasks_2(size, npartitions, groups):
    func = lambda x: x % groups
    b = db.range(size, npartitions=npartitions).groupby(func, shuffle="tasks")
    result = b.compute(scheduler="sync")
    assert dict(result) == groupby(func, range(size))


def test_groupby_tasks_3():
    func = lambda x: x % 10
    b = db.range(20, npartitions=5).groupby(func, shuffle="tasks", max_branch=2)
    result = b.compute(scheduler="sync")
    assert dict(result) == groupby(func, range(20))
    # assert b.npartitions == 5


def test_to_textfiles_empty_partitions():
    with tmpdir() as d:
        b = db.range(5, npartitions=5).filter(lambda x: x == 1).map(str)
        b.to_textfiles(os.path.join(d, "*.txt"))
        assert len(os.listdir(d)) == 5


def test_reduction_empty():
    b = db.from_sequence(range(10), npartitions=100)
    assert_eq(b.filter(lambda x: x % 2 == 0).max(), 8)
    assert_eq(b.filter(lambda x: x % 2 == 0).min(), 0)


@pytest.mark.parametrize("npartitions", [1, 2, 4])
def test_reduction_empty_aggregate(npartitions):
    b = db.from_sequence([0, 0, 0, 1], npartitions=npartitions).filter(None)
    assert_eq(b.min(split_every=2), 1)
    vals = db.compute(b.min(split_every=2), b.max(split_every=2), scheduler="sync")
    assert vals == (1, 1)
    with pytest.raises(ValueError):
        b = db.from_sequence([0, 0, 0, 0], npartitions=npartitions)
        b.filter(None).min(split_every=2).compute(scheduler="sync")


class StrictReal(int):
    def __eq__(self, other):
        assert isinstance(other, StrictReal)
        return self.real == other.real

    def __ne__(self, other):
        assert isinstance(other, StrictReal)
        return self.real != other.real


def test_reduction_with_non_comparable_objects():
    b = db.from_sequence([StrictReal(x) for x in range(10)], partition_size=2)
    assert_eq(b.fold(max, max), StrictReal(9))


@pytest.mark.parametrize("container", ["array", "matrix"])
def test_reduction_with_sparse_matrices(container):
    sp = pytest.importorskip("scipy.sparse")
    if container == "array" and not hasattr(sp, "csr_array"):
        pytest.skip("scipy<1.11 has no sparray")
    cls = sp.csr_matrix if container == "matrix" else sp.csr_array

    b = db.from_sequence([cls([0]) for x in range(4)], partition_size=2)

    def sp_reduce(a, b):
        return sp.vstack([a, b])

    assert b.fold(sp_reduce, sp_reduce).compute(scheduler="sync").shape == (4, 1)


def test_empty():
    assert list(db.from_sequence([])) == []


def test_bag_picklable():
    from pickle import dumps, loads

    b = db.from_sequence(range(100))
    b2 = loads(dumps(b))
    assert b.compute() == b2.compute()

    s = b.sum()
    s2 = loads(dumps(s))
    assert s.compute() == s2.compute()


def test_msgpack_unicode():
    b = db.from_sequence([{"a": 1}]).groupby("a")
    result = b.compute(scheduler="sync")
    assert dict(result) == {1: [{"a": 1}]}


def test_bag_with_single_callable():
    f = lambda: None
    b = db.from_sequence([f])
    assert_eq(b, [f])


def test_optimize_fuse_keys():
    x = db.range(10, npartitions=2)
    y = x.map(inc)
    z = y.map(inc)

    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__())
    assert not y.dask.keys() & dsk.keys()

    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__(), fuse_keys=y.__dask_keys__())
    assert all(k in dsk for k in y.__dask_keys__())


def test_reductions_are_lazy():
    current = [None]

    def part():
        for i in range(10):
            current[0] = i
            yield i

    def func(part):
        assert current[0] == 0
        return sum(part)

    b = Bag({("foo", 0): part()}, "foo", 1)

    res = b.reduction(func, sum)

    assert_eq(res, sum(range(10)))


def test_repeated_groupby():
    b = db.range(10, npartitions=4)
    c = b.groupby(lambda x: x % 3)
    assert valmap(len, dict(c)) == valmap(len, dict(c))


def test_temporary_directory(tmpdir):
    b = db.range(10, npartitions=4)

    # We use a pool to avoid a race condition between the pool close
    # cleaning up files, and the assert below.
    with ProcessPoolExecutor(4) as pool:
        with dask.config.set(temporary_directory=str(tmpdir), pool=pool):
            b2 = b.groupby(lambda x: x % 2)
            b2.compute()
            assert any(fn.endswith(".partd") for fn in os.listdir(str(tmpdir)))


def test_empty_bag():
    b = db.from_sequence([])
    assert_eq(b.map(inc).all(), True)
    assert_eq(b.map(inc).any(), False)
    assert_eq(b.map(inc).sum(), False)
    assert_eq(b.map(inc).count(), False)


def test_bag_paths():
    b = db.from_sequence(["abc", "123", "xyz"], npartitions=2)
    paths = b.to_textfiles("foo*")
    assert paths[0].endswith("foo0")
    assert paths[1].endswith("foo1")

    os.remove("foo0")
    os.remove("foo1")


def test_map_partitions_arg():
    def append_str(partition, s):
        return [x + s for x in partition]

    mybag = db.from_sequence(["a", "b", "c"])

    assert_eq(mybag.map_partitions(append_str, "foo"), ["afoo", "bfoo", "cfoo"])
    assert_eq(
        mybag.map_partitions(append_str, dask.delayed("foo")), ["afoo", "bfoo", "cfoo"]
    )


def test_map_keynames():
    b = db.from_sequence([1, 2, 3])
    d = dict(b.map(inc).__dask_graph__())
    assert "inc" in map(dask.utils.key_split, d)

    assert set(b.map(inc).__dask_graph__()) != set(
        b.map_partitions(inc).__dask_graph__()
    )


def test_map_releases_element_references_as_soon_as_possible():
    # Ensure that Bag.map doesn't keep *element* references longer than
    # necessary. Previous map implementations used ``yield``, which would keep
    # a reference to the yielded element until the yielded method resumed (this
    # is just how generator functions work in CPython).
    #
    # See https://github.com/dask/dask/issues/5189
    #
    # We test 2 variant of potential extra references here:
    # 1. Within an element of a partition:
    #    At the time of the second `f_create` for each element, the `C` from
    #    the first `f_create` should be dropped.
    # 2. Within a partition:
    #    When the second item within a partition is processed, `C` from the
    #    first item should already be dropped.
    class C:
        def __init__(self, i):
            self.i = i

    # keep a weakref to all existing instances of `C`
    in_memory = weakref.WeakSet()

    def f_create(i):
        # check that there are no instances of `C` left
        assert len(in_memory) == 0

        # create new instance
        o = C(i)
        in_memory.add(o)

        return o

    def f_drop(o):
        # o reference dropped on return, should collect
        return o.i + 100

    b = (
        db.from_sequence(range(2), npartitions=1)
        .map(f_create)
        .map(f_drop)
        .map(f_create)
        .map(f_drop)
        .sum()
    )
    try:
        # Disable gc to ensure refcycles don't matter here
        gc.disable()
        b.compute(scheduler="sync")
    finally:
        gc.enable()


def test_bagged_array_delayed():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")

    obj = da.ones(10, chunks=5).to_delayed()[0]
    bag = db.from_delayed(obj)
    b = bag.compute()
    assert_eq(b, [1.0, 1.0, 1.0, 1.0, 1.0])


def test_dask_layers():
    a = db.from_sequence([1, 2], npartitions=2)
    assert a.__dask_layers__() == (a.name,)
    assert a.dask.layers.keys() == {a.name}
    assert a.dask.dependencies == {a.name: set()}
    i = a.min()
    assert i.__dask_layers__() == (i.key,)
    assert i.dask.layers.keys() == {a.name, i.key}
    assert i.dask.dependencies == {a.name: set(), i.key: {a.name}}


@pytest.mark.parametrize("optimize", [False, True])
def test_dask_layers_to_delayed(optimize):
    # `da.Array.to_delayed` causes the layer name to not match the key.
    # Ensure the layer name is propagated between `Delayed` and `Item`.
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    i = db.Item.from_delayed(da.ones(1).to_delayed()[0])
    name = i.key[0]
    assert i.key[1:] == (0,)
    assert i.dask.layers.keys() == {"delayed-" + name}
    assert i.dask.dependencies == {"delayed-" + name: set()}
    assert i.__dask_layers__() == ("delayed-" + name,)

    arr = da.ones(1) + 1
    delayed = arr.to_delayed(optimize_graph=optimize)[0]
    i = db.Item.from_delayed(delayed)
    assert i.key == delayed.key
    assert i.dask is delayed.dask
    assert i.__dask_layers__() == delayed.__dask_layers__()

    back = i.to_delayed(optimize_graph=optimize)
    assert back.__dask_layers__() == i.__dask_layers__()

    if not optimize:
        assert back.dask is arr.dask
        # When not optimized, the key is not a layer in the graph, so using it should fail
        with pytest.raises(ValueError, match="not in"):
            db.Item(back.dask, back.key)

    with pytest.raises(ValueError, match="not in"):
        db.Item(arr.dask, (arr.name,), layer="foo")


def test_to_dataframe_optimize_graph():
    pytest.importorskip("pandas")
    pytest.importorskip("dask.dataframe")

    from dask.dataframe.utils import assert_eq as assert_eq_df

    x = db.from_sequence(
        [{"name": "test1", "v1": 1}, {"name": "test2", "v1": 2}], npartitions=2
    )

    # linear `map` tasks will be fused by graph optimization
    with dask.annotate(foo=True):
        y = x.map(lambda a: dict(**a, v2=a["v1"] + 1))
        y = y.map(lambda a: dict(**a, v3=a["v2"] + 1))
        y = y.map(lambda a: dict(**a, v4=a["v3"] + 1))

    # verifying the maps are not fused yet
    assert len(y.dask) == y.npartitions * 4

    # with optimizations
    d = y.to_dataframe()
    # no optimizations
    d2 = y.to_dataframe(optimize_graph=False)
    assert_eq_df(d, d2)


@pytest.mark.parametrize("nworkers", [100, 250, 500, 1000])
def test_default_partitioning_worker_saturation(nworkers):
    # Ensure that Dask Bag can saturate any number of workers with concurrent tasks.
    # The default partitioning scheme partitions items to keep the task to item ratio sensible
    # but it should always be possible to saturate any number of workers given enough items in the bag.
    ntasks = 0
    nitems = 1
    while ntasks < nworkers:
        ntasks = len(db.from_sequence(range(nitems)).dask)
        nitems += math.floor(max(1, nworkers / 10))
        assert nitems < 20_000


@pytest.mark.parametrize("nworkers", [100, 250, 500, 1000])
def test_npartitions_saturation(nworkers):
    # If npartitions is set the bag should always contain at least that number of tasks
    for nitems in range(nworkers, 10 * nworkers, max(1, math.floor(nworkers / 10))):
        assert (
            len(db.from_sequence(range(nitems), npartitions=nworkers).dask) >= nworkers
        )


def test_map_total_mem_usage():
    """https://github.com/dask/dask/issues/10338"""
    b = db.from_sequence(range(1, 100), npartitions=3)
    total_mem_b = sum(b.map_partitions(total_mem_usage).compute())
    c = b.map(lambda x: x)
    total_mem_c = sum(c.map_partitions(total_mem_usage).compute())
    assert total_mem_b == total_mem_c
