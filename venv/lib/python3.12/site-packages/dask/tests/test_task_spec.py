from __future__ import annotations

import itertools
import pickle
import sys
from collections import namedtuple
from collections.abc import Mapping

import pytest

from dask._task_spec import (
    Alias,
    DataNode,
    DependenciesMapping,
    Dict,
    List,
    Set,
    Task,
    TaskRef,
    Tuple,
    _get_dependencies,
    convert_legacy_graph,
    execute_graph,
    fuse_linear_task_spec,
    parse_input,
    resolve_aliases,
)
from dask.core import keys_in_tasks, reverse_dict
from dask.sizeof import sizeof
from dask.tokenize import tokenize
from dask.utils import funcname


def convert_and_verify_keys(dsk):
    new_dsk = convert_legacy_graph(dsk)
    vals = [(k, v.key) for k, v in new_dsk.items()]
    assert all(v[0] == v[1] for v in vals), vals
    return new_dsk


def identity(x):
    return x


def func(*args):
    try:
        return "-".join(args)
    except TypeError:
        return "-".join(map(str, args))


def func2(*args):
    return "=".join(args)


def func3(*args, **kwargs):
    return "+".join(args) + "//" + "+".join(f"{k}={v}" for k, v in kwargs.items())


def test_convert_legacy_dsk_skip_new():
    dsk = {
        "key-1": Task("key-1", func, "a", "b"),
    }
    converted = convert_and_verify_keys(dsk)
    assert converted["key-1"] is dsk["key-1"]
    assert converted == dsk


def test_repr():
    t = Task("key", func, "a", "b")
    assert repr(t) == "<Task 'key' func('a', 'b')>"

    t = Task("nested", func2, t, t.ref())
    assert repr(t) == "<Task 'nested' func2(<Task 'key' func('a', 'b')>, Alias('key'))>"

    def long_function_name_longer_even_longer(a, b):
        return a + b

    t = Task("long", long_function_name_longer_even_longer, t, t.ref())
    assert repr(t) == "<Task 'long' long_function_name_longer_even_longer(...)>"

    def use_kwargs(a, kwarg=None):
        return a + kwarg

    t = Task("kwarg", use_kwargs, "foo", kwarg="kwarg_value")
    assert repr(t) == "<Task 'kwarg' use_kwargs('foo', kwarg='kwarg_value')>"


def test_task_eq():
    assert Alias(("x", 0, 0)) == Alias(("x", 0, 0))
    assert Alias(("x", 0, 0)) != Alias(("x", 0, 1))

    assert List(Alias(("x", 0, 0))) == List(Alias(("x", 0, 0)))
    assert List(Alias(("x", 0, 0))) == List(Alias("a", "a")).substitute(
        {"a": ("x", 0, 0)}
    )
    assert List(Alias(("x", 0, 0)), Alias(("x", 0, 1))) != List(
        Alias(("x", 0, 0)), Alias(("x", 0, 2))
    )

    obj = Task(
        ("z", 0, 0),
        func,
        List(TaskRef(("x", 0, 0)), TaskRef(("x", 0, 1))),
        List(TaskRef(("y", 0, 0)), TaskRef(("y", 1, 0))),
    )
    assert obj == obj
    assert obj == Task(
        ("z", 0, 0),
        func,
        List(TaskRef(("x", 0, 0)), TaskRef(("x", 0, 1))),
        List(TaskRef(("y", 0, 0)), TaskRef(("y", 1, 0))),
    )


def long_function_name_longer_even_longer(a, b):
    return a + b


def use_kwargs(a, kwarg=None):
    return a + kwarg


def test_convert_legacy_dsk():
    def func(*args):
        return "-".join(args)

    def func2(*args):
        return "=".join(args)

    dsk = {
        "key-1": (func, "a", "b"),
        "key-2": (func2, "key-1", "c"),
        "key-3": (func, (func2, "c", "key-1"), "key-2"),
        "const": "foo",
        "key-4": [
            (func, "key-1", "b"),
            (func, "c", "key-2"),
            (func, "key-3", "f"),
            (func, "const", "bar"),
        ],
    }
    new_dsk = convert_and_verify_keys(dsk)
    t1 = new_dsk["key-1"]
    assert isinstance(t1, Task)
    assert t1.func == func
    v1 = t1()
    assert v1 == func("a", "b")

    t2 = new_dsk["key-2"]
    assert isinstance(t2, Task)
    v2 = t2({"key-1": v1})
    assert v2 == func2(func("a", "b"), "c")

    t3 = new_dsk["key-3"]
    assert isinstance(t3, Task)
    assert len(t3.dependencies) == 2
    v3 = t3({"key-1": v1, "key-2": v2})
    assert v3 == func(func2("c", func("a", "b")), func2(func("a", "b"), "c"))
    t4 = new_dsk["key-4"]
    assert isinstance(t4, Task)
    assert t4.dependencies == {"key-1", "key-2", "key-3", "const"}
    assert t4(
        {
            "key-1": v1,
            "key-2": v2,
            "key-3": v3,
            "const": "foo",
        }
    ) == [
        func(v1, "b"),
        func("c", v2),
        func(v3, "f"),
        func("foo", "bar"),
    ]


def test_task_executable():
    t1 = Task("key-1", func, "a", "b")
    assert t1() == func("a", "b")


def test_task_nested_sequence():
    # We do **not** recurse into structures since this kills performance
    # we're only considering top-level task objects as arguments or kwargs
    t1 = Task("key-1", func, "a", "b")
    t2 = Task("key-2", func2, "c", "d")

    tseq = Task("seq", identity, [t1, t2])
    assert tseq() == [t1, t2]

    tseq = Task("set", identity, {t1, t2})
    assert tseq() == {t1, t2}

    tseq = Task("dict", identity, {"a": t1, "b": t2})
    assert tseq() == {"a": t1, "b": t2}


def test_reference_remote():
    t1 = Task("key-1", func, "a", "b")
    t2 = Task("key-2", func2, TaskRef("key-1"), "c")
    with pytest.raises(RuntimeError):
        t2()
    assert t2({"key-1": t1()}) == func2(func("a", "b"), "c")


def test_reference_remote_twice_same():
    t1 = Task("key-1", func, "a", "b")
    t2 = Task("key-2", func2, TaskRef("key-1"), TaskRef("key-1"))
    with pytest.raises(RuntimeError):
        t2()
    assert t2({"key-1": t1()}) == func2(func("a", "b"), func("a", "b"))


def test_reference_remote_twice_different():
    t1 = Task("key-1", func, "a", "b")
    t2 = Task("key-2", func, "c", "d")
    t3 = Task("key-3", func2, TaskRef("key-1"), TaskRef("key-2"))
    with pytest.raises(RuntimeError):
        t3()
    assert t3({"key-1": t1(), "key-2": t2()}) == func2(func("a", "b"), func("c", "d"))


def test_task_nested():
    t1 = Task("key-1", func, "a", "b")
    t2 = Task("key-2", func2, t1, "c")
    assert t2() == func2(func("a", "b"), "c")


class SerializeOnlyOnce:
    deserialized = False
    serialized = False

    def __getstate__(self):
        if SerializeOnlyOnce.serialized:
            raise RuntimeError()
        SerializeOnlyOnce.serialized = True
        return {}

    def __setstate__(self, state):
        if SerializeOnlyOnce.deserialized:
            raise RuntimeError()
        SerializeOnlyOnce.deserialized = True

    def __call__(self, a, b):
        return a + b


def test_pickle():

    def assert_slots_equal(a, b):
        def get_all_slots(obj):
            slots = set()
            for cls in obj.__class__.mro():
                slots.update(getattr(cls, "__slots__", ()))
            return slots

        all_slots = get_all_slots(a) | get_all_slots(a)
        assert all_slots == get_all_slots(a) == get_all_slots(a)
        assert all(getattr(a, slot) == getattr(b, slot) for slot in all_slots)
        assert not hasattr(a, "__dict__")
        assert not hasattr(b, "__dict__")

    t1 = Task("key-1", func, "a", "b")
    t2 = Task("key-2", func, "c", "d")

    rtt1 = pickle.loads(pickle.dumps(t1))
    assert repr(rtt1) == repr(t1)
    rtt2 = pickle.loads(pickle.dumps(t2))

    assert_slots_equal(t1, rtt1)
    assert_slots_equal(t2, rtt2)
    assert t1 == rtt1
    assert t1.func == rtt1.func
    assert t1.func is rtt1.func
    assert t1.func is rtt2.func

    l = Tuple(t1, t2)
    rtl = pickle.loads(pickle.dumps(l))
    assert l == rtl
    assert l() == rtl()
    d = Dict(key=t1)
    rtd = pickle.loads(pickle.dumps(d))
    assert d == rtd
    assert d() == rtd()


def test_pickle_size():
    # We will serialize many of these objects which drives both memory usage and
    # serialization runtime performance.
    # Reducing pickle size is beneficial but the numbers below are determined
    # empirically
    # Analyzing the output with pickletools.dis is useful to debug memoization
    # and serialization by value

    a = Alias("a", "b")
    # We cannot shrink it to nothing
    assert len(pickle.dumps(a)) < 55
    b = Alias("b", "c")
    # But most of it should be overhead that is memoized
    assert len(pickle.dumps((a, b))) <= 70

    # Pickle should be able to memoize this. On py3.10 that's 2 additional bytes
    assert len(pickle.dumps((a, b, b))) <= len(pickle.dumps((a, b))) + 10

    t1 = Task("key-1", func, "a", "b")
    assert len(pickle.dumps(t1)) < 120

    t2 = Task("key-2", func, TaskRef("key-1"), "c")
    assert len(pickle.dumps(t2)) < 140

    assert len(pickle.dumps((t1, t2))) < 170

    l = List(t1, t2)
    assert len(pickle.dumps(l)) <= 272

    sizes = []
    growth = []
    inner = List(t1, t2)
    for depth in range(20):
        inner = List(inner, t1)
        size = len(pickle.dumps(inner))
        if len(sizes) > 0:
            growth.append(size - sizes[-1][1])
        sizes.append((depth, size))
    growth = set(growth)
    # If this breaks, something cannot be memoized. That's very concerning
    assert len(growth) == 1
    # If this goes up, that's not great but not a disaster
    assert growth.pop() <= 32


def test_tokenize():
    t = Task("key-1", func, "a", "b")
    assert tokenize(t) == tokenize(t)

    t2 = Task("key-1", func, "a", "b")
    assert tokenize(t) == tokenize(t2)

    tokenize(t)

    # Literals are often generated with random/anom names but that should not
    # impact hashing. Otherwise identical submits would end up with different
    # tokens
    l = DataNode("key-1", "a")
    l2 = DataNode("key-2", "a")
    assert tokenize(l) == tokenize(l2)


async def afunc(a, b):
    return a + b


def test_async_func():
    pytest.importorskip("distributed")

    from distributed.utils_test import gen_test

    @gen_test()
    async def _():

        t = Task("key-1", afunc, "a", "b")
        assert t.is_coro
        assert await t() == "ab"
        assert await pickle.loads(pickle.dumps(t))() == "ab"

    _()


def test_parse_curry():
    def curry(func, *args, **kwargs):
        return func(*args, **kwargs) + "c"

    dsk = {
        "key-1": (curry, func, "a", "b"),
    }
    converted = convert_and_verify_keys(dsk)
    t = Task("key-1", curry, func, "a", "b")
    assert converted["key-1"]() == t()
    assert t() == "a-bc"


def test_curry():
    def curry(func, *args, **kwargs):
        return func(*args, **kwargs) + "c"

    t = Task("key-1", curry, func, "a", "b")
    assert t() == "a-bc"


def test_avoid_cycles():
    pytest.importorskip("distributed")
    from dask._task_spec import TaskRef

    dsk = {
        "key": TaskRef("key"),  # e.g. a persisted key
    }
    new_dsk = convert_and_verify_keys(dsk)
    assert not new_dsk


def test_runnable_as_kwarg():
    def func_kwarg(a, b, c=""):
        return a + b + str(c)

    t = Task(
        "key-1",
        func_kwarg,
        "a",
        "b",
        c=Task("key-2", sum, [1, 2]),
    )
    assert t() == "ab3"


def test_dependency_as_kwarg():
    def func_kwarg(a, b, c=""):
        return a + b + str(c)

    t1 = Task("key-1", sum, [1, 2])
    t2 = Task(
        "key-2",
        func_kwarg,
        "a",
        "b",
        c=t1.ref(),
    )
    with pytest.raises(RuntimeError, match="missing"):
        t2()
    # It isn't sufficient to raise. We also rely on the attribute being set
    # properly since distribute will use this to infer actual dependencies The
    # exception may be raised recursively
    assert t2.dependencies
    assert t2({"key-1": t1()}) == "ab3"


def test_array_as_argument():
    np = pytest.importorskip("numpy")
    t = Task("key-1", func, np.array([1, 2]), "b")
    assert t() == "[1 2]-b"

    # This will **not** work since we do not want to recurse into an array!
    t2 = Task("key-2", func, np.array([1, t.ref()]), "b")
    assert t2({"key-1": "foo"}) != "[1 foo]-b"
    assert not _get_dependencies(np.array([1, t.ref()]))


@pytest.mark.parametrize(
    "inst",
    [
        Task("key-1", func, "a", "b"),
        Alias("key-1"),
        DataNode("key-1", 1),
    ],
)
def test_ensure_slots(inst):
    assert not hasattr(inst, "__dict__")


class PlainNamedTuple(namedtuple("PlainNamedTuple", "value")):
    """Namedtuple with a default constructor."""


class NewArgsNamedTuple(namedtuple("NewArgsNamedTuple", "ab, c")):
    """Namedtuple with a custom constructor."""

    def __new__(cls, a, b, c):
        return super().__new__(cls, f"{a}-{b}", c)

    def __getnewargs__(self):
        return *self.ab.split("-"), self.c


class NewArgsExNamedTuple(namedtuple("NewArgsExNamedTuple", "ab, c, k, v")):
    """Namedtuple with a custom constructor including keywords-only arguments."""

    def __new__(cls, a, b, c, **kw):
        return super().__new__(cls, f"{a}-{b}", c, tuple(kw.keys()), tuple(kw.values()))

    def __getnewargs_ex__(self):
        return (*self.ab.split("-"), self.c), dict(zip(self.k, self.v))


@pytest.mark.parametrize(
    "typ, args, kwargs",
    [
        (PlainNamedTuple, ["some-data"], {}),
        (NewArgsNamedTuple, ["some", "data", "more"], {}),
        (NewArgsExNamedTuple, ["some", "data", "more"], {"another": "data"}),
    ],
)
def test_parse_graph_namedtuple_legacy(typ, args, kwargs):
    def func(x):
        return x

    dsk = {"foo": (func, typ(*args, **kwargs))}
    new_dsk = convert_and_verify_keys(dsk)

    assert new_dsk["foo"]() == typ(*args, **kwargs)


@pytest.mark.parametrize(
    "typ, args, kwargs",
    [
        (PlainNamedTuple, ["some-data"], {}),
        (NewArgsNamedTuple, ["some", "data", "more"], {}),
        (NewArgsExNamedTuple, ["some", "data", "more"], {"another": "data"}),
    ],
)
def test_parse_namedtuple(typ, args, kwargs):
    def func(x):
        return x

    obj = typ(*args, **kwargs)
    t = Task("foo", func, parse_input(obj))

    assert t() == obj

    # The other test tuple do weird things to their input
    if typ is PlainNamedTuple:
        args = tuple([TaskRef("b")] + list(args)[1:])
        obj = typ(*args, **kwargs)
        t = Task("foo", func, parse_input(obj))
        assert t.dependencies == {"b"}
        assert t({"b": "foo"}) == typ(*tuple(["foo"] + list(args)[1:]), **kwargs)


def test_pickle_literals():
    np = pytest.importorskip("numpy")
    obj = DataNode("foo", np.transpose)
    roundtripped = pickle.loads(pickle.dumps(obj))
    assert roundtripped == obj


@pytest.mark.parametrize("obj", [set, {0}, [], [1], {}, {2: 3}, (), (4,)])
def test_parse_non_task_inputs(obj):
    assert parse_input(obj) == obj


def test_resolve_aliases():
    tasks = [
        Alias("bar", "foo"),
        Task("foo", func, "a", "b"),
        Alias("baz", "bar"),
    ]
    dsk = {t.key: t for t in tasks}
    assert len(dsk) == 3

    optimized = resolve_aliases(dsk, {"baz"}, reverse_dict(DependenciesMapping(dsk)))
    assert len(optimized) == 1
    expected = dsk["foo"].copy()
    expected.key = "baz"
    assert optimized["baz"] == expected

    optimized = resolve_aliases(
        dsk, {"baz", "bar"}, reverse_dict(DependenciesMapping(dsk))
    )
    assert len(optimized) == 2
    expected = dsk["foo"].copy()
    expected.key = "bar"
    assert optimized["bar"] == expected

    tasks = [
        bar := Alias("bar", "foo"),
        Task("foo", func, "a", "b"),
        Alias("baz", bar.ref()),
        Task("foo2", func, bar.ref(), "c"),
    ]
    dsk = {t.key: t for t in tasks}
    optimized = resolve_aliases(
        dsk, {"baz", "foo2"}, reverse_dict(DependenciesMapping(dsk))
    )
    assert len(optimized) == 3
    # FIXME: Ideally, the above example would optimize to this but this isn't
    # implemented. Instead, we'll block to not mess up anything
    # assert sorted(optimized.values(), key=lambda t: t.key) == sorted(
    #     [
    #         Task("baz", func, "a", "b"),
    #         Task("foo", func, TaskRef("baz"), "c"),
    #     ],
    #     key=lambda t: t.key,
    # )
    # `bar` won't be inlined because it's used in `foo2` AND `baz`
    assert "bar" in optimized
    assert optimized["bar"].key == "bar"
    assert "foo" not in optimized

    # Handle cases with external dependencies
    foo = Task("foo", func, "a", TaskRef("b"))
    dsk = {t.key: t for t in [foo]}
    optimized = resolve_aliases(dsk, {"foo"}, reverse_dict(DependenciesMapping(dsk)))
    assert optimized == dsk


def test_resolve_multiple_aliases():

    tasks = [
        Task("first", func, 10),
        Alias("second", "first"),
        Task("third", func, TaskRef("second")),
        Alias("fourth", "third"),
        Task("fifth", func, TaskRef("fourth")),
    ]
    dsk = {t.key: t for t in tasks}
    assert len(dsk) == 5

    optimized = resolve_aliases(dsk, {"fifth"}, reverse_dict(DependenciesMapping(dsk)))
    assert len(optimized) == 3
    expected = dsk["third"].copy()
    expected.key = "fourth"
    assert optimized["fourth"] == expected

    expected = dsk["first"].copy()
    expected.key = "second"
    assert optimized["second"] == expected


def test_convert_resolve():
    dsk = {
        "first": (func, 10),
        "second": "first",
        "third": (func, "second"),
        "fourth": "third",
        "fifth": (func, "fourth"),
    }
    dsk = convert_and_verify_keys(dsk)
    assert len(dsk) == 5

    optimized = resolve_aliases(dsk, {"fifth"}, reverse_dict(DependenciesMapping(dsk)))
    assert len(optimized) == 3
    expected = dsk["third"].copy()
    expected.key = "fourth"
    assert optimized["fourth"] == expected

    expected = dsk["first"].copy()
    expected.key = "second"
    assert optimized["second"] == expected


def test_parse_nested():
    t = Task(
        "key",
        func3,
        x=TaskRef("y"),
    )

    assert t({"y": "y"}) == "//x=y"


class CountSerialization:
    serialization = 0
    deserialization = 0

    def __getstate__(self):
        CountSerialization.serialization += 1
        return "foo"

    def __setstate__(self, state):
        CountSerialization.deserialization += 1
        pass

    def __call__(self):
        return 1


class RaiseOnSerialization:
    def __getstate__(self):
        raise ValueError("Nope")

    def __call__(self):
        return 1


class RaiseOnDeSerialization:
    def __getstate__(self):
        return "Nope"

    def __setstate__(self, state):
        raise ValueError(state)

    def __call__(self):
        return 1


# This is duplicated from distributed/utils_test.py
def _get_gc_overhead():
    class _CustomObject:
        def __sizeof__(self):
            return 0

    return sys.getsizeof(_CustomObject())


_size_obj = _get_gc_overhead()


class SizeOf:
    """
    An object that returns exactly nbytes when inspected by dask.sizeof.sizeof
    """

    def __init__(self, nbytes: int) -> None:
        if not isinstance(nbytes, int):
            raise TypeError(f"Expected integer for nbytes but got {type(nbytes)}")
        if nbytes < _size_obj:
            raise ValueError(
                f"Expected a value larger than {_size_obj} integer but got {nbytes}."
            )
        self._nbytes = nbytes - _size_obj

    def __sizeof__(self) -> int:
        return self._nbytes


def test_sizeof():
    t = Task("key", func, "a", "b")
    assert sizeof(t) >= sizeof(Task) + sizeof(func) + 2 * sizeof("a")

    t = Task("key", func, SizeOf(100_000))
    assert sizeof(t) > 100_000
    t = DataNode("key", SizeOf(100_000))
    assert sizeof(t) > 100_000


def test_execute_tasks_in_graph():
    dsk = [
        t1 := Task("key-1", func, "a", "b"),
        t2 := Task("key-2", func2, t1.ref(), "c"),
        t3 := Task("key-3", func, "foo", "bar"),
        Task("key-4", func, t3.ref(), t2.ref()),
    ]
    res = execute_graph(dsk, keys=["key-4"])
    assert len(res) == 1
    assert res["key-4"] == "foo-bar-a-b=c"


def test_deterministic_tokenization_respected():
    with pytest.raises(RuntimeError, match="deterministic"):
        tokenize(Task("key", func, object()), ensure_deterministic=True)


def test_keys_in_tasks():
    b = Task("b", func, "1", "2")
    b_legacy = Task("b", func, "1", "2")

    a = Task("a", func, "1", b.ref())
    a_legacy = (func, "1", "b")

    for task in [a, a_legacy]:
        assert not keys_in_tasks({"a"}, [task])

        assert keys_in_tasks({"a", "b"}, [task]) == {"b"}

        assert keys_in_tasks({"a", "b", "c"}, [task]) == {"b"}

    for tasks in [[a, b], [a_legacy, b_legacy]]:
        assert keys_in_tasks({"a", "b"}, tasks) == {"b"}


def test_dependencies_mapping_doesnt_mutate_task():

    t = Task("key", func, "a", "b")
    t2 = Task("key2", func, "a", t.ref())

    assert t2.dependencies == {"key"}

    dsk = {t.key: t, t2.key: t2}
    deps = DependenciesMapping(dsk)
    del deps[t.key]
    # The getitem is doing weird stuff and could mutate the task state
    deps[t2.key]
    assert t2.dependencies == {"key"}


def test_fuse_tasks_key():
    a = Task("key-1", func, "a", "b")
    b = Task("key-2", func2, a.ref(), "d")
    for t1, t2 in itertools.permutations((a, b)):
        fused = Task.fuse(t2, t1)
        assert fused.key == b.key

        fused = Task.fuse(t2, t1, key="new-key")
        assert fused.key == "new-key"


def test_fuse_tasks():
    a = Task("key-1", func, "a", "b")
    b = Task("key-2", func2, a.ref(), "d")
    c = Task("key-3", func3, b.ref(), "e")
    for t1, t2, t3 in itertools.permutations((a, b, c)):
        fused = Task.fuse(t3, t2, t1)
        assert fused() == func3(func2(func("a", "b"), "d"), "e")
        t1 = Task("key-1", func, TaskRef("dependency"), "b")
        t2 = Task("key-2", func2, t1.ref(), "d")
        t3 = Task("key-3", func3, t2.ref(), "e")
        fused = Task.fuse(t3, t2, t1)
        assert fused.dependencies == {"dependency"}
        assert fused({"dependency": "dep"}) == func3(func2(func("dep", "b"), "d"), "e")


def test_fuse_reject_multiple_outputs():
    a = Task("key-1", func, "a", "b")
    b = Task("key-2", func2, "a", "d")
    for t1, t2 in itertools.permutations((a, b)):
        with pytest.raises(ValueError, match="multiple outputs"):
            Task.fuse(t1, t2)


def test_fused_ensure_only_executed_once():
    counter = []

    def counter_func(a, b):
        counter.append(None)
        return func(a, b)

    a = Task("key-1", counter_func, "a", "a")
    b = Task("key-2", func2, a.ref(), "b")
    c = Task("key-3", func2, a.ref(), "c")
    d = Task("key-4", func, b.ref(), c.ref())
    for perm in itertools.permutations([a, b, c, d]):
        fused = Task.fuse(*perm)
        counter.clear()
        assert fused() == func(func2(func("a", "a"), "b"), func2(func("a", "a"), "c"))
        assert len(counter) == 1


def test_fused_dont_hold_in_memory_too_long():
    tasks = []
    prev = None

    # If we execute a fused task we want to release objects as quickly as
    # possible. If every task generates this object, we must at most hold two of
    # them in memory
    class OnlyTwice:
        counter = 0
        total = 0

        def __init__(self):
            OnlyTwice.counter += 1
            OnlyTwice.total += 1
            if OnlyTwice.counter > 2:
                raise ValueError("Didn't release as expected")

        def __del__(self):
            OnlyTwice.counter -= 1

    def generate_object(arg):
        return OnlyTwice()

    prev = None
    for ix in range(10):
        prev = t = Task(
            f"key-{ix}", generate_object, prev.ref() if prev is not None else ix
        )
        tasks.append(t)
    fuse = Task.fuse(*tasks)
    assert fuse()
    assert OnlyTwice.total == 10


def test_linear_fusion():
    tasks = [
        io := Task("foo", func, 1),
        second := Task("second", func, io.ref()),
        Task("third", func, second.ref()),
    ]
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"third"})
    assert len(result) == 2
    assert isinstance(result["third"], Alias)
    assert (
        isinstance(result["foo-second-third"], Task)
        and funcname(result["foo-second-third"].func) == "_execute_subgraph"
    )

    # Data Nodes don't get fused

    tasks = [
        io := DataNode("foo", 1),
        second := Task("second", func, io.ref()),
        third := Task("third", func, second.ref()),
        Task("fourth", func, third.ref()),
    ]
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert len(result) == 2
    assert isinstance(result["fourth"], Alias)
    assert (
        isinstance(result["foo-second-third-fourth"], Task)
        and funcname(result["foo-second-third-fourth"].func) == "_execute_subgraph"
    )
    assert "foo" not in result

    # Branch, so no fusion
    tasks.append(Task("branch", func, third.ref()))
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert len(result) == 4
    assert "foo-second-third" in result
    assert isinstance(result["third"], Alias)

    # Branch, so no fusion at all
    tasks.append(Task("branch2", func, second.ref()))
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert len(result) == 6
    assert not any("-" in k for k in dsk)


def test_linear_fusion_intermediate_branch():
    tasks = [
        io := DataNode("foo", 1),
        second := Task("second", func, io.ref()),
        third := Task("third", func, second.ref()),
        other := Task("other", func, 1),
        Task("fourth", func, third.ref(), other.ref()),
    ]
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert "foo-second-third" in result
    assert isinstance(result["fourth"], Task)


def test_linear_fusion_two_branches():
    tasks = [
        left_one := DataNode("left_one", "a"),
        left_two := Task("left_two", func, left_one.ref()),
        right_one := DataNode("right_one", "a"),
        right_two := Task("right_two", func, right_one.ref()),
        middle := Task("middle", func, right_two.ref(), left_two.ref()),
        third := Task("third", func, middle.ref()),
        Task("fourth", func, third.ref()),
    ]
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert "left_one-left_two" in result
    assert "right_one-right_two" in result
    assert isinstance(result["right_two"], Alias)
    assert isinstance(result["left_two"], Alias)
    assert isinstance(result["fourth"], Alias)
    assert isinstance(result["third-fourth"], Task)


def test_linear_fusion_multiple_outputs():
    tasks = [
        first := DataNode("first", "a"),
        second := Task("second", func, first.ref()),
        third := Task("third", func, second.ref()),
        Task("fourth", func, third.ref()),
    ]
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth", "second"})
    assert "first-second" in result
    assert "second" in result
    assert isinstance(result["second"], Alias)
    assert isinstance(result["fourth"], Alias)
    assert isinstance(result["first-second"], Task)
    assert isinstance(result["third-fourth"], Task)


def test_nested_containers():
    t = List(Task("key-1", func, "a", "b"), Task("key-2", func, "c", "d"))
    assert t() == ["a-b", "c-d"]

    t = List(Task("key-1", func, "a", TaskRef("b")), Task("key-2", func, "c", "d"))
    assert t.dependencies == {"b"}
    assert t({"b": "b"}) == ["a-b", "c-d"]


def test_dict_class():
    t = Dict(k=Task("key-1", func, "a", "b"), v=Task("key-2", func, "c", "d"))
    assert not t.dependencies
    assert t() == {"k": "a-b", "v": "c-d"}
    t2 = Dict({"k": Task("key-1", func, "a", "b"), "v": Task("key-2", func, "c", "d")})
    assert t == t2
    assert t() == t2()
    t = Dict(k=Task("key-1", func, "a", TaskRef("b")), v=Task("key-2", func, "c", "d"))
    assert t.dependencies == {"b"}
    assert t({"b": "b"}) == {"k": "a-b", "v": "c-d"}

    t = Dict(
        {
            "k": Task("key-1", func, "a", TaskRef("b")),
            # Use a hashable key that isn't a string. This cannot be used as
            # kwargs, for instance
            ("v", 1): Task("key-2", func, "c", "d"),
        }
    )
    assert t.dependencies == {"b"}
    assert t({"b": "b"}) == {"k": "a-b", ("v", 1): "c-d"}

    t = Dict(
        k=Task("key-1", func, "a", "b"),
        v=Task("key-2", func, "c", "d"),
    )
    t2 = Dict(
        v=Task("key-2", func, "c", "d"),
        k=Task("key-1", func, "a", "b"),
    )
    assert t == t2
    assert tokenize(t) == tokenize(t2)

    d = Dict(
        [
            ["k", Task("key-1", func, "a", TaskRef("b"))],
            ["v", Task("key-2", func, TaskRef("c"), "d")],
        ]
    )
    assert d.dependencies == {"c", "b"}
    assert d({"b": "b", "c": "c"}) == {"k": "a-b", "v": "c-d"}

    d = Dict([("columns", ["a", "b"])])
    assert d() == {"columns": ["a", "b"]}

    d = Dict([["columns", ["a", "b"]]])
    assert d() == {"columns": ["a", "b"]}
    # Can be converted to a dict, e.g. also used as **kwargs
    assert isinstance(t, Mapping)
    assert dict(t) == {
        "k": Task("key-1", func, "a", "b"),
        "v": Task("key-2", func, "c", "d"),
    }
    assert len(t) == len(dict(t)) == 2

    def as_kwargs(**kwargs):
        return kwargs == {
            "k": Task("key-1", func, "a", "b"),
            "v": Task("key-2", func, "c", "d"),
        }

    assert as_kwargs(**t)


def test_block_io_fusion():

    class SubTask(Task):

        @property
        def block_fusion(self) -> bool:
            return True

    tasks = [
        io := DataNode("foo", 1),
        second := Task("second", func, io.ref()),
        third := SubTask("third", func, second.ref()),
        Task("fourth", func, third.ref()),
    ]
    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert "foo-second" in result
    assert isinstance(result["fourth"], Task)
    assert isinstance(result["third"], SubTask)


def test_data_producer():
    tasks = [
        io := DataNode("foo", 1),
        second := Task("second", func, io.ref(), _data_producer=True),
        third := Task("third", func, io.ref(), _data_producer=True),
        fourth := Task("fourth", func, second.ref()),
        fifth := Task("fifth", func, third.ref()),
        Task("sixth", func, fourth.ref(), fifth.ref()),
    ]

    dsk = {t.key: t for t in tasks}
    result = fuse_linear_task_spec(dsk, {"fourth"})
    assert "second-fourth" in result
    assert "third-fifth" in result
    assert result["second-fourth"].data_producer
    assert result["third-fifth"].data_producer
    assert not result["sixth"].data_producer
    assert result["foo"].data_producer


@pytest.mark.parametrize(
    "task_type, python_type",
    [
        (Tuple, tuple),
        (List, list),
        (Set, set),
    ],
)
def test_nested_containers_different_types(task_type, python_type):
    t = task_type(Task("key-1", func, "a", TaskRef("b")), Task("key-2", func, "c", "d"))
    assert t.dependencies == {"b"}
    assert t({"b": "b"}) == python_type(("a-b", "c-d"))


def test_substitute():
    t1 = Task("key-1", func, TaskRef("a"), "b")
    assert t1.substitute({"a": "a"}) is t1
    assert t1.substitute({"a": "a"}, key="foo") is not t1
    t2 = t1.substitute({"a": "c"})
    assert t2 is not t1
    assert t2.dependencies == {"c"}

    assert t1({"a": "a"}) == t2({"c": "a"})


def test_alias_task_ref_key():
    t = Alias(TaskRef("a"), "b")
    assert t.key == "a"
    assert isinstance(t.key, str)


def test_substitute_nested():
    def func(alist):
        return alist[0] + alist[1]["foo"]

    t1 = Task(
        "key-1",
        func,
        List(
            TaskRef("a"),
            Dict(
                {
                    "foo": TaskRef("b"),
                }
            ),
        ),
    )
    assert t1.dependencies == {"a", "b"}
    t2 = t1.substitute({"a": "c", "b": "d"})
    assert t2.dependencies == {"c", "d"}
    assert t1({"a": "a", "b": "b"}) == t2({"c": "a", "d": "b"})


@pytest.mark.parametrize("Container", [Dict, List, Set, Tuple])
def test_nested_containers_empty(Container):
    assert Container(Container.klass())() == Container.klass()


class MySubclass(Task):
    __slots__ = ("custom_kwarg_only",)

    def __init__(self, key, func, /, *args, custom_kwarg_only, **kwargs):
        self.custom_kwarg_only = custom_kwarg_only
        super().__init__(key, func, *args, **kwargs)


def test_substitute_subclasses():
    t = MySubclass("key", func, "a", TaskRef("b"), custom_kwarg_only="foo")
    t2 = t.substitute({"b": "c"})
    assert t2.custom_kwarg_only == "foo"
    assert t2({"a": "a", "c": "b"}) == t({"a": "a", "b": "b"})
