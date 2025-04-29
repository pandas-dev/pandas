from __future__ import annotations

import dataclasses
import datetime
import decimal
import operator
import pathlib
import pickle
import random
import subprocess
import sys
import textwrap
from concurrent.futures import ThreadPoolExecutor
from enum import Enum, Flag, IntEnum, IntFlag
from typing import Union

import cloudpickle
import pytest
from packaging.version import Version
from tlz import compose, curry, partial

import dask
from dask._compatibility import PY_VERSION
from dask.core import flatten, literal
from dask.tokenize import TokenizationError, normalize_token, tokenize
from dask.utils import tmpfile
from dask.utils_test import import_or_none

da = import_or_none("dask.array")
dd = import_or_none("dask.dataframe")
np = import_or_none("numpy")
sp = import_or_none("scipy.sparse")
pa = import_or_none("pyarrow")
pd = import_or_none("pandas")
numba = import_or_none("numba")


@pytest.fixture(autouse=True)
def check_clean_state():
    """Test that tokenize() and normalize_token() properly clean up state"""
    from dask.tokenize import _ENSURE_DETERMINISTIC, _SEEN

    assert not _SEEN
    with pytest.raises(LookupError):
        _ENSURE_DETERMINISTIC.get()
    yield
    assert not _SEEN
    with pytest.raises(LookupError):
        _ENSURE_DETERMINISTIC.get()


def check_tokenize(*args, **kwargs):
    with dask.config.set({"tokenize.ensure-deterministic": True}):
        before = tokenize(*args, **kwargs)

        # Test idempotency (the same object tokenizes to the same value)
        after = tokenize(*args, **kwargs)

        assert before == after, (args, kwargs)

        # Test same-interpreter determinism (two identical objects tokenize to the
        # same value as long as you do it on the same interpreter) We are not
        # particularly interested in a class that's never been pickled vs. one
        # that's already been pickled already (cloudpickle can introduce artifacts
        # on the first round-trip). We do care however about classes that have both
        # been through a serialization roundtrip at least once (not necessarily the
        # same amount of times).
        args2, kwargs2 = cloudpickle.loads(cloudpickle.dumps((args, kwargs)))
        args3, kwargs3 = cloudpickle.loads(cloudpickle.dumps((args, kwargs)))
        args3, kwargs3 = cloudpickle.loads(cloudpickle.dumps((args3, kwargs3)))

        tok2 = tokenize(*args2, **kwargs2)
        assert tok2 == before, (args, kwargs)

        tok3 = tokenize(*args3, **kwargs3)
        assert tok2 == tok3, (args, kwargs)

        # Skip: different interpreter determinism

    return before


def test_check_tokenize():
    # idempotent and deterministic
    check_tokenize(123)
    # returns tokenized value
    assert tokenize(123) == check_tokenize(123)

    # idempotent, but not deterministic
    class A:
        def __init__(self):
            self.tok = random.random()

        def __reduce__(self):
            return A, ()

        def __dask_tokenize__(self):
            return self.tok

    a = A()
    with pytest.raises(AssertionError):
        check_tokenize(a)

    # Not idempotent
    class B:
        def __dask_tokenize__(self):
            return random.random()

    b = B()
    with pytest.raises(AssertionError):
        check_tokenize(b)


def test_tokenize():
    a = (1, 2, 3)
    assert isinstance(tokenize(a), str)


@pytest.mark.skipif("not np")
def test_tokenize_scalar():
    assert check_tokenize(np.int64(1)) != check_tokenize(1)
    assert check_tokenize(np.int64(1)) != check_tokenize(np.int32(1))
    assert check_tokenize(np.int64(1)) != check_tokenize(np.uint32(1))
    assert check_tokenize(np.int64(1)) != check_tokenize("1")
    assert check_tokenize(np.int64(1)) != check_tokenize(np.float64(1))


@pytest.mark.skipif("not np")
def test_tokenize_numpy_array_consistent_on_values():
    assert check_tokenize(
        np.random.RandomState(1234).random_sample(1000)
    ) == check_tokenize(np.random.RandomState(1234).random_sample(1000))


@pytest.mark.skipif("not np")
def test_tokenize_numpy_array_supports_uneven_sizes():
    check_tokenize(np.random.random(7).astype(dtype="i2"))


@pytest.mark.skipif("not np")
def test_tokenize_discontiguous_numpy_array():
    arr = np.random.random(8)
    assert check_tokenize(arr[::2]) != check_tokenize(arr[::3])


@pytest.mark.skipif("not np")
def test_tokenize_numpy_datetime():
    check_tokenize(np.array(["2000-01-01T12:00:00"], dtype="M8[ns]"))


@pytest.mark.skipif("not np")
def test_tokenize_numpy_scalar():
    assert check_tokenize(np.array(1.0, dtype="f8")) == check_tokenize(
        np.array(1.0, dtype="f8")
    )
    assert check_tokenize(
        np.array([(1, 2)], dtype=[("a", "i4"), ("b", "i8")])[0]
    ) == check_tokenize(np.array([(1, 2)], dtype=[("a", "i4"), ("b", "i8")])[0])


@pytest.mark.skipif("not np")
def test_tokenize_numpy_scalar_string_rep():
    # Test tokenizing numpy scalars doesn't depend on their string representation
    with np.printoptions(formatter={"all": lambda x: "foo"}):
        assert check_tokenize(np.array(1)) != check_tokenize(np.array(2))


@pytest.mark.skipif("not np")
def test_tokenize_numpy_array_on_object_dtype():
    a = np.array(["a", "aa", "aaa"], dtype=object)
    assert check_tokenize(a) == check_tokenize(a)
    assert check_tokenize(np.array(["a", None, "aaa"], dtype=object)) == check_tokenize(
        np.array(["a", None, "aaa"], dtype=object)
    )
    assert check_tokenize(
        np.array([(1, "a"), (1, None), (1, "aaa")], dtype=object)
    ) == check_tokenize(np.array([(1, "a"), (1, None), (1, "aaa")], dtype=object))

    class NeedsCloudPickle:
        pass

    x = np.array(["a", None, NeedsCloudPickle()], dtype=object)
    check_tokenize(x)


@pytest.mark.skipif("not np")
def test_empty_numpy_array():
    arr = np.array([])
    assert arr.strides
    assert check_tokenize(arr) == check_tokenize(arr)
    arr2 = np.array([], dtype=np.int64())
    assert check_tokenize(arr) != check_tokenize(arr2)


@pytest.mark.skipif("not np")
def test_tokenize_numpy_memmap_offset(tmpdir):
    # Test two different memmaps into the same numpy file
    fn = str(tmpdir.join("demo_data"))

    with open(fn, "wb") as f:
        f.write(b"ashekwicht")

    with open(fn, "rb") as f:
        mmap1 = np.memmap(f, dtype=np.uint8, mode="r", offset=0, shape=5)
        mmap2 = np.memmap(f, dtype=np.uint8, mode="r", offset=5, shape=5)
        mmap3 = np.memmap(f, dtype=np.uint8, mode="r", offset=0, shape=5)
        assert check_tokenize(mmap1) == check_tokenize(mmap1)
        assert check_tokenize(mmap1) == check_tokenize(mmap3)
        assert check_tokenize(mmap1) != check_tokenize(mmap2)
        # also make sure that they tokenize correctly when taking sub-arrays
        assert check_tokenize(mmap1[1:-1]) == check_tokenize(mmap1[1:-1])
        assert check_tokenize(mmap1[1:-1]) == check_tokenize(mmap3[1:-1])
        assert check_tokenize(mmap1[1:2]) == check_tokenize(mmap3[1:2])
        assert check_tokenize(mmap1[1:2]) != check_tokenize(mmap1[1:3])
        assert check_tokenize(mmap1[1:2]) != check_tokenize(mmap3[1:3])
        sub1 = mmap1[1:-1]
        sub2 = mmap2[1:-1]
        assert check_tokenize(sub1) != check_tokenize(sub2)


@pytest.mark.skipif("not np")
def test_tokenize_numpy_memmap():
    with tmpfile(".npy") as fn:
        x1 = np.arange(5)
        np.save(fn, x1)
        y = check_tokenize(np.load(fn, mmap_mode="r"))

    with tmpfile(".npy") as fn:
        x2 = np.arange(5)
        np.save(fn, x2)
        z = check_tokenize(np.load(fn, mmap_mode="r"))

    assert check_tokenize(x1) == check_tokenize(x2)
    assert y == z

    with tmpfile(".npy") as fn:
        x = np.random.normal(size=(10, 10))
        np.save(fn, x)
        mm = np.load(fn, mmap_mode="r")
        mm2 = np.load(fn, mmap_mode="r")
        a = check_tokenize(mm[0, :])
        b = check_tokenize(mm[1, :])
        c = check_tokenize(mm[0:3, :])
        d = check_tokenize(mm[:, 0])
        assert len({a, b, c, d}) == 4
        assert check_tokenize(mm) == check_tokenize(mm2)
        assert check_tokenize(mm[1, :]) == check_tokenize(mm2[1, :])


@pytest.mark.skipif("not np")
def test_tokenize_numpy_memmap_no_filename():
    # GH 1562:
    with tmpfile(".npy") as fn1, tmpfile(".npy") as fn2:
        x = np.arange(5)
        np.save(fn1, x)
        np.save(fn2, x)

        a = np.load(fn1, mmap_mode="r")
        b = a + a
        assert check_tokenize(b) == check_tokenize(b)


@pytest.mark.skipif("not np")
def test_tokenize_numpy_ufunc():
    assert check_tokenize(np.sin) != check_tokenize(np.cos)

    np_ufunc = np.sin
    np_ufunc2 = np.cos
    assert isinstance(np_ufunc, np.ufunc)
    assert isinstance(np_ufunc2, np.ufunc)
    assert check_tokenize(np_ufunc) != check_tokenize(np_ufunc2)

    # for this we'll need the dask.array equivalent
    inc = da.ufunc.frompyfunc(lambda x: x + 1, 1, 1)
    inc2 = da.ufunc.frompyfunc(lambda x: x + 1, 1, 1)
    inc3 = da.ufunc.frompyfunc(lambda x: x + 2, 1, 1)
    assert check_tokenize(inc) != check_tokenize(inc2)
    assert check_tokenize(inc) != check_tokenize(inc3)


@pytest.mark.skipif("not np")
def test_normalize_numpy_ufunc_unserializable():
    # Make a ufunc that isn't in the numpy namespace and can't be serialized
    inc = np.frompyfunc(lambda x: x + 1, 1, 1)
    with dask.config.set({"tokenize.ensure-deterministic": False}):
        # Not idempotent
        assert tokenize(inc) != tokenize(inc)
        # You can call normalize_token directly
        assert normalize_token(inc) != normalize_token(inc)

    with dask.config.set({"tokenize.ensure-deterministic": True}):
        with pytest.raises(
            TokenizationError, match=r"Cannot tokenize.*dask\.array\.ufunc.*instead"
        ):
            tokenize(inc)

    # Test env override
    assert tokenize(inc, ensure_deterministic=False) != tokenize(
        inc, ensure_deterministic=False
    )
    with pytest.raises(TokenizationError, match="Cannot tokenize"):
        tokenize(inc, ensure_deterministic=True)


def test_normalize_object_unserializable():
    class C:
        def __reduce__(self):
            assert False

    c = C()

    with dask.config.set({"tokenize.ensure-deterministic": False}):
        # Not idempotent
        assert tokenize(c) != tokenize(c)
        # You can call normalize_token directly
        assert normalize_token(c) != normalize_token(c)

    with dask.config.set({"tokenize.ensure-deterministic": True}):
        with pytest.raises(
            TokenizationError, match="cannot be deterministically hashed"
        ):
            tokenize(c)

    # Test env override
    assert tokenize(c, ensure_deterministic=False) != tokenize(
        c, ensure_deterministic=False
    )
    with pytest.raises(TokenizationError, match="cannot be deterministically hashed"):
        tokenize(c, ensure_deterministic=True)


def test_tokenize_partial_func_args_kwargs_consistent():
    f = partial(f3, f2, c=f1)
    g = partial(f3, f2, c=f1)
    h = partial(f3, f2, c=5)
    assert check_tokenize(f) == check_tokenize(g)
    assert check_tokenize(f) != check_tokenize(h)


def test_normalize_base():
    for i in [
        1,
        1.1,
        "1",
        slice(1, 2, 3),
        decimal.Decimal("1.1"),
        datetime.date(2021, 6, 25),
        pathlib.PurePath("/this/that"),
    ]:
        assert normalize_token(i) is i


def test_tokenize_object():
    with dask.config.set({"tokenize.ensure-deterministic": False}):
        # object() tokenization is idempotent...
        o = object()
        assert tokenize(o) == tokenize(o)
        # ...but not deterministic
        assert tokenize(object()) != tokenize(object())

        # Two objects don't tokenize to the same token even if their pickle output is
        # identical. Stress id reuse by creating and dereferencing many objects in quick
        # succession.
        assert len({tokenize(object()) for _ in range(100)}) == 100

        # You can call normalize_token even if the _ensure_deterministic context
        # variable hasn't been set
        assert normalize_token(o) == normalize_token(o)

    with dask.config.set({"tokenize.ensure-deterministic": True}):
        with pytest.raises(TokenizationError, match="deterministic"):
            tokenize(o)
        with pytest.raises(TokenizationError, match="deterministic"):
            normalize_token(o)

    # Test env override
    assert tokenize(o, ensure_deterministic=False) == tokenize(
        o, ensure_deterministic=False
    )
    with pytest.raises(TokenizationError, match="deterministic"):
        tokenize(o, ensure_deterministic=True)


def nested_tokenize_ensure_deterministic():
    """Test that the ensure_deterministic override is not lost if tokenize() is
    called recursively
    """

    class C:
        def __dask_tokenize__(self):
            return tokenize(object())

    assert tokenize(C(), ensure_deterministic=False) != tokenize(
        C(), ensure_deterministic=False
    )
    with pytest.raises(TokenizationError):
        tokenize(C())


_GLOBAL = 1


def _local_functions():
    all_funcs = [
        lambda x: x,
        lambda x: x + 1,
        lambda y: y,
        lambda y: y + 1,
    ]
    a, b = all_funcs[:2]

    def func(x):
        return x

    def f2(x):  # Differs by function name
        return x

    # Suppress token differences determined by function name
    all_funcs += [func, f2]

    local_scope = 1

    def func():
        nonlocal local_scope
        local_scope += 1
        return a(local_scope)

    all_funcs.append(func)

    def func():
        global _GLOBAL
        _GLOBAL += 1
        return _GLOBAL

    all_funcs.append(func)

    # These functions differ only by the parameter defaults, which are also lambdas
    # Parameter defaults are lambdas

    def func(x, c=a):
        return c(x)

    all_funcs.append(func)

    def func(x, c=b):
        return c(x)

    all_funcs.append(func)

    # These functions differ only by the constants, which are also lambdas
    def func(x):
        c = lambda x: x + 2
        return c(x)

    all_funcs.append(func)

    def func(x):
        c = lambda x: x + 3
        return c(x)

    all_funcs.append(func)

    # These functions differ only by the imported names, which are also lambdas
    def func(x):
        c = a
        return c(x)

    all_funcs.append(func)

    def func(x):
        c = b
        return c(x)

    all_funcs.append(func)
    return all_funcs


class WithClassMethod:
    def f(self):
        pass

    @classmethod
    def g(cls):
        pass


_special_callables = [
    getattr,
    str.join,
    "foo".join,
    WithClassMethod.__str__,
    WithClassMethod().__str__,
    WithClassMethod.f,
    WithClassMethod().f,
    WithClassMethod.g,
]


@pytest.mark.parametrize("func", _local_functions())
def test_tokenize_local_functions(func):
    check_tokenize(func)


@pytest.mark.parametrize("func", _special_callables)
def test_tokenize_special_callables(func):
    check_tokenize(func)


def test_tokenize_functions_unique_token():
    all_funcs = _local_functions() + _special_callables
    tokens = [check_tokenize(func) for func in all_funcs]
    assert len(set(tokens)) == len(tokens)


@pytest.mark.xfail(reason="https://github.com/cloudpipe/cloudpickle/issues/453")
@pytest.mark.parametrize("instance", [False, True])
def test_tokenize_local_classes_from_different_contexts(instance):
    def f():
        class C:
            pass

        return C() if instance else C

    assert check_tokenize(f()) == check_tokenize(f())


def test_tokenize_local_functions_from_different_contexts():
    def f():
        def g():
            return 123

        return g

    assert check_tokenize(f()) == check_tokenize(f())


def f1(a, b, c=1):
    pass


def f2(a, b=1, c=2):
    pass


def f3(a):
    pass


def test_tokenize_callable():
    assert check_tokenize(f1) != check_tokenize(f2)


def test_tokenize_composite_functions():
    assert check_tokenize(partial(f2, b=2)) != check_tokenize(partial(f2, b=3))
    assert check_tokenize(partial(f1, b=2)) != check_tokenize(partial(f2, b=2))
    assert check_tokenize(compose(f2, f3)) != check_tokenize(compose(f2, f1))
    assert check_tokenize(curry(f2)) != check_tokenize(curry(f1))
    assert check_tokenize(curry(f2, b=1)) != check_tokenize(curry(f2, b=2))


@pytest.mark.skipif("not pd")
def test_tokenize_pandas():
    a = pd.DataFrame({"x": [1, 2, 3], "y": ["4", "asd", None]}, index=[1, 2, 3])
    b = pd.DataFrame({"x": [1, 2, 3], "y": ["4", "asd", None]}, index=[1, 2, 3])

    assert check_tokenize(a) == check_tokenize(b)
    b.index.name = "foo"
    assert check_tokenize(a) != check_tokenize(b)

    a = pd.DataFrame({"x": [1, 2, 3], "y": ["a", "b", "a"]})
    b = pd.DataFrame({"x": [1, 2, 3], "y": ["a", "b", "a"]})
    a["z"] = a.y.astype("category")
    assert check_tokenize(a) != check_tokenize(b)
    b["z"] = a.y.astype("category")
    assert check_tokenize(a) == check_tokenize(b)


@pytest.mark.skipif("not pd")
def test_tokenize_pandas_invalid_unicode():
    # see https://github.com/dask/dask/issues/2713
    df = pd.DataFrame(
        {"x\ud83d": [1, 2, 3], "y\ud83d": ["4", "asd\ud83d", None]}, index=[1, 2, 3]
    )
    check_tokenize(df)


@pytest.mark.skipif("not pd")
def test_tokenize_pandas_mixed_unicode_bytes():
    df = pd.DataFrame(
        {"ö".encode(): [1, 2, 3], "ö": ["ö", "ö".encode(), None]},
        index=[1, 2, 3],
    )
    check_tokenize(df)


@pytest.mark.skipif("not pd")
def test_tokenize_pandas_cloudpickle():
    class NeedsCloudPickle:
        # pickling not supported because it is a local class
        pass

    df = pd.DataFrame({"x": ["foo", None, NeedsCloudPickle()]})
    check_tokenize(df)


@pytest.mark.skipif("not dd")
def test_tokenize_pandas_extension_array():
    arrays = [
        pd.array([1, 0, None], dtype="Int64"),
        pd.array(["2000"], dtype="Period[D]"),
        pd.array([1, 0, 0], dtype="Sparse[int]"),
        pd.array([pd.Timestamp("2000")], dtype="datetime64[ns]"),
        pd.array([pd.Timestamp("2000", tz="CET")], dtype="datetime64[ns, CET]"),
        pd.array(
            ["a", "b"],
            dtype=pd.api.types.CategoricalDtype(["a", "b", "c"], ordered=False),
        ),
    ]

    arrays.extend(
        [
            pd.array(["a", "b", None], dtype="string"),
            pd.array([True, False, None], dtype="boolean"),
        ]
    )

    for arr in arrays:
        check_tokenize(arr)


@pytest.mark.skipif("not pd")
def test_tokenize_na():
    check_tokenize(pd.NA)


@pytest.mark.skipif("not pd")
def test_tokenize_offset():
    for offset in [
        pd.offsets.Second(1),
        pd.offsets.MonthBegin(2),
        pd.offsets.Day(1),
        pd.offsets.BQuarterEnd(2),
        pd.DateOffset(years=1),
        pd.DateOffset(months=7),
        pd.DateOffset(days=10),
    ]:
        check_tokenize(offset)


@pytest.mark.skipif("not pd")
def test_tokenize_pandas_index():
    idx = pd.Index(["a", "b"])
    check_tokenize(idx)

    idx = pd.MultiIndex.from_product([["a", "b"], [0, 1]])
    check_tokenize(idx)


def test_tokenize_kwargs():
    check_tokenize(5, x=1)
    assert check_tokenize(5) != check_tokenize(5, x=1)
    assert check_tokenize(5, x=1) != check_tokenize(5, x=2)
    assert check_tokenize(5, x=1) != check_tokenize(5, y=1)
    assert check_tokenize(5, foo="bar") != check_tokenize(5, {"foo": "bar"})


def test_tokenize_same_repr():
    class Foo:
        def __init__(self, x):
            self.x = x

        def __repr__(self):
            return "a foo"

    assert check_tokenize(Foo(1)) != check_tokenize(Foo(2))


def test_tokenize_slotted():
    class Foo:
        __slots__ = ("x",)

        def __init__(self, x):
            self.x = x

    assert check_tokenize(Foo(1)) != check_tokenize(Foo(2))


def test_tokenize_slotted_no_value():
    class Foo:
        __slots__ = ("x", "y")

        def __init__(self, x=None, y=None):
            if x is not None:
                self.x = x
            if y is not None:
                self.y = y

    assert check_tokenize(Foo(x=1)) != check_tokenize(Foo(y=1))
    check_tokenize(Foo())


def test_tokenize_slots_and_dict():
    class Foo:
        __slots__ = ("x",)

    class Bar(Foo):
        def __init__(self, x, y):
            self.x = x
            if y is not None:
                self.y = y

    assert Bar(1, 2).__dict__ == {"y": 2}

    tokens = [
        check_tokenize(Bar(1, 2)),
        check_tokenize(Bar(1, 3)),
        check_tokenize(Bar(1, None)),
        check_tokenize(Bar(2, 2)),
    ]
    assert len(set(tokens)) == len(tokens)


def test_tokenize_method():
    class Foo:
        def __init__(self, x):
            self.x = x

        def __dask_tokenize__(self):
            return self.x

        def hello(self):
            return "Hello world"

    a, b = Foo(1), Foo(2)
    assert check_tokenize(a) != check_tokenize(b)

    assert check_tokenize(a.hello) != check_tokenize(b.hello)

    # dispatch takes precedence
    before = check_tokenize(a)
    normalize_token.register(Foo, lambda self: self.x + 1)
    after = check_tokenize(a)
    assert before != after
    del normalize_token._lookup[Foo]


def test_tokenize_callable_class():
    class C:
        def __init__(self, x):
            self.x = x

        def __call__(self):
            return self.x

    class D(C):
        pass

    a, b, c = C(1), C(2), D(1)
    assert check_tokenize(a) != check_tokenize(b)
    assert check_tokenize(a) != check_tokenize(c)


def test_tokenize_callable_class_with_tokenize_method():
    """Always use ___dask_tokenize__ if present"""

    class C:
        def __init__(self, x, y):
            self.x = x
            self.y = y

        def __dask_tokenize__(self):
            return self.x

        def __call__(self):
            pass

    assert check_tokenize(C(1, 2)) == check_tokenize(C(1, 3))
    assert check_tokenize(C(1, 2)) != check_tokenize(C(2, 2))


class HasStaticMethods:
    def __init__(self, x):
        self.x = x

    def __dask_tokenize__(self):
        return normalize_token(type(self)), self.x

    def normal_method(self):
        pass

    @staticmethod
    def static_method():
        pass

    @classmethod
    def class_method(cls):
        pass


class HasStaticMethods2(HasStaticMethods):
    pass


def test_staticmethods():
    a, b, c = HasStaticMethods(1), HasStaticMethods(2), HasStaticMethods2(1)
    # These invoke HasStaticMethods.__dask_tokenize__()
    assert check_tokenize(a.normal_method) != check_tokenize(b.normal_method)
    assert check_tokenize(a.normal_method) != check_tokenize(c.normal_method)
    # These don't
    assert check_tokenize(a.static_method) == check_tokenize(b.static_method)
    assert check_tokenize(a.static_method) == check_tokenize(c.static_method)
    assert check_tokenize(a.class_method) == check_tokenize(b.class_method)
    assert check_tokenize(a.class_method) != check_tokenize(c.class_method)


def test_tokenize_sequences():
    assert check_tokenize([1]) != check_tokenize([2])
    assert check_tokenize([1]) != check_tokenize((1,))
    assert check_tokenize([1]) == check_tokenize([1])

    # You can call normalize_token directly.
    x = (1, 2)
    y = [x, x, [x, (2, 3)]]
    assert normalize_token(y)


def test_nested_tokenize_seen():
    """Test that calling tokenize() recursively doesn't alter the output due to
    memoization of already-seen objects
    """
    o = [1, 2, 3]

    class C:
        def __init__(self, x):
            self.x = x
            self.tok = None

        def __dask_tokenize__(self):
            if not self.tok:
                self.tok = tokenize(self.x)
            return self.tok

    c1, c2 = C(o), C(o)
    check_tokenize(o, c1, o)
    assert c1.tok
    assert check_tokenize(c1) == check_tokenize(c2)


def test_tokenize_dict():
    # Insertion order is ignored. Keys can be an unsortable mix of types.
    assert check_tokenize({"x": 1, 1: "x"}) == check_tokenize({1: "x", "x": 1})
    assert check_tokenize({"x": 1, 1: "x"}) != check_tokenize({"x": 1, 2: "x"})
    assert check_tokenize({"x": 1, 1: "x"}) != check_tokenize({"x": 2, 1: "x"})


def test_tokenize_set():
    assert check_tokenize({1, 2, "x", (1, "x")}) == check_tokenize(
        {2, "x", (1, "x"), 1}
    )


def test_tokenize_ordered_dict():
    from collections import OrderedDict

    a = OrderedDict([("a", 1), ("b", 2)])
    b = OrderedDict([("a", 1), ("b", 2)])
    c = OrderedDict([("b", 2), ("a", 1)])

    assert check_tokenize(a) == check_tokenize(b)
    assert check_tokenize(a) != check_tokenize(c)


def test_tokenize_dict_doesnt_call_str_on_values():
    class C:
        def __dask_tokenize__(self):
            return "C"

        def __repr__(self):
            assert False

    check_tokenize({1: C(), "2": C()})


def test_tokenize_sorts_dict_before_seen_map():
    v = (1, 2, 3)
    d1 = {1: v, 2: v}
    d2 = {2: v, 1: v}
    assert check_tokenize(d1) == check_tokenize(d2)


def test_tokenize_sorts_set_before_seen_map():
    v = (1, 2, 3)
    s1 = {(i, v) for i in range(100)}
    s2 = {(i, v) for i in reversed(range(100))}
    assert check_tokenize(s1) == check_tokenize(s2)


def test_tokenize_timedelta():
    assert check_tokenize(datetime.timedelta(days=1)) == check_tokenize(
        datetime.timedelta(days=1)
    )
    assert check_tokenize(datetime.timedelta(days=1)) != check_tokenize(
        datetime.timedelta(days=2)
    )


@pytest.mark.parametrize("enum_type", [Enum, IntEnum, IntFlag, Flag])
def test_tokenize_enum(enum_type):
    class Color(enum_type):
        RED = 1
        BLUE = 2

    assert check_tokenize(Color.RED) == check_tokenize(Color.RED)
    assert check_tokenize(Color.RED) != check_tokenize(Color.BLUE)


@dataclasses.dataclass
class ADataClass:
    a: int


@dataclasses.dataclass
class BDataClass:
    a: float


@dataclasses.dataclass
class NoValueDataClass:
    a: int = dataclasses.field(init=False)


class GlobalClass:
    def __init__(self, val) -> None:
        self.val = val


def test_local_objects():
    class LocalType:
        foo = "bar"

    class LocalReducible:
        def __reduce__(self):
            return LocalReducible, ()

    class LocalDaskTokenize:
        def __dask_tokenize__(self):
            return "foo"

    class LocalChild(GlobalClass):
        pass

    check_tokenize(GlobalClass(1))
    assert check_tokenize(GlobalClass(1)) != check_tokenize(GlobalClass(2))
    check_tokenize(LocalType())
    check_tokenize(LocalChild(1))

    assert check_tokenize(LocalDaskTokenize()) != check_tokenize(LocalReducible())


@pytest.mark.skipif(
    PY_VERSION >= Version("3.13"), reason="https://github.com/dask/dask/issues/11457"
)
def test_tokenize_dataclass():
    a1 = ADataClass(1)
    a2 = ADataClass(2)
    check_tokenize(a1)
    assert check_tokenize(a1) != check_tokenize(a2)

    # Same field names and values, but dataclass types are different
    b1 = BDataClass(1)
    assert check_tokenize(ADataClass) != check_tokenize(BDataClass)
    assert check_tokenize(a1) != check_tokenize(b1)

    class SubA(ADataClass):
        pass

    assert dataclasses.is_dataclass(SubA)
    assert check_tokenize(ADataClass) != check_tokenize(SubA)
    assert check_tokenize(SubA(1)) != check_tokenize(a1)

    # Same name, same values, new definition: tokenize differently
    ADataClassRedefinedDifferently = dataclasses.make_dataclass(
        "ADataClass", [("a", Union[int, str])]
    )
    assert check_tokenize(a1) != check_tokenize(ADataClassRedefinedDifferently(1))

    # Dataclass with unpopulated value
    nv = NoValueDataClass()
    check_tokenize(nv)


@pytest.mark.parametrize(
    "other",
    [
        (1, 10, 2),  # Different start
        (5, 15, 2),  # Different stop
        (5, 10, 1),  # Different step
    ],
)
def test_tokenize_range(other):
    assert check_tokenize(range(5, 10, 2)) != check_tokenize(range(*other))


@pytest.mark.skipif("not np")
def test_tokenize_numpy_array():
    x = np.arange(2000)  # long enough to drop information in repr
    y = np.arange(2000)
    y[1000] = 0  # middle isn't printed in repr
    assert check_tokenize([x]) != check_tokenize([y])


@pytest.mark.skipif("not np")
def test_tokenize_object_array_with_nans():
    a = np.array(["foo", "Jos\xe9", np.nan], dtype="O")
    check_tokenize(a)


@pytest.mark.parametrize(
    "x", [1, True, "a", b"a", 1.0, 1j, 1.0j, [], (), {}, None, str, int]
)
def test_tokenize_base_types(x):
    check_tokenize(x)


def test_tokenize_literal():
    assert check_tokenize(literal(["x", 1])) != check_tokenize(literal(["x", 2]))


@pytest.mark.skipif("not np")
@pytest.mark.filterwarnings("ignore:the matrix:PendingDeprecationWarning")
def test_tokenize_numpy_matrix():
    rng = np.random.RandomState(1234)
    a = np.asmatrix(rng.rand(100))
    b = a.copy()
    assert check_tokenize(a) == check_tokenize(b)

    b[:10] = 1
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.skipif("not sp")
@pytest.mark.parametrize("cls_name", ("dok",))
def test_tokenize_dense_sparse_array(cls_name):
    rng = np.random.RandomState(1234)

    a = sp.rand(10, 100, random_state=rng).asformat(cls_name)
    b = a.copy()

    assert check_tokenize(a) == check_tokenize(b)

    # modifying the data values
    if hasattr(b, "data"):
        b.data[:10] = 1
    elif cls_name == "dok":
        b[3, 3] = 1
    else:
        raise ValueError
    check_tokenize(b)
    assert check_tokenize(a) != check_tokenize(b)

    # modifying the data indices
    b = a.copy().asformat("coo")
    b.row[:10] = np.arange(10)
    b = b.asformat(cls_name)
    assert check_tokenize(a) != check_tokenize(b)


def test_tokenize_circular_recursion():
    a = [1, 2]
    a[0] = a

    # Test that tokenization doesn't stop as soon as you hit the circular recursion
    b = [1, 3]
    b[0] = b

    assert check_tokenize(a) != check_tokenize(b)

    # Different circular recursions tokenize differently
    c = [[], []]
    c[0].append(c[0])
    c[1].append(c[1])

    d = [[], []]
    d[0].append(d[1])
    d[1].append(d[0])
    assert check_tokenize(c) != check_tokenize(d)

    # For dicts, the dict itself is not passed to _normalize_seq_func
    e = {}
    e[0] = e
    check_tokenize(e)


@pytest.mark.parametrize(
    "other",
    [
        (2002, 6, 25),  # Different year
        (2021, 7, 25),  # Different month
        (2021, 6, 26),  # Different day
    ],
)
def test_tokenize_datetime_date(other):
    a = datetime.date(2021, 6, 25)
    b = datetime.date(*other)
    assert check_tokenize(a) != check_tokenize(b)


def test_tokenize_datetime_time():
    # Same time
    check_tokenize(datetime.time(1, 2, 3, 4, datetime.timezone.utc))
    check_tokenize(datetime.time(1, 2, 3, 4))
    check_tokenize(datetime.time(1, 2, 3))
    check_tokenize(datetime.time(1, 2))
    # Different hour
    assert check_tokenize(
        datetime.time(1, 2, 3, 4, datetime.timezone.utc)
    ) != check_tokenize(datetime.time(2, 2, 3, 4, datetime.timezone.utc))
    # Different minute
    assert check_tokenize(
        datetime.time(1, 2, 3, 4, datetime.timezone.utc)
    ) != check_tokenize(datetime.time(1, 3, 3, 4, datetime.timezone.utc))
    # Different second
    assert check_tokenize(
        datetime.time(1, 2, 3, 4, datetime.timezone.utc)
    ) != check_tokenize(datetime.time(1, 2, 4, 4, datetime.timezone.utc))
    # Different micros
    assert check_tokenize(
        datetime.time(1, 2, 3, 4, datetime.timezone.utc)
    ) != check_tokenize(datetime.time(1, 2, 3, 5, datetime.timezone.utc))
    # Different tz
    assert check_tokenize(
        datetime.time(1, 2, 3, 4, datetime.timezone.utc)
    ) != check_tokenize(datetime.time(1, 2, 3, 4))


def test_tokenize_datetime_datetime():
    # Same datetime
    required = [1, 2, 3]  # year, month, day
    optional = [4, 5, 6, 7, datetime.timezone.utc]
    for i in range(len(optional) + 1):
        args = required + optional[:i]
        check_tokenize(datetime.datetime(*args))

    # Different year
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(2, 2, 3, 4, 5, 6, 7, datetime.timezone.utc))
    # Different month
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 1, 3, 4, 5, 6, 7, datetime.timezone.utc))
    # Different day
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 2, 2, 4, 5, 6, 7, datetime.timezone.utc))
    # Different hour
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 2, 3, 3, 5, 6, 7, datetime.timezone.utc))
    # Different minute
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 2, 3, 4, 4, 6, 7, datetime.timezone.utc))
    # Different second
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 2, 3, 4, 5, 5, 7, datetime.timezone.utc))
    # Different micros
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 2, 3, 4, 5, 6, 6, datetime.timezone.utc))
    # Different tz
    assert check_tokenize(
        datetime.datetime(1, 2, 3, 4, 5, 6, 7, datetime.timezone.utc)
    ) != check_tokenize(datetime.datetime(1, 2, 3, 4, 5, 6, 7, None))


def test_tokenize_functions_main():
    script = """
    def inc(x):
        return x + 1

    inc2 = inc
    def sum(x, y):
        return x + y

    from dask.base import tokenize
    assert tokenize(inc) != tokenize(sum)
    # That this is an alias shouldn't matter
    assert tokenize(inc) == tokenize(inc2)

    def inc(x):
        return x + 1

    assert tokenize(inc2) != tokenize(inc)

    def inc(y):
        return y + 1

    assert tokenize(inc2) != tokenize(inc)

    def inc(x):
        y = x
        return y + 1

    assert tokenize(inc2) != tokenize(inc)
    """
    proc = subprocess.run([sys.executable, "-c", textwrap.dedent(script)])
    proc.check_returncode()


def test_tokenize_dataclass_field_no_repr():
    A = dataclasses.make_dataclass(
        "A",
        [("param", float, dataclasses.field(repr=False))],
        namespace={"__dask_tokenize__": lambda self: self.param},
    )

    a1, a2 = A(1), A(2)

    assert check_tokenize(a1) != check_tokenize(a2)


def test_tokenize_operator():
    """Top-level functions in the operator module have a __self__ attribute, which is
    the module itself
    """
    assert check_tokenize(operator.add) != check_tokenize(operator.mul)


def test_tokenize_random_state():
    a = random.Random(123)
    b = random.Random(123)
    c = random.Random(456)
    assert check_tokenize(a) == check_tokenize(b)
    assert check_tokenize(a) != check_tokenize(c)
    a.random()
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.skipif("not np")
def test_tokenize_random_state_numpy():
    a = np.random.RandomState(123)
    b = np.random.RandomState(123)
    c = np.random.RandomState(456)
    assert check_tokenize(a) == check_tokenize(b)
    assert check_tokenize(a) != check_tokenize(c)
    a.random()
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.parametrize(
    "module",
    ["random", pytest.param("np.random", marks=pytest.mark.skipif("not np"))],
)
def test_tokenize_random_functions(module):
    """random.random() and other methods of the global random state do not compare as
    equal to themselves after a pickle roundtrip"""
    module = eval(module)
    module.seed(2)

    a = module.random
    b = pickle.loads(pickle.dumps(a))
    assert check_tokenize(a) == check_tokenize(b)

    # Drawing elements or reseeding changes the global state
    a()
    c = pickle.loads(pickle.dumps(a))
    assert check_tokenize(a) == check_tokenize(c)
    assert check_tokenize(a) != check_tokenize(b)

    module.seed(123)
    d = pickle.loads(pickle.dumps(a))
    assert check_tokenize(a) == check_tokenize(d)
    assert check_tokenize(a) != check_tokenize(c)


def test_tokenize_random_functions_with_state():
    a = random.Random(123).random
    b = random.Random(456).random
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.skipif("not np")
def test_tokenize_random_functions_with_state_numpy():
    a = np.random.RandomState(123).random
    b = np.random.RandomState(456).random
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.skipif("not pa")
def test_tokenize_pyarrow_datatypes_simple():
    a = pa.int64()
    b = pa.float64()
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.skipif("not pa")
def test_tokenize_pyarrow_datatypes_complex():
    a = pa.struct({"x": pa.int32(), "y": pa.string()})
    b = pa.struct({"x": pa.float64(), "y": pa.int16()})
    assert check_tokenize(a) != check_tokenize(b)


@pytest.mark.skipif("not pa")
def test_pyarrow_table():
    a = pa.table({"x": [1, 2, 3], "y": ["a", "b", "c"]})
    b = pa.table({"x": [1, 2, 3], "y": ["a", "b", "c"]})
    c = pa.table({"x": [1, 2, 3], "y": ["a", "b", "d"]})
    assert check_tokenize(a) == check_tokenize(b)
    assert check_tokenize(a) != check_tokenize(c)


@pytest.mark.skipif("not np")
def test_tokenize_opaque_object_with_buffers():
    # pickle will extract PickleBuffer objects out of this
    class C:
        def __init__(self, x):
            self.x = np.array(x)

    assert check_tokenize(C([1, 2])) != check_tokenize(C([1, 3]))


if not numba:

    class NumbaDummy:
        def __bool__(self):
            return False

        def _dummy_decorator(self, *args, **kwargs):
            def wrapper(func):
                return func

            return wrapper

        jit = vectorize = guvectorize = _dummy_decorator

    numba = NumbaDummy()


@numba.jit(nopython=True)
def numba_jit(x, y):
    return x + y


@numba.jit("f8(f8, f8)", nopython=True)
def numba_jit_with_signature(x, y):
    return x + y


@numba.vectorize(nopython=True)
def numba_vectorize(x, y):
    return x + y


@numba.vectorize("f8(f8, f8)", nopython=True)
def numba_vectorize_with_signature(x, y):
    return x + y


@numba.guvectorize(["f8,f8,f8[:]"], "(),()->()")
def numba_guvectorize(x, y, out):
    out[0] = x + y


all_numba_funcs = [
    numba_jit,
    numba_jit_with_signature,
    numba_vectorize,
    numba_vectorize_with_signature,
    numba_guvectorize,
]


@pytest.mark.skipif("not numba")
@pytest.mark.parametrize("func", all_numba_funcs)
def test_tokenize_numba(func):
    assert func(1, 2) == 3
    check_tokenize(func)
    for func in all_numba_funcs:
        tokens = normalize_token(func)

        # Ensure that we attempt to tokenize it instead of dumping it into pickle
        assert isinstance(tokens, tuple)
        assert isinstance(tokens[1], tuple)


@pytest.mark.skipif("not numba")
def test_tokenize_numba_unique_token():
    tokens = [check_tokenize(func) for func in all_numba_funcs]
    assert len(tokens) == len(set(tokens))


@pytest.mark.skipif("not numba")
def test_numba_local():
    @numba.jit(nopython=True)
    def local_jit(x, y):
        return x + y

    @numba.jit("f8(f8, f8)", nopython=True)
    def local_jit_with_signature(x, y):
        return x + y

    @numba.vectorize(nopython=True)
    def local_vectorize(x, y):
        return x + y

    @numba.vectorize("f8(f8, f8)", nopython=True)
    def local_vectorize_with_signature(x, y):
        return x + y

    @numba.guvectorize(["f8,f8,f8[:]"], "(),()->()")
    def local_guvectorize(x, y, out):
        out[0] = x + y

    all_funcs = [
        local_jit,
        local_jit_with_signature,
        local_vectorize,
        local_vectorize_with_signature,
        local_guvectorize,
    ]
    tokens = [check_tokenize(func) for func in all_funcs]
    assert len(tokens) == len(set(tokens))


@pytest.mark.skipif("not np")
def test_tokenize_np_dtype():
    arr = np.array([1, 2, 3], dtype=np.int64)
    arr2 = np.array([1, 2, 3], dtype=np.int32)
    assert check_tokenize(arr.dtype) != check_tokenize(arr2.dtype)


@pytest.mark.skipif("not pd")
def test_tokenize_pandas_arrow_strings():
    ser = pd.Series(["a", "b"], dtype="string[pyarrow]")
    check_tokenize(ser)
    tokens = normalize_token(ser)
    # Maybe a little brittle but will do for now
    assert any(str(tok) == "string" for tok in flatten(tokens))


class APickleable:
    counter = 0

    def __reduce__(self):
        APickleable.counter += 1
        return APickleable, ()


def test_normalize_pickle():
    a = APickleable()
    tokenize(a)
    # We're pickling multiple times because pickle is caching things on
    # instances such that subsequent pickles can yield different results.
    # For a trivial object like this, this should only happen twice
    assert APickleable.counter <= 2


def test_tokenize_recursive_respects_ensure_deterministic():
    class Foo:
        def __dask_tokenize__(self):
            return tokenize(object())

    with pytest.raises(RuntimeError):
        tokenize(Foo(), ensure_deterministic=True)


def test_tokenize_nested_sequence_thread_safe():
    nested_list = []
    for ix in range(100):
        nested_list = [ix, nested_list]

    nested_dict = {}
    for ix in range(50):
        nested_dict = {ix: nested_dict}

    normalize_token(nested_list)
    with ThreadPoolExecutor() as pool:
        futures = [pool.submit(normalize_token, nested_list) for _ in range(1000)]
        assert len({f.result() for f in futures}) == 1
        futures = [pool.submit(normalize_token, nested_dict) for _ in range(1000)]
        assert len({f.result() for f in futures}) == 1
