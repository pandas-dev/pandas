from __future__ import annotations

import pickle
from collections import namedtuple

import pytest

from dask.core import (
    flatten,
    get,
    get_dependencies,
    get_deps,
    getcycle,
    ishashable,
    iskey,
    istask,
    literal,
    preorder_traversal,
    quote,
    subs,
    validate_key,
)
from dask.utils_test import GetFunctionTestMixin, add, inc


def contains(a, b):
    """

    >>> contains({'x': 1, 'y': 2}, {'x': 1})
    True
    >>> contains({'x': 1, 'y': 2}, {'z': 3})
    False
    """
    return all(a.get(k) == v for k, v in b.items())


def test_ishashable():
    class C:
        pass

    assert ishashable("x")
    assert ishashable(1)
    assert ishashable(C())
    assert ishashable((1, 2))
    assert not ishashable([1, 2])
    assert not ishashable({1: 2})


def test_iskey():
    class C:
        pass

    assert iskey("x")
    assert iskey(1)
    assert not iskey(C())  # Custom hashables can't be dask keys
    assert not iskey((C(),))
    assert iskey((1, 2))
    assert iskey(())
    assert iskey(("x",))
    assert not iskey([1, 2])
    assert not iskey({1, 2})
    assert not iskey({1: 2})


def test_iskey_numpy_types():
    np = pytest.importorskip("numpy")
    one = np.int64(1)
    assert not iskey(one)
    assert not iskey(("foo", one))


def test_validate_key():
    validate_key(1)
    validate_key(("x", 1))
    with pytest.raises(TypeError, match="Unexpected key type.*list"):
        validate_key(["x", 1])

    with pytest.raises(TypeError, match="unexpected key type at index=1"):
        validate_key((2, int))


def test_istask():
    assert istask((inc, 1))
    assert not istask(1)
    assert not istask((1, 2))
    f = namedtuple("f", ["x", "y"])
    assert not istask(f(sum, 2))


def test_preorder_traversal():
    t = (add, 1, 2)
    assert list(preorder_traversal(t)) == [add, 1, 2]
    t = (add, (add, 1, 2), (add, 3, 4))
    assert list(preorder_traversal(t)) == [add, add, 1, 2, add, 3, 4]
    t = (add, (sum, [1, 2]), 3)
    assert list(preorder_traversal(t)) == [add, sum, list, 1, 2, 3]


class TestGet(GetFunctionTestMixin):
    get = staticmethod(get)


def test_GetFunctionTestMixin_class():
    class TestCustomGetFail(GetFunctionTestMixin):
        get = staticmethod(lambda x, y: 1)

    custom_testget = TestCustomGetFail()
    pytest.raises(AssertionError, custom_testget.test_get)

    class TestCustomGetPass(GetFunctionTestMixin):
        get = staticmethod(get)

    custom_testget = TestCustomGetPass()
    custom_testget.test_get()


def test_get_dependencies_nested():
    dsk = {"x": 1, "y": 2, "z": (add, (inc, [["x"]]), "y")}

    assert get_dependencies(dsk, "z") == {"x", "y"}
    assert sorted(get_dependencies(dsk, "z", as_list=True)) == ["x", "y"]


def test_get_dependencies_empty():
    dsk = {"x": (inc,)}
    assert get_dependencies(dsk, "x") == set()
    assert get_dependencies(dsk, "x", as_list=True) == []


def test_get_dependencies_list():
    dsk = {"x": 1, "y": 2, "z": ["x", [(inc, "y")]]}
    assert get_dependencies(dsk, "z") == {"x", "y"}
    assert sorted(get_dependencies(dsk, "z", as_list=True)) == ["x", "y"]


def test_get_dependencies_task():
    dsk = {"x": 1, "y": 2, "z": ["x", [(inc, "y")]]}
    assert get_dependencies(dsk, task=(inc, "x")) == {"x"}
    assert get_dependencies(dsk, task=(inc, "x"), as_list=True) == ["x"]


def test_get_dependencies_nothing():
    with pytest.raises(ValueError):
        get_dependencies({})


def test_get_dependencies_many():
    dsk = {
        "a": [1, 2, 3],
        "b": "a",
        "c": [1, (inc, 1)],
        "d": [(sum, "c")],
        "e": ["a", "b", "zzz"],
        "f": [["a", "b"], 2, 3],
    }

    tasks = [dsk[k] for k in ("d", "f")]
    s = get_dependencies(dsk, task=tasks)
    assert s == {"a", "b", "c"}
    s = get_dependencies(dsk, task=tasks, as_list=True)
    assert sorted(s) == ["a", "b", "c"]

    s = get_dependencies(dsk, task=[])
    assert s == set()
    s = get_dependencies(dsk, task=[], as_list=True)
    assert s == []


def test_get_dependencies_task_none():
    # Regression test for https://github.com/dask/distributed/issues/2756
    dsk = {"foo": None}
    assert get_dependencies(dsk, task=dsk["foo"]) == set()


def test_get_deps():
    dsk = {
        "a": [1, 2, 3],
        "b": "a",
        "c": [1, (inc, 1)],
        "d": [(sum, "c")],
        "e": ["b", "zzz", "b"],
        "f": [["a", "b"], 2, 3],
    }
    dependencies, dependents = get_deps(dsk)
    assert dependencies == {
        "a": set(),
        "b": {"a"},
        "c": set(),
        "d": {"c"},
        "e": {"b"},
        "f": {"a", "b"},
    }
    assert dependents == {
        "a": {"b", "f"},
        "b": {"e", "f"},
        "c": {"d"},
        "d": set(),
        "e": set(),
        "f": set(),
    }


def test_flatten():
    assert list(flatten(())) == []
    assert list(flatten("foo")) == ["foo"]


def test_subs():
    assert subs((sum, [1, "x"]), "x", 2) == (sum, [1, 2])
    assert subs((sum, [1, ["x"]]), "x", 2) == (sum, [1, [2]])


class MutateOnEq:
    hit_eq = 0

    def __eq__(self, other):
        self.hit_eq += 1
        return False


def test_subs_no_key_data_eq():
    # Numpy throws a deprecation warning on bool(array == scalar), which
    # pollutes the terminal. This test checks that `subs` never tries to
    # compare keys (scalars) with values (which could be arrays)`subs` never
    # tries to compare keys (scalars) with values (which could be arrays).
    a = MutateOnEq()
    subs(a, "x", 1)
    assert a.hit_eq == 0
    subs((add, a, "x"), "x", 1)
    assert a.hit_eq == 0


def test_subs_with_unfriendly_eq():
    try:
        import numpy as np
    except ImportError:
        return
    else:
        task = (np.sum, np.array([1, 2]))
        assert (subs(task, (4, 5), 1) == task) is True

    class MyException(Exception):
        pass

    class F:
        def __eq__(self, other):
            raise MyException()

    task = F()
    assert subs(task, 1, 2) is task


def test_subs_with_surprisingly_friendly_eq():
    try:
        import pandas as pd
    except ImportError:
        return
    else:
        df = pd.DataFrame()
        assert subs(df, "x", 1) is df


def test_subs_arbitrary_key():
    key = (1.2, "foo", (3,))
    assert subs((id, key), key, 1) == (id, 1)


def test_quote():
    literals = [[1, 2, 3], (add, 1, 2), [1, [2, 3]], (add, 1, (add, 2, 3)), {"x": "x"}]

    for l in literals:
        assert get({"x": quote(l)}, "x") == l


def test_literal_serializable():
    l = literal((add, 1, 2))
    assert pickle.loads(pickle.dumps(l)).data == (add, 1, 2)


def test_getcycle():
    dsk = {
        0: [7, 13, 7, 9, 13, 3, 9, 18, 18, 17],
        1: [14, 14, 12, 1, 9, 16, 4, 5, 9, 8],
        2: [3, 1, 7, 7, 2, 0, 0, 6, 3, 2],
        3: [4, 8, 3, 14, 15, 19, 14, 1, 9, 1],
        4: [9, 13, 19, 4, 16, 8, 11, 1, 16, 1],
        5: [9, 8, 12, 13, 10, 13, 19, 3, 18, 18],
        6: [10, 2, 13, 16, 3, 12, 7, 16, 5, 17],
        7: [16, 8, 6, 4, 4, 10, 3, 1, 13, 10],
        8: [11, 4, 12, 10, 14, 6, 18, 15, 16, 12],
        9: [17, 18, 5, 16, 19, 16, 3, 6, 16, 18],
    }
    assert len(getcycle(dsk, list(dsk))) <= 3  # 7->6->7
    dsk = {
        0: [1, 27, 6, 25, 30, 16, 1, 1, 35, 17],
        1: [35, 22, 21, 31, 2, 28, 36, 0, 29, 29],
        2: [14, 27, 10, 1, 38, 18, 28, 28, 6, 0],
        3: [0, 37, 7, 16, 38, 28, 34, 13, 30, 10],
        4: [22, 22, 27, 13, 29, 36, 22, 9, 39, 19],
        5: [38, 7, 18, 17, 33, 5, 29, 11, 23, 30],
        6: [3, 30, 28, 38, 4, 13, 14, 27, 29, 38],
        7: [22, 27, 12, 2, 22, 16, 34, 15, 18, 16],
        8: [36, 21, 24, 22, 28, 38, 21, 2, 4, 24],
        9: [38, 32, 38, 7, 31, 34, 39, 20, 30, 18],
    }
    assert len(getcycle(dsk, list(dsk))) <= 4  # 0->1->2->0
