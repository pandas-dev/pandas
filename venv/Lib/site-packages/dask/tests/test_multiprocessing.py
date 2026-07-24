from __future__ import annotations

import multiprocessing
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor
from operator import add

import pytest

import dask
from dask import compute, delayed
from dask.multiprocessing import _dumps, _loads, get, get_context, remote_exception
from dask.system import CPU_COUNT
from dask.utils_test import inc


def unrelated_function_global(a):
    np = pytest.importorskip("numpy")
    return np.array([a])


def my_small_function_global(a, b):
    return a + b


def test_pickle_globals():
    """Unrelated globals should not be included in serialized bytes"""
    b = _dumps(my_small_function_global)
    assert b"my_small_function_global" in b
    assert b"unrelated_function_global" not in b
    assert b"numpy" not in b


def test_pickle_locals():
    """Unrelated locals should not be included in serialized bytes"""
    np = pytest.importorskip("numpy")

    def unrelated_function_local(a):
        return np.array([a])

    def my_small_function_local(a, b):
        return a + b

    b = _dumps(my_small_function_local)
    assert b"my_small_function_global" not in b
    assert b"my_small_function_local" in b
    assert b"unrelated_function_local" not in b


@pytest.mark.skipif(pickle.HIGHEST_PROTOCOL < 5, reason="requires pickle protocol 5")
def test_out_of_band_pickling():
    """Test that out-of-band pickling works"""
    np = pytest.importorskip("numpy")
    pytest.importorskip("cloudpickle", minversion="1.3.0")

    a = np.arange(5)

    l = []
    b = _dumps(a, buffer_callback=l.append)
    assert len(l) == 1
    assert isinstance(l[0], pickle.PickleBuffer)
    assert memoryview(l[0]) == memoryview(a)

    a2 = _loads(b, buffers=l)
    assert np.all(a == a2)


def bad():
    raise ValueError("12345")


def test_errors_propagate():
    dsk = {"x": (bad,)}

    with pytest.raises(ValueError) as e:
        get(dsk, "x")
    assert "12345" in str(e.value)


def test_remote_exception():
    e = TypeError("hello")
    a = remote_exception(e, "traceback-body")
    b = remote_exception(e, "traceback-body")

    assert type(a) == type(b)
    assert isinstance(a, TypeError)
    assert "hello" in str(a)
    assert "Traceback" in str(a)
    assert "traceback-body" in str(a)


def test_lambda_with_cloudpickle():
    dsk = {"x": 2, "y": (lambda x: x + 1, "x")}
    assert get(dsk, "y") == 3


def lambda_result():
    return lambda x: x + 1


def test_lambda_results_with_cloudpickle():
    dsk = {"x": (lambda_result,)}
    f = get(dsk, "x")
    assert f(2) == 3


class NotUnpickleable:
    def __getstate__(self):
        return ()

    def __setstate__(self, state):
        raise ValueError("Can't unpickle me")


def test_unpicklable_args_generate_errors():
    a = NotUnpickleable()

    dsk = {"x": (bool, a)}

    with pytest.raises(ValueError):
        get(dsk, "x")

    dsk = {"x": (bool, "a"), "a": a}

    with pytest.raises(ValueError):
        get(dsk, "x")


@pytest.mark.parametrize("pool_typ", [multiprocessing.Pool, ProcessPoolExecutor])
def test_reuse_pool(pool_typ):
    with pool_typ(CPU_COUNT) as pool:
        with dask.config.set(pool=pool):
            assert get({"x": (inc, 1)}, "x") == 2
            assert get({"x": (inc, 1)}, "x") == 2


def test_dumps_loads():
    with dask.config.set(func_dumps=pickle.dumps, func_loads=pickle.loads):
        assert get({"x": 1, "y": (add, "x", 2)}, "y") == 3


def test_fuse_doesnt_clobber_intermediates():
    d = {"x": 1, "y": (inc, "x"), "z": (add, 10, "y")}
    assert get(d, ["y", "z"]) == (2, 12)


def test_optimize_graph_false():
    from dask.callbacks import Callback

    d = {"x": 1, "y": (inc, "x"), "z": (add, 10, "y")}
    keys = []
    with Callback(pretask=lambda key, *args: keys.append(key)):
        get(d, "z", optimize_graph=False)
    assert len(keys) == 2


def test_works_with_highlevel_graph():
    """Previously `dask.multiprocessing.get` would accidentally forward
    `HighLevelGraph` graphs through the dask optimization/scheduling routines,
    resulting in odd errors. One way to trigger this was to have a
    non-indexable object in a task. This is just a smoketest to ensure that
    things work properly even if `HighLevelGraph` objects get passed to
    `dask.multiprocessing.get`. See https://github.com/dask/dask/issues/7190.
    """

    class NoIndex:
        def __init__(self, x):
            self.x = x

        def __getitem__(self, key):
            raise Exception("Oh no!")

    x = delayed(lambda x: x)(NoIndex(1))
    (res,) = get(x.dask, x.__dask_keys__())
    assert isinstance(res, NoIndex)
    assert res.x == 1


@pytest.mark.parametrize("random", ["numpy", "random"])
def test_random_seeds(random):
    if random == "numpy":
        np = pytest.importorskip("numpy")
        random = np.random
    else:
        import random

    @delayed(pure=False)
    def f():
        return tuple(random.randint(0, 10000) for i in range(5))

    N = 10
    with dask.config.set(scheduler="processes"):
        (results,) = compute([f() for _ in range(N)])

    assert len(set(results)) == N


class global_:
    value = 0


def proc_init():
    global_.value = 1


@pytest.mark.parametrize(
    "scheduler, initializer, expected_results",
    [
        ("threading", None, [1] * 10),
        ("processes", None, [0] * 10),
        ("processes", proc_init, [1] * 10),
    ],
)
def test_process_initializer(scheduler, initializer, expected_results):
    @delayed(pure=False)
    def f():
        return global_.value

    global_.value = 1

    with dask.config.set(
        {"scheduler": scheduler, "multiprocessing.initializer": initializer}
    ):
        (results,) = compute([f() for _ in range(10)])
    assert results == expected_results

    (results2,) = compute(
        [f() for _ in range(10)],
        scheduler=scheduler,
        initializer=initializer,
    )
    assert results2 == expected_results


def check_for_pytest():
    """We check for spawn by ensuring subprocess doesn't have modules only
    parent process should have:
    """
    import sys

    return "FAKE_MODULE_FOR_TEST" in sys.modules


@pytest.mark.skipif(
    sys.platform == "win32", reason="Windows doesn't support different contexts"
)
def test_custom_context_used_python3_posix():
    """The 'multiprocessing.context' config is used to create the pool.

    We assume default is 'spawn', and therefore test for 'fork'.
    """
    # We check for 'fork' by ensuring subprocess doesn't have modules only
    # parent process should have:

    def check_for_pytest():
        import sys

        return "FAKE_MODULE_FOR_TEST" in sys.modules

    import sys

    sys.modules["FAKE_MODULE_FOR_TEST"] = 1
    try:
        with dask.config.set({"multiprocessing.context": "fork"}):
            result = get({"x": (check_for_pytest,)}, "x")
        assert result
    finally:
        del sys.modules["FAKE_MODULE_FOR_TEST"]


@pytest.mark.skipif(
    sys.platform == "win32", reason="Windows doesn't support different contexts"
)
def test_get_context_using_python3_posix():
    """get_context() respects configuration.

    If default context is changed this test will need to change too.
    """
    assert get_context() is multiprocessing.get_context("spawn")
    with dask.config.set({"multiprocessing.context": "forkserver"}):
        assert get_context() is multiprocessing.get_context("forkserver")
    with dask.config.set({"multiprocessing.context": "fork"}):
        assert get_context() is multiprocessing.get_context("fork")


@pytest.mark.skipif(sys.platform != "win32", reason="POSIX supports different contexts")
def test_custom_context_ignored_elsewhere():
    """On Windows, setting 'multiprocessing.context' doesn't explode.

    Presumption is it's not used since it's unsupported, but mostly we care about
    not breaking anything.
    """
    assert get({"x": (inc, 1)}, "x") == 2
    with pytest.warns(UserWarning):
        with dask.config.set({"multiprocessing.context": "forkserver"}):
            assert get({"x": (inc, 1)}, "x") == 2


@pytest.mark.skipif(sys.platform != "win32", reason="POSIX supports different contexts")
def test_get_context_always_default():
    """On Python 2/Windows, get_context() always returns same context."""
    assert get_context() is multiprocessing
    with pytest.warns(UserWarning):
        with dask.config.set({"multiprocessing.context": "forkserver"}):
            assert get_context() is multiprocessing
