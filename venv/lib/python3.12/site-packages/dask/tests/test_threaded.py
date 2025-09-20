from __future__ import annotations

import signal
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import ThreadPool
from time import sleep, time

import pytest

import dask
from dask.system import CPU_COUNT
from dask.threaded import get
from dask.utils_test import add, inc


def test_get():
    dsk = {"x": 1, "y": 2, "z": (inc, "x"), "w": (add, "z", "y")}
    assert get(dsk, "w") == 4
    assert get(dsk, ["w", "z"]) == (4, 2)


def test_nested_get():
    dsk = {"x": 1, "y": 2, "a": (add, "x", "y"), "b": (sum, ["x", "y"])}
    assert get(dsk, ["a", "b"]) == (3, 3)


def test_get_without_computation():
    dsk = {"x": 1}
    assert get(dsk, "x") == 1


def test_broken_callback():
    from dask.callbacks import Callback

    def _f_ok(*args, **kwargs):
        pass

    def _f_broken(*args, **kwargs):
        raise ValueError("my_exception")

    dsk = {"x": 1}

    with Callback(start=_f_broken, finish=_f_ok):
        with Callback(start=_f_ok, finish=_f_ok):
            with pytest.raises(ValueError, match="my_exception"):
                get(dsk, "x")


def bad(x):
    raise ValueError()


def test_exceptions_rise_to_top():
    dsk = {"x": 1, "y": (bad, "x")}
    pytest.raises(ValueError, lambda: get(dsk, "y"))


@pytest.mark.parametrize("pool_typ", [ThreadPool, ThreadPoolExecutor])
def test_reuse_pool(pool_typ):
    with pool_typ(CPU_COUNT) as pool:
        with dask.config.set(pool=pool):
            assert get({"x": (inc, 1)}, "x") == 2
            assert get({"x": (inc, 1)}, "x") == 2


@pytest.mark.parametrize("pool_typ", [ThreadPool, ThreadPoolExecutor])
def test_pool_kwarg(pool_typ):
    def f():
        sleep(0.01)
        return threading.get_ident()

    dsk = {("x", i): (f,) for i in range(30)}
    dsk["x"] = (len, (set, [("x", i) for i in range(len(dsk))]))

    with pool_typ(3) as pool:
        assert get(dsk, "x", pool=pool) == 3


def test_threaded_within_thread():
    L = []

    def f(i):
        result = get({"x": (lambda: i,)}, "x", num_workers=2)
        L.append(result)

    before = threading.active_count()

    for _ in range(20):
        t = threading.Thread(target=f, args=(1,))
        t.daemon = True
        t.start()
        t.join()
        assert L == [1]
        del L[:]

    start = time()  # wait for most threads to join
    while threading.active_count() > before + 10:
        sleep(0.01)
        assert time() < start + 5


def test_dont_spawn_too_many_threads():
    before = threading.active_count()

    dsk = {("x", i): (lambda i=i: i,) for i in range(10)}
    dsk["x"] = (sum, list(dsk))
    for _ in range(20):
        get(dsk, "x", num_workers=4)

    after = threading.active_count()

    assert after <= before + 8


def test_dont_spawn_too_many_threads_CPU_COUNT():
    before = threading.active_count()

    dsk = {("x", i): (lambda i=i: i,) for i in range(10)}
    dsk["x"] = (sum, list(dsk))
    for _ in range(20):
        get(dsk, "x")

    after = threading.active_count()

    assert after <= before + CPU_COUNT * 2


def test_thread_safety():
    def f(x):
        return 1

    dsk = {"x": (sleep, 0.05), "y": (f, "x")}

    L = []

    def test_f():
        L.append(get(dsk, "y"))

    threads = []
    for _ in range(20):
        t = threading.Thread(target=test_f)
        t.daemon = True
        t.start()
        threads.append(t)

    for thread in threads:
        thread.join()

    assert L == [1] * 20


def test_interrupt():
    # Windows implements `queue.get` using polling,
    # which means we can set an exception to interrupt the call to `get`.
    # Python 3 on other platforms requires sending SIGINT to the main thread.
    if sys.platform == "win32":
        from _thread import interrupt_main
    else:
        main_thread = threading.get_ident()

        def interrupt_main() -> None:
            signal.pthread_kill(main_thread, signal.SIGINT)

    in_clog_event = threading.Event()
    clog_event = threading.Event()

    def clog(in_clog_event: threading.Event, clog_event: threading.Event) -> None:
        in_clog_event.set()
        clog_event.wait()

    def interrupt(in_clog_event: threading.Event) -> None:
        in_clog_event.wait()
        interrupt_main()

    dsk = {("x", i): (clog, in_clog_event, clog_event) for i in range(20)}
    dsk["x"] = (len, list(dsk.keys()))

    interrupter = threading.Thread(target=interrupt, args=(in_clog_event,))
    interrupter.start()

    # Use explicitly created ThreadPoolExecutor to avoid leaking threads after
    # the KeyboardInterrupt
    with ThreadPoolExecutor(CPU_COUNT) as pool:
        with pytest.raises(KeyboardInterrupt):
            get(dsk, "x", pool=pool)
        clog_event.set()
    interrupter.join()
