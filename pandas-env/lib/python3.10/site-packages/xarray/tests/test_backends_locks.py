from __future__ import annotations

import threading

from xarray.backends import locks


def test_threaded_lock() -> None:
    lock1 = locks._get_threaded_lock("foo")
    assert isinstance(lock1, type(threading.Lock()))
    lock2 = locks._get_threaded_lock("foo")
    assert lock1 is lock2

    lock3 = locks._get_threaded_lock("bar")
    assert lock1 is not lock3
