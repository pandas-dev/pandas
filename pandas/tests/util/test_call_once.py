import time

from pandas.util._call_once import call_once


class CallbackMock:
    # Should be MagicMock but uNitTesT.mOck is not allowed in pandas :(
    # (see #24648)

    def __init__(self):
        self.call_count = 0

    def __call__(self):
        self.call_count += 1

    def reset_mock(self):
        self.call_count = 0


def test_basic():
    callback1, callback2 = CallbackMock(), CallbackMock()
    with call_once(callback1):
        with call_once(callback1):
            with call_once(callback2):
                with call_once(callback1):
                    with call_once(callback2):
                        pass

    assert callback1.call_count == 1
    assert callback2.call_count == 1


def test_with_key():
    cb1, cb2 = CallbackMock(), CallbackMock()
    with call_once(cb1, key="callback"):
        with call_once(cb2, key="callback"):
            pass

    assert cb1.call_count == 1
    assert cb2.call_count == 0


def test_across_stack_frames():
    callback = CallbackMock()

    def f():
        with call_once(callback):
            pass

    def g():
        with call_once(callback):
            f()

    f()
    assert callback.call_count == 1
    callback.reset_mock()

    g()
    assert callback.call_count == 1


def test_concurrent_threading():
    import threading

    sleep_time = 0.01
    callback = CallbackMock()

    def run(initial_sleep=0):
        time.sleep(initial_sleep)
        with call_once(callback):
            with call_once(callback):
                time.sleep(2 * sleep_time)

    thread1 = threading.Thread(target=run)
    thread2 = threading.Thread(target=run, kwargs={"initial_sleep": sleep_time})
    thread2.start()
    thread1.start()
    thread1.join()
    thread2.join()
    assert callback.call_count == 2


def test_concurrent_asyncio():
    import asyncio

    callback = CallbackMock()

    async def task():
        with call_once(callback):
            with call_once(callback):
                await asyncio.sleep(0.01)

    async def main():
        await asyncio.gather(task(), task())

    asyncio.run(main())
    assert callback.call_count == 2
