import threading
from numba import cuda
from numba.cuda.cudadrv.driver import driver
from numba.cuda.testing import unittest, ContextResettingTestCase
from queue import Queue


class TestResetDevice(ContextResettingTestCase):
    def test_reset_device(self):

        def newthread(exception_queue):
            try:
                devices = range(driver.get_device_count())
                for _ in range(2):
                    for d in devices:
                        cuda.select_device(d)
                        cuda.close()
            except Exception as e:
                exception_queue.put(e)

        # Do test on a separate thread so that we don't affect
        # the current context in the main thread.

        exception_queue = Queue()
        t = threading.Thread(target=newthread, args=(exception_queue,))
        t.start()
        t.join()

        exceptions = []
        while not exception_queue.empty():
            exceptions.append(exception_queue.get())
        self.assertEqual(exceptions, [])


if __name__ == '__main__':
    unittest.main()
