#
# Test does not work on some cards.
#
import threading
from queue import Queue

import numpy as np
from numba import cuda
from numba.cuda.testing import unittest, ContextResettingTestCase


def newthread(exception_queue):
    try:
        cuda.select_device(0)
        stream = cuda.stream()
        A = np.arange(100)
        dA = cuda.to_device(A, stream=stream)
        stream.synchronize()
        del dA
        del stream
        cuda.close()
    except Exception as e:
        exception_queue.put(e)


class TestSelectDevice(ContextResettingTestCase):
    def test_select_device(self):
        exception_queue = Queue()
        for i in range(10):
            t = threading.Thread(target=newthread, args=(exception_queue,))
            t.start()
            t.join()

        exceptions = []
        while not exception_queue.empty():
            exceptions.append(exception_queue.get())
        self.assertEqual(exceptions, [])


if __name__ == '__main__':
    unittest.main()
