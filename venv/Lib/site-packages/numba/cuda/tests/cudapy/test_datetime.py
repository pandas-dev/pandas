import numpy as np

from numba import cuda, vectorize, guvectorize
from numba.np.numpy_support import from_dtype
from numba.cuda.testing import CUDATestCase, skip_on_cudasim
import unittest


class TestCudaDateTime(CUDATestCase):
    def test_basic_datetime_kernel(self):
        @cuda.jit
        def foo(start, end, delta):
            for i in range(cuda.grid(1), delta.size, cuda.gridsize(1)):
                delta[i] = end[i] - start[i]

        arr1 = np.arange('2005-02', '2006-02', dtype='datetime64[D]')
        arr2 = arr1 + np.random.randint(0, 10000, arr1.size)
        delta = np.zeros_like(arr1, dtype='timedelta64[D]')

        foo[1, 32](arr1, arr2, delta)

        self.assertPreciseEqual(delta, arr2 - arr1)

    def test_scalar_datetime_kernel(self):
        @cuda.jit
        def foo(dates, target, delta, matches, outdelta):
            for i in range(cuda.grid(1), matches.size, cuda.gridsize(1)):
                matches[i] = dates[i] == target
                outdelta[i] = dates[i] - delta
        arr1 = np.arange('2005-02', '2006-02', dtype='datetime64[D]')
        target = arr1[5]           # datetime
        delta = arr1[6] - arr1[5]  # timedelta
        matches = np.zeros_like(arr1, dtype=np.bool_)
        outdelta = np.zeros_like(arr1, dtype='datetime64[D]')

        foo[1, 32](arr1, target, delta, matches, outdelta)
        where = matches.nonzero()

        self.assertEqual(list(where), [5])
        self.assertPreciseEqual(outdelta, arr1 - delta)

    @skip_on_cudasim('ufunc API unsupported in the simulator')
    def test_ufunc(self):
        datetime_t = from_dtype(np.dtype('datetime64[D]'))

        @vectorize([(datetime_t, datetime_t)], target='cuda')
        def timediff(start, end):
            return end - start

        arr1 = np.arange('2005-02', '2006-02', dtype='datetime64[D]')
        arr2 = arr1 + np.random.randint(0, 10000, arr1.size)

        delta = timediff(arr1, arr2)

        self.assertPreciseEqual(delta, arr2 - arr1)

    @skip_on_cudasim('ufunc API unsupported in the simulator')
    def test_gufunc(self):
        datetime_t = from_dtype(np.dtype('datetime64[D]'))
        timedelta_t = from_dtype(np.dtype('timedelta64[D]'))

        @guvectorize([(datetime_t, datetime_t, timedelta_t[:])], '(),()->()',
                     target='cuda')
        def timediff(start, end, out):
            out[0] = end - start

        arr1 = np.arange('2005-02', '2006-02', dtype='datetime64[D]')
        arr2 = arr1 + np.random.randint(0, 10000, arr1.size)

        delta = timediff(arr1, arr2)

        self.assertPreciseEqual(delta, arr2 - arr1)

    @skip_on_cudasim('no .copy_to_host() in the simulator')
    def test_datetime_view_as_int64(self):
        arr = np.arange('2005-02', '2006-02', dtype='datetime64[D]')
        darr = cuda.to_device(arr)
        viewed = darr.view(np.int64)
        self.assertPreciseEqual(arr.view(np.int64), viewed.copy_to_host())
        self.assertEqual(viewed.gpu_data, darr.gpu_data)

    @skip_on_cudasim('no .copy_to_host() in the simulator')
    def test_timedelta_view_as_int64(self):
        arr = np.arange('2005-02', '2006-02', dtype='datetime64[D]')
        arr = arr - (arr - 1)
        self.assertEqual(arr.dtype, np.dtype('timedelta64[D]'))
        darr = cuda.to_device(arr)
        viewed = darr.view(np.int64)
        self.assertPreciseEqual(arr.view(np.int64), viewed.copy_to_host())
        self.assertEqual(viewed.gpu_data, darr.gpu_data)


if __name__ == '__main__':
    unittest.main()
