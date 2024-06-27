import unittest
import itertools
import numpy as np
from numba.misc.dummyarray import Array


class TestSlicing(unittest.TestCase):

    def assertSameContig(self, arr, nparr):
        attrs = 'C_CONTIGUOUS', 'F_CONTIGUOUS'
        for attr in attrs:
            if arr.flags[attr] != nparr.flags[attr]:
                if arr.size == 0 and nparr.size == 0:
                    # numpy <=1.7 bug that some empty array are contiguous and
                    # some are not
                    pass
                else:
                    self.fail("contiguous flag mismatch:\ngot=%s\nexpect=%s" %
                              (arr.flags, nparr.flags))

    #### 1D

    def test_slice0_1d(self):
        nparr = np.empty(4)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        self.assertSameContig(arr, nparr)
        xx = -2, -1, 0, 1, 2
        for x in xx:
            expect = nparr[x:]
            got = arr[x:]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    def test_slice1_1d(self):
        nparr = np.empty(4)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, -1, 0, 1, 2
        for x in xx:
            expect = nparr[:x]
            got = arr[:x]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    def test_slice2_1d(self):
        nparr = np.empty(4)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, -1, 0, 1, 2
        for x, y in itertools.product(xx, xx):
            expect = nparr[x:y]
            got = arr[x:y]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    #### 2D

    def test_slice0_2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, 0, 1, 2
        for x in xx:
            expect = nparr[x:]
            got = arr[x:]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

        for x, y in itertools.product(xx, xx):
            expect = nparr[x:, y:]
            got = arr[x:, y:]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    def test_slice1_2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, 0, 2
        for x in xx:
            expect = nparr[:x]
            got = arr[:x]
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)
            self.assertSameContig(got, expect)

        for x, y in itertools.product(xx, xx):
            expect = nparr[:x, :y]
            got = arr[:x, :y]
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)
            self.assertSameContig(got, expect)

    def test_slice2_2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, 0, 2
        for s, t, u, v in itertools.product(xx, xx, xx, xx):
            expect = nparr[s:t, u:v]
            got = arr[s:t, u:v]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

        for x, y in itertools.product(xx, xx):
            expect = nparr[s:t, u:v]
            got = arr[s:t, u:v]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    #### Strided

    def test_strided_1d(self):
        nparr = np.empty(4)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, -1, 1, 2
        for x in xx:
            expect = nparr[::x]
            got = arr[::x]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    def test_strided_2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, -1, 1, 2
        for a, b in itertools.product(xx, xx):
            expect = nparr[::a, ::b]
            got = arr[::a, ::b]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    def test_strided_3d(self):
        nparr = np.empty((4, 5, 6))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        xx = -2, -1, 1, 2
        for a, b, c in itertools.product(xx, xx, xx):
            expect = nparr[::a, ::b, ::c]
            got = arr[::a, ::b, ::c]
            self.assertSameContig(got, expect)
            self.assertEqual(got.shape, expect.shape)
            self.assertEqual(got.strides, expect.strides)

    def test_issue_2766(self):
        z = np.empty((1, 2, 3))
        z = np.transpose(z, axes=(2, 0, 1))
        arr = Array.from_desc(0, z.shape, z.strides, z.itemsize)
        self.assertEqual(z.flags['C_CONTIGUOUS'], arr.flags['C_CONTIGUOUS'])
        self.assertEqual(z.flags['F_CONTIGUOUS'], arr.flags['F_CONTIGUOUS'])


class TestReshape(unittest.TestCase):
    def test_reshape_2d2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(5, 4)
        got = arr.reshape(5, 4)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_2d1d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(5 * 4)
        got = arr.reshape(5 * 4)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_3d3d(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(5, 3, 4)
        got = arr.reshape(5, 3, 4)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_3d2d(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(3 * 4, 5)
        got = arr.reshape(3 * 4, 5)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_3d1d(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(3 * 4 * 5)
        got = arr.reshape(3 * 4 * 5)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_infer2d2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(-1, 4)
        got = arr.reshape(-1, 4)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_infer2d1d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(-1)
        got = arr.reshape(-1)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_infer3d3d(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(5, -1, 4)
        got = arr.reshape(5, -1, 4)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_infer3d2d(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(3, -1)
        got = arr.reshape(3, -1)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_infer3d1d(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        expect = nparr.reshape(-1)
        got = arr.reshape(-1)[0]
        self.assertEqual(got.shape, expect.shape)
        self.assertEqual(got.strides, expect.strides)

    def test_reshape_infer_two_unknowns(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)

        with self.assertRaises(ValueError) as raises:
            arr.reshape(-1, -1, 3)
        self.assertIn('can only specify one unknown dimension', str(raises.exception))

    def test_reshape_infer_invalid_shape(self):
        nparr = np.empty((3, 4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)

        with self.assertRaises(ValueError) as raises:
            arr.reshape(-1, 7)
        self.assertIn('cannot infer valid shape for unknown dimension', str(raises.exception))


class TestSqueeze(unittest.TestCase):
    def test_squeeze(self):
        nparr = np.empty((1, 2, 1, 4, 1, 3))
        arr = Array.from_desc(
            0, nparr.shape, nparr.strides, nparr.dtype.itemsize
        )
        def _assert_equal_shape_strides(arr1, arr2):
            self.assertEqual(arr1.shape, arr2.shape)
            self.assertEqual(arr1.strides, arr2.strides)
        _assert_equal_shape_strides(arr, nparr)
        _assert_equal_shape_strides(arr.squeeze()[0], nparr.squeeze())
        for axis in (0, 2, 4, (0, 2), (0, 4), (2, 4), (0, 2, 4)):
            _assert_equal_shape_strides(
                arr.squeeze(axis=axis)[0], nparr.squeeze(axis=axis)
            )

    def test_squeeze_invalid_axis(self):
        nparr = np.empty((1, 2, 1, 4, 1, 3))
        arr = Array.from_desc(
            0, nparr.shape, nparr.strides, nparr.dtype.itemsize
        )
        with self.assertRaises(ValueError):
            arr.squeeze(axis=1)
        with self.assertRaises(ValueError):
            arr.squeeze(axis=(2, 3))


class TestExtent(unittest.TestCase):
    def test_extent_1d(self):
        nparr = np.empty(4)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        s, e = arr.extent
        self.assertEqual(e - s, nparr.size * nparr.dtype.itemsize)

    def test_extent_2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        s, e = arr.extent
        self.assertEqual(e - s, nparr.size * nparr.dtype.itemsize)

    def test_extent_iter_1d(self):
        nparr = np.empty(4)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        [ext] = list(arr.iter_contiguous_extent())
        self.assertEqual(ext, arr.extent)

    def test_extent_iter_2d(self):
        nparr = np.empty((4, 5))
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)
        [ext] = list(arr.iter_contiguous_extent())
        self.assertEqual(ext, arr.extent)

        self.assertEqual(len(list(arr[::2].iter_contiguous_extent())), 2)


class TestIterate(unittest.TestCase):
    def test_for_loop(self):
        # for #4201
        N = 5
        nparr = np.empty(N)
        arr = Array.from_desc(0, nparr.shape, nparr.strides,
                              nparr.dtype.itemsize)

        x = 0  # just a placeholder
        # this loop should not raise AssertionError
        for val in arr:
            x = val


if __name__ == '__main__':
    unittest.main()

