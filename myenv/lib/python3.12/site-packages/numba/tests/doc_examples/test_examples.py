# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import sys
import unittest

from numba.tests.support import TestCase, captured_stdout
from numba.core.config import IS_WIN32
from numba.np.numpy_support import numpy_version


class MatplotlibBlocker:
    '''Blocks the import of matplotlib, so that doc examples that attempt to
    plot the output don't result in plots popping up and blocking testing.'''

    def find_spec(self, fullname, path, target=None):
        if fullname == 'matplotlib':
            msg = 'Blocked import of matplotlib for test suite run'
            raise ImportError(msg)


class DocsExamplesTest(TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mpl_blocker = MatplotlibBlocker()

    def setUp(self):
        sys.meta_path.insert(0, self._mpl_blocker)

    def tearDown(self):
        sys.meta_path.remove(self._mpl_blocker)

    def test_mandelbrot(self):
        with captured_stdout():
            # magictoken.ex_mandelbrot.begin
            from timeit import default_timer as timer
            try:
                from matplotlib.pylab import imshow, show
                have_mpl = True
            except ImportError:
                have_mpl = False
            import numpy as np
            from numba import jit

            @jit(nopython=True)
            def mandel(x, y, max_iters):
                """
                Given the real and imaginary parts of a complex number,
                determine if it is a candidate for membership in the Mandelbrot
                set given a fixed number of iterations.
                """
                i = 0
                c = complex(x,y)
                z = 0.0j
                for i in range(max_iters):
                    z = z * z + c
                    if (z.real * z.real + z.imag * z.imag) >= 4:
                        return i

                return 255

            @jit(nopython=True)
            def create_fractal(min_x, max_x, min_y, max_y, image, iters):
                height = image.shape[0]
                width = image.shape[1]

                pixel_size_x = (max_x - min_x) / width
                pixel_size_y = (max_y - min_y) / height
                for x in range(width):
                    real = min_x + x * pixel_size_x
                    for y in range(height):
                        imag = min_y + y * pixel_size_y
                        color = mandel(real, imag, iters)
                        image[y, x] = color

                return image

            image = np.zeros((500 * 2, 750 * 2), dtype=np.uint8)
            s = timer()
            create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20)
            e = timer()
            print(e - s)
            if have_mpl:
                imshow(image)
                show()
            # magictoken.ex_mandelbrot.end

    def test_moving_average(self):
        with captured_stdout():
            # magictoken.ex_moving_average.begin
            import numpy as np

            from numba import guvectorize

            @guvectorize(['void(float64[:], intp[:], float64[:])'],
                         '(n),()->(n)')
            def move_mean(a, window_arr, out):
                window_width = window_arr[0]
                asum = 0.0
                count = 0
                for i in range(window_width):
                    asum += a[i]
                    count += 1
                    out[i] = asum / count
                for i in range(window_width, len(a)):
                    asum += a[i] - a[i - window_width]
                    out[i] = asum / count

            arr = np.arange(20, dtype=np.float64).reshape(2, 10)
            print(arr)
            print(move_mean(arr, 3))
            # magictoken.ex_moving_average.end

    def test_nogil(self):
        with captured_stdout():
            # magictoken.ex_no_gil.begin
            import math
            import threading
            from timeit import repeat

            import numpy as np
            from numba import jit

            nthreads = 4
            size = 10**6

            def func_np(a, b):
                """
                Control function using Numpy.
                """
                return np.exp(2.1 * a + 3.2 * b)

            @jit('void(double[:], double[:], double[:])', nopython=True,
                 nogil=True)
            def inner_func_nb(result, a, b):
                """
                Function under test.
                """
                for i in range(len(result)):
                    result[i] = math.exp(2.1 * a[i] + 3.2 * b[i])

            def timefunc(correct, s, func, *args, **kwargs):
                """
                Benchmark *func* and print out its runtime.
                """
                print(s.ljust(20), end=" ")
                # Make sure the function is compiled before the benchmark is
                # started
                res = func(*args, **kwargs)
                if correct is not None:
                    assert np.allclose(res, correct), (res, correct)
                # time it
                print('{:>5.0f} ms'.format(min(repeat(
                    lambda: func(*args, **kwargs), number=5, repeat=2)) * 1000))
                return res

            def make_singlethread(inner_func):
                """
                Run the given function inside a single thread.
                """
                def func(*args):
                    length = len(args[0])
                    result = np.empty(length, dtype=np.float64)
                    inner_func(result, *args)
                    return result
                return func

            def make_multithread(inner_func, numthreads):
                """
                Run the given function inside *numthreads* threads, splitting
                its arguments into equal-sized chunks.
                """
                def func_mt(*args):
                    length = len(args[0])
                    result = np.empty(length, dtype=np.float64)
                    args = (result,) + args
                    chunklen = (length + numthreads - 1) // numthreads
                    # Create argument tuples for each input chunk
                    chunks = [[arg[i * chunklen:(i + 1) * chunklen] for arg in
                               args] for i in range(numthreads)]
                    # Spawn one thread per chunk
                    threads = [threading.Thread(target=inner_func, args=chunk)
                               for chunk in chunks]
                    for thread in threads:
                        thread.start()
                    for thread in threads:
                        thread.join()
                    return result
                return func_mt

            func_nb = make_singlethread(inner_func_nb)
            func_nb_mt = make_multithread(inner_func_nb, nthreads)

            a = np.random.rand(size)
            b = np.random.rand(size)

            correct = timefunc(None, "numpy (1 thread)", func_np, a, b)
            timefunc(correct, "numba (1 thread)", func_nb, a, b)
            timefunc(correct, "numba (%d threads)" % nthreads, func_nb_mt, a, b)
            # magictoken.ex_no_gil.end

    def test_vectorize_one_signature(self):
        with captured_stdout():
            # magictoken.ex_vectorize_one_signature.begin
            from numba import vectorize, float64

            @vectorize([float64(float64, float64)])
            def f(x, y):
                return x + y
            # magictoken.ex_vectorize_one_signature.end

    def test_vectorize_multiple_signatures(self):
        with captured_stdout():
            # magictoken.ex_vectorize_multiple_signatures.begin
            from numba import vectorize, int32, int64, float32, float64
            import numpy as np

            @vectorize([int32(int32, int32),
                        int64(int64, int64),
                        float32(float32, float32),
                        float64(float64, float64)])
            def f(x, y):
                return x + y
            # magictoken.ex_vectorize_multiple_signatures.end

            # magictoken.ex_vectorize_return_call_one.begin
            a = np.arange(6)
            result = f(a, a)
            # result == array([ 0,  2,  4,  6,  8, 10])
            # magictoken.ex_vectorize_return_call_one.end

            self.assertIsInstance(result, np.ndarray)
            correct = np.array([0, 2, 4, 6, 8, 10])
            np.testing.assert_array_equal(result, correct)

            # magictoken.ex_vectorize_return_call_two.begin
            a = np.linspace(0, 1, 6)
            result = f(a, a)
            # Now, result == array([0. , 0.4, 0.8, 1.2, 1.6, 2. ])
            # magictoken.ex_vectorize_return_call_two.end

            self.assertIsInstance(result, np.ndarray)
            correct = np.array([0., 0.4, 0.8, 1.2, 1.6, 2. ])
            np.testing.assert_allclose(result, correct)

            # magictoken.ex_vectorize_return_call_three.begin
            a = np.arange(12).reshape(3, 4)
            # a == array([[ 0,  1,  2,  3],
            #             [ 4,  5,  6,  7],
            #             [ 8,  9, 10, 11]])

            result1 = f.reduce(a, axis=0)
            # result1 == array([12, 15, 18, 21])

            result2 = f.reduce(a, axis=1)
            # result2 == array([ 6, 22, 38])

            result3 = f.accumulate(a)
            # result3 == array([[ 0,  1,  2,  3],
            #                   [ 4,  6,  8, 10],
            #                   [12, 15, 18, 21]])

            result4 = f.accumulate(a, axis=1)
            # result3 == array([[ 0,  1,  3,  6],
            #                   [ 4,  9, 15, 22],
            #                   [ 8, 17, 27, 38]])
            # magictoken.ex_vectorize_return_call_three.end

            self.assertIsInstance(result1, np.ndarray)
            correct = np.array([12, 15, 18, 21])
            np.testing.assert_array_equal(result1, correct)

            self.assertIsInstance(result2, np.ndarray)
            correct = np.array([6, 22, 38])
            np.testing.assert_array_equal(result2, correct)

            self.assertIsInstance(result3, np.ndarray)
            correct = np.array([
                [0, 1, 2, 3],
                [4, 6, 8, 10],
                [12, 15, 18, 21]
            ])
            np.testing.assert_array_equal(result3, correct)

            self.assertIsInstance(result4, np.ndarray)
            correct = np.array([
                [0, 1, 3, 6],
                [4, 9, 15, 22],
                [8, 17, 27, 38]
            ])
            np.testing.assert_array_equal(result4, correct)

    def test_guvectorize(self):
        with captured_stdout():
            # magictoken.ex_guvectorize.begin
            from numba import guvectorize, int64
            import numpy as np

            @guvectorize([(int64[:], int64, int64[:])], '(n),()->(n)')
            def g(x, y, res):
                for i in range(x.shape[0]):
                    res[i] = x[i] + y
            # magictoken.ex_guvectorize.end

            # magictoken.ex_guvectorize_call_one.begin
            a = np.arange(5)
            result = g(a, 2)
            # result == array([2, 3, 4, 5, 6])
            # magictoken.ex_guvectorize_call_one.end

            self.assertIsInstance(result, np.ndarray)
            correct = np.array([2, 3, 4, 5, 6])
            np.testing.assert_array_equal(result, correct)

            # magictoken.ex_guvectorize_call_two.begin
            a = np.arange(6).reshape(2, 3)
            # a == array([[0, 1, 2],
            #             [3, 4, 5]])

            result1 = g(a, 10)
            # result1 == array([[10, 11, 12],
            #                   [13, 14, 15]])

            result2 = g(a, np.array([10, 20]))
            g(a, np.array([10, 20]))
            # result2 == array([[10, 11, 12],
            #                   [23, 24, 25]])
            # magictoken.ex_guvectorize_call_two.end

            self.assertIsInstance(result1, np.ndarray)
            correct = np.array([[10, 11, 12], [13, 14, 15]])
            np.testing.assert_array_equal(result1, correct)

            self.assertIsInstance(result2, np.ndarray)
            correct = np.array([[10, 11, 12], [23, 24, 25]])
            np.testing.assert_array_equal(result2, correct)

    def test_guvectorize_scalar_return(self):
        with captured_stdout():
            # magictoken.ex_guvectorize_scalar_return.begin
            from numba import guvectorize, int64
            import numpy as np

            @guvectorize([(int64[:], int64, int64[:])], '(n),()->()')
            def g(x, y, res):
                acc = 0
                for i in range(x.shape[0]):
                    acc += x[i] + y
                res[0] = acc
            # magictoken.ex_guvectorize_scalar_return.end

            # magictoken.ex_guvectorize_scalar_return_call.begin
            a = np.arange(5)
            result = g(a, 2)
            # At this point, result == 20.
            # magictoken.ex_guvectorize_scalar_return_call.end

            self.assertIsInstance(result, np.integer)
            self.assertEqual(result, 20)

    def test_guvectorize_jit(self):
        with captured_stdout():
            # magictoken.gufunc_jit.begin
            import numpy as np

            from numba import jit, guvectorize

            @guvectorize('(n)->(n)')
            def copy(x, res):
                for i in range(x.shape[0]):
                    res[i] = x[i]

            @jit(nopython=True)
            def jit_fn(x, res):
                copy(x, res)
            # magictoken.gufunc_jit.end

            # magictoken.gufunc_jit_call.begin
            x = np.arange(5, dtype='i4')
            res = np.zeros_like(x)
            jit_fn(x, res)
            # At this point, res == np.array([0, 1, 2, 3, 4], 'i4').
            # magictoken.gufunc_jit_call.end
            self.assertPreciseEqual(x, res)

    def test_guvectorize_jit_fail(self):
        with captured_stdout():
            # magictoken.gufunc_jit_fail.begin
            import numpy as np
            from numba import jit, guvectorize

            @guvectorize('(n)->(n)')
            def copy(x, res):
                for i in range(x.shape[0]):
                    res[i] = x[i]

            @jit(nopython=True)
            def jit_fn(x, res):
                copy(x, res)

            x = np.ones((1, 5))
            res = np.empty((5,))
            with self.assertRaises(ValueError) as raises:
                jit_fn(x, res)
            # magictoken.gufunc_jit_fail.end
            self.assertIn('Loop and array shapes are incompatible',
                          str(raises.exception))

    def test_guvectorize_overwrite(self):
        with captured_stdout():
            # magictoken.ex_guvectorize_overwrite.begin
            from numba import guvectorize, float64
            import numpy as np

            @guvectorize([(float64[:], float64[:])], '()->()')
            def init_values(invals, outvals):
                invals[0] = 6.5
                outvals[0] = 4.2
            # magictoken.ex_guvectorize_overwrite.end

            # magictoken.ex_guvectorize_overwrite_call_one.begin
            invals = np.zeros(shape=(3, 3), dtype=np.float64)
            # invals == array([[6.5, 6.5, 6.5],
            #                  [6.5, 6.5, 6.5],
            #                  [6.5, 6.5, 6.5]])

            outvals = init_values(invals)
            # outvals == array([[4.2, 4.2, 4.2],
            #                   [4.2, 4.2, 4.2],
            #                   [4.2, 4.2, 4.2]])
            # magictoken.ex_guvectorize_overwrite_call_one.end

            self.assertIsInstance(invals, np.ndarray)
            correct = np.array([
                [6.5, 6.5, 6.5],
                [6.5, 6.5, 6.5],
                [6.5, 6.5, 6.5]])
            np.testing.assert_array_equal(invals, correct)

            self.assertIsInstance(outvals, np.ndarray)
            correct = np.array([
                [4.2, 4.2, 4.2],
                [4.2, 4.2, 4.2],
                [4.2, 4.2, 4.2]])
            np.testing.assert_array_equal(outvals, correct)

            # magictoken.ex_guvectorize_overwrite_call_two.begin
            invals = np.zeros(shape=(3, 3), dtype=np.float32)
            # invals == array([[0., 0., 0.],
            #                  [0., 0., 0.],
            #                  [0., 0., 0.]], dtype=float32)
            outvals = init_values(invals)
            # outvals == array([[4.2, 4.2, 4.2],
            #                   [4.2, 4.2, 4.2],
            #                   [4.2, 4.2, 4.2]])
            print(invals)
            # invals == array([[0., 0., 0.],
            #                  [0., 0., 0.],
            #                  [0., 0., 0.]], dtype=float32)
            # magictoken.ex_guvectorize_overwrite_call_two.end

            self.assertIsInstance(invals, np.ndarray)
            correct = np.array([
                [0., 0., 0.],
                [0., 0., 0.],
                [0., 0., 0.]], dtype=np.float32)
            np.testing.assert_array_equal(invals, correct)

            self.assertIsInstance(outvals, np.ndarray)
            correct = np.array([
                [4.2, 4.2, 4.2],
                [4.2, 4.2, 4.2],
                [4.2, 4.2, 4.2]])
            np.testing.assert_array_equal(outvals, correct)

            # magictoken.ex_guvectorize_overwrite_call_three.begin
            @guvectorize(
                [(float64[:], float64[:])],
                '()->()',
                writable_args=('invals',)
            )
            def init_values(invals, outvals):
                invals[0] = 6.5
                outvals[0] = 4.2

            invals = np.zeros(shape=(3, 3), dtype=np.float32)
            # invals == array([[0., 0., 0.],
            #                  [0., 0., 0.],
            #                  [0., 0., 0.]], dtype=float32)
            outvals = init_values(invals)
            # outvals == array([[4.2, 4.2, 4.2],
            #                   [4.2, 4.2, 4.2],
            #                   [4.2, 4.2, 4.2]])
            print(invals)
            # invals == array([[6.5, 6.5, 6.5],
            #                  [6.5, 6.5, 6.5],
            #                  [6.5, 6.5, 6.5]], dtype=float32)
            # magictoken.ex_guvectorize_overwrite_call_three.end

            self.assertIsInstance(invals, np.ndarray)
            correct = np.array([
                [6.5, 6.5, 6.5],
                [6.5, 6.5, 6.5],
                [6.5, 6.5, 6.5]])
            np.testing.assert_array_equal(invals, correct)

            self.assertIsInstance(outvals, np.ndarray)
            correct = np.array([
                [4.2, 4.2, 4.2],
                [4.2, 4.2, 4.2],
                [4.2, 4.2, 4.2]])
            np.testing.assert_array_equal(outvals, correct)

    def test_vectorize_dynamic(self):
        with captured_stdout():
            # magictoken.ex_vectorize_dynamic.begin
            from numba import vectorize

            @vectorize
            def f(x, y):
                return x * y
            # magictoken.ex_vectorize_dynamic.end

            # magictoken.ex_vectorize_dynamic_call_one.begin
            result = f(3,4)
            # result == 12

            print(f.types)
            # ['ll->l']
            # magictoken.ex_vectorize_dynamic_call_one.end

            self.assertEqual(result, 12)
            if IS_WIN32:
                if numpy_version < (2, 0):
                    correct = ['ll->q']
                else:
                    correct = ['qq->q']
            else:
                correct = ['ll->l']
            self.assertEqual(f.types, correct)

            # magictoken.ex_vectorize_dynamic_call_two.begin
            result = f(1.,2.)
            # result == 2.0

            print(f.types)
            # ['ll->l', 'dd->d']
            # magictoken.ex_vectorize_dynamic_call_two.end

            self.assertEqual(result, 2.0)
            if IS_WIN32:
                if numpy_version < (2, 0):
                    correct = ['ll->q', 'dd->d']
                else:
                    correct = ['qq->q', 'dd->d']
            else:
                correct = ['ll->l', 'dd->d']
            self.assertEqual(f.types, correct)

            # magictoken.ex_vectorize_dynamic_call_three.begin
            result = f(1,2.)
            # result == 2.0

            print(f.types)
            # ['ll->l', 'dd->d']
            # magictoken.ex_vectorize_dynamic_call_three.end

            self.assertEqual(result, 2.0)
            if IS_WIN32:
                if numpy_version < (2, 0):
                    correct = ['ll->q', 'dd->d']
                else:
                    correct = ['qq->q', 'dd->d']
            else:
                correct = ['ll->l', 'dd->d']
            self.assertEqual(f.types, correct)

            # magictoken.ex_vectorize_dynamic_call_four.begin
            @vectorize
            def g(a, b):
                return a / b

            print(g(2.,3.))
            # 0.66666666666666663

            print(g(2,3))
            # 0.66666666666666663

            print(g.types)
            # ['dd->d']
            # magictoken.ex_vectorize_dynamic_call_four.end

            correct = ['dd->d']
            self.assertEqual(g.types, correct)

    def test_guvectorize_dynamic(self):
        with captured_stdout():
            # magictoken.ex_guvectorize_dynamic.begin
            from numba import guvectorize
            import numpy as np

            @guvectorize('(n),()->(n)')
            def g(x, y, res):
                for i in range(x.shape[0]):
                    res[i] = x[i] + y
            # magictoken.ex_guvectorize_dynamic.end

            # magictoken.ex_guvectorize_dynamic_call_one.begin
            x = np.arange(5, dtype=np.int64)
            y = 10
            res = np.zeros_like(x)
            g(x, y, res)
            # res == array([10, 11, 12, 13, 14])
            print(g.types)
            # ['ll->l']
            # magictoken.ex_guvectorize_dynamic_call_one.end

            correct = np.array([10, 11, 12, 13, 14])
            np.testing.assert_array_equal(res, correct)
            if IS_WIN32:
                correct = ['qq->q']
            else:
                correct = ['ll->l']
            self.assertEqual(g.types, correct)

            # magictoken.ex_guvectorize_dynamic_call_two.begin
            x = np.arange(5, dtype=np.double)
            y = 2.2
            res = np.zeros_like(x)
            g(x, y, res)
            # res == array([2.2, 3.2, 4.2, 5.2, 6.2])
            # magictoken.ex_guvectorize_dynamic_call_two.end

            # magictoken.ex_guvectorize_dynamic_call_three.begin
            print(g.types) # shorthand for g.ufunc.types
            # ['ll->l', 'dd->d']
            # magictoken.ex_guvectorize_dynamic_call_three.end

            if IS_WIN32:
                correct = ['qq->q', 'dd->d']
            else:
                correct = ['ll->l', 'dd->d']
            self.assertEqual(g.types, correct)

            # magictoken.ex_guvectorize_dynamic_call_four.begin
            x = np.arange(5, dtype=np.int64)
            y = 2
            res = np.zeros_like(x)
            g(x, y, res)
            print(res)
            # res == array([2, 3, 4, 5, 6])
            # magictoken.ex_guvectorize_dynamic_call_four.end

            correct = np.array([2, 3, 4, 5, 6])
            np.testing.assert_array_equal(res, correct)


if __name__ == '__main__':
    unittest.main()
