from numba import cuda, njit, types
from numba.core.errors import TypingError
from numba.core.extending import overload, overload_attribute
from numba.core.typing.typeof import typeof
from numba.cuda.testing import CUDATestCase, skip_on_cudasim, unittest
import numpy as np


# Dummy function definitions to overload

def generic_func_1():
    pass


def cuda_func_1():
    pass


def generic_func_2():
    pass


def cuda_func_2():
    pass


def generic_calls_generic():
    pass


def generic_calls_cuda():
    pass


def cuda_calls_generic():
    pass


def cuda_calls_cuda():
    pass


def target_overloaded():
    pass


def generic_calls_target_overloaded():
    pass


def cuda_calls_target_overloaded():
    pass


def target_overloaded_calls_target_overloaded():
    pass


# To recognise which functions are resolved for a call, we identify each with a
# prime number. Each function called multiplies a value by its prime (starting
# with the value 1), and we can check that the result is as expected based on
# the final value after all multiplications.

GENERIC_FUNCTION_1 = 2
CUDA_FUNCTION_1 = 3
GENERIC_FUNCTION_2 = 5
CUDA_FUNCTION_2 = 7
GENERIC_CALLS_GENERIC = 11
GENERIC_CALLS_CUDA = 13
CUDA_CALLS_GENERIC = 17
CUDA_CALLS_CUDA = 19
GENERIC_TARGET_OL = 23
CUDA_TARGET_OL = 29
GENERIC_CALLS_TARGET_OL = 31
CUDA_CALLS_TARGET_OL = 37
GENERIC_TARGET_OL_CALLS_TARGET_OL = 41
CUDA_TARGET_OL_CALLS_TARGET_OL = 43


# Overload implementations

@overload(generic_func_1, target='generic')
def ol_generic_func_1(x):
    def impl(x):
        x[0] *= GENERIC_FUNCTION_1
    return impl


@overload(cuda_func_1, target='cuda')
def ol_cuda_func_1(x):
    def impl(x):
        x[0] *= CUDA_FUNCTION_1
    return impl


@overload(generic_func_2, target='generic')
def ol_generic_func_2(x):
    def impl(x):
        x[0] *= GENERIC_FUNCTION_2
    return impl


@overload(cuda_func_2, target='cuda')
def ol_cuda_func(x):
    def impl(x):
        x[0] *= CUDA_FUNCTION_2
    return impl


@overload(generic_calls_generic, target='generic')
def ol_generic_calls_generic(x):
    def impl(x):
        x[0] *= GENERIC_CALLS_GENERIC
        generic_func_1(x)
    return impl


@overload(generic_calls_cuda, target='generic')
def ol_generic_calls_cuda(x):
    def impl(x):
        x[0] *= GENERIC_CALLS_CUDA
        cuda_func_1(x)
    return impl


@overload(cuda_calls_generic, target='cuda')
def ol_cuda_calls_generic(x):
    def impl(x):
        x[0] *= CUDA_CALLS_GENERIC
        generic_func_1(x)
    return impl


@overload(cuda_calls_cuda, target='cuda')
def ol_cuda_calls_cuda(x):
    def impl(x):
        x[0] *= CUDA_CALLS_CUDA
        cuda_func_1(x)
    return impl


@overload(target_overloaded, target='generic')
def ol_target_overloaded_generic(x):
    def impl(x):
        x[0] *= GENERIC_TARGET_OL
    return impl


@overload(target_overloaded, target='cuda')
def ol_target_overloaded_cuda(x):
    def impl(x):
        x[0] *= CUDA_TARGET_OL
    return impl


@overload(generic_calls_target_overloaded, target='generic')
def ol_generic_calls_target_overloaded(x):
    def impl(x):
        x[0] *= GENERIC_CALLS_TARGET_OL
        target_overloaded(x)
    return impl


@overload(cuda_calls_target_overloaded, target='cuda')
def ol_cuda_calls_target_overloaded(x):
    def impl(x):
        x[0] *= CUDA_CALLS_TARGET_OL
        target_overloaded(x)
    return impl


@overload(target_overloaded_calls_target_overloaded, target='generic')
def ol_generic_calls_target_overloaded_generic(x):
    def impl(x):
        x[0] *= GENERIC_TARGET_OL_CALLS_TARGET_OL
        target_overloaded(x)
    return impl


@overload(target_overloaded_calls_target_overloaded, target='cuda')
def ol_generic_calls_target_overloaded_cuda(x):
    def impl(x):
        x[0] *= CUDA_TARGET_OL_CALLS_TARGET_OL
        target_overloaded(x)
    return impl


@skip_on_cudasim('Overloading not supported in cudasim')
class TestOverload(CUDATestCase):
    def check_overload(self, kernel, expected):
        x = np.ones(1, dtype=np.int32)
        cuda.jit(kernel)[1, 1](x)
        self.assertEqual(x[0], expected)

    def check_overload_cpu(self, kernel, expected):
        x = np.ones(1, dtype=np.int32)
        njit(kernel)(x)
        self.assertEqual(x[0], expected)

    def test_generic(self):
        def kernel(x):
            generic_func_1(x)

        expected = GENERIC_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_cuda(self):
        def kernel(x):
            cuda_func_1(x)

        expected = CUDA_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_generic_and_cuda(self):
        def kernel(x):
            generic_func_1(x)
            cuda_func_1(x)

        expected = GENERIC_FUNCTION_1 * CUDA_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_call_two_generic_calls(self):
        def kernel(x):
            generic_func_1(x)
            generic_func_2(x)

        expected = GENERIC_FUNCTION_1 * GENERIC_FUNCTION_2
        self.check_overload(kernel, expected)

    def test_call_two_cuda_calls(self):
        def kernel(x):
            cuda_func_1(x)
            cuda_func_2(x)

        expected = CUDA_FUNCTION_1 * CUDA_FUNCTION_2
        self.check_overload(kernel, expected)

    def test_generic_calls_generic(self):
        def kernel(x):
            generic_calls_generic(x)

        expected = GENERIC_CALLS_GENERIC * GENERIC_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_generic_calls_cuda(self):
        def kernel(x):
            generic_calls_cuda(x)

        expected = GENERIC_CALLS_CUDA * CUDA_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_cuda_calls_generic(self):
        def kernel(x):
            cuda_calls_generic(x)

        expected = CUDA_CALLS_GENERIC * GENERIC_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_cuda_calls_cuda(self):
        def kernel(x):
            cuda_calls_cuda(x)

        expected = CUDA_CALLS_CUDA * CUDA_FUNCTION_1
        self.check_overload(kernel, expected)

    def test_call_target_overloaded(self):
        def kernel(x):
            target_overloaded(x)

        expected = CUDA_TARGET_OL
        self.check_overload(kernel, expected)

    def test_generic_calls_target_overloaded(self):
        def kernel(x):
            generic_calls_target_overloaded(x)

        expected = GENERIC_CALLS_TARGET_OL * CUDA_TARGET_OL
        self.check_overload(kernel, expected)

    def test_cuda_calls_target_overloaded(self):
        def kernel(x):
            cuda_calls_target_overloaded(x)

        expected = CUDA_CALLS_TARGET_OL * CUDA_TARGET_OL
        self.check_overload(kernel, expected)

    def test_target_overloaded_calls_target_overloaded(self):
        def kernel(x):
            target_overloaded_calls_target_overloaded(x)

        # Check the CUDA overloads are used on CUDA
        expected = CUDA_TARGET_OL_CALLS_TARGET_OL * CUDA_TARGET_OL
        self.check_overload(kernel, expected)

        # Also check that the CPU overloads are used on the CPU
        expected = GENERIC_TARGET_OL_CALLS_TARGET_OL * GENERIC_TARGET_OL
        self.check_overload_cpu(kernel, expected)

    def test_overload_attribute_target(self):
        MyDummy, MyDummyType = self.make_dummy_type()
        mydummy_type = typeof(MyDummy())

        @overload_attribute(MyDummyType, 'cuda_only', target='cuda')
        def ov_dummy_cuda_attr(obj):
            def imp(obj):
                return 42

            return imp

        # Ensure that we cannot use the CUDA target-specific attribute on the
        # CPU, and that an appropriate typing error is raised
        with self.assertRaisesRegex(TypingError,
                                    "Unknown attribute 'cuda_only'"):
            @njit(types.int64(mydummy_type))
            def illegal_target_attr_use(x):
                return x.cuda_only

        # Ensure that the CUDA target-specific attribute is usable and works
        # correctly when the target is CUDA - note eager compilation via
        # signature
        @cuda.jit(types.void(types.int64[::1], mydummy_type))
        def cuda_target_attr_use(res, dummy):
            res[0] = dummy.cuda_only


if __name__ == '__main__':
    unittest.main()
