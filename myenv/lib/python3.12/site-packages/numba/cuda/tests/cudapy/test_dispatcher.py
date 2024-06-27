import numpy as np
import threading

from numba import boolean, config, cuda, float32, float64, int32, int64, void
from numba.core.errors import TypingError
from numba.cuda.testing import skip_on_cudasim, unittest, CUDATestCase
import math


def add(x, y):
    return x + y


def add_kernel(r, x, y):
    r[0] = x + y


@skip_on_cudasim('Specialization not implemented in the simulator')
class TestDispatcherSpecialization(CUDATestCase):
    def _test_no_double_specialize(self, dispatcher, ty):

        with self.assertRaises(RuntimeError) as e:
            dispatcher.specialize(ty)

        self.assertIn('Dispatcher already specialized', str(e.exception))

    def test_no_double_specialize_sig_same_types(self):
        # Attempting to specialize a kernel jitted with a signature is illegal,
        # even for the same types the kernel is already specialized for.
        @cuda.jit('void(float32[::1])')
        def f(x):
            pass

        self._test_no_double_specialize(f, float32[::1])

    def test_no_double_specialize_no_sig_same_types(self):
        # Attempting to specialize an already-specialized kernel is illegal,
        # even for the same types the kernel is already specialized for.
        @cuda.jit
        def f(x):
            pass

        f_specialized = f.specialize(float32[::1])
        self._test_no_double_specialize(f_specialized, float32[::1])

    def test_no_double_specialize_sig_diff_types(self):
        # Attempting to specialize a kernel jitted with a signature is illegal.
        @cuda.jit('void(int32[::1])')
        def f(x):
            pass

        self._test_no_double_specialize(f, float32[::1])

    def test_no_double_specialize_no_sig_diff_types(self):
        # Attempting to specialize an already-specialized kernel is illegal.
        @cuda.jit
        def f(x):
            pass

        f_specialized = f.specialize(int32[::1])
        self._test_no_double_specialize(f_specialized, float32[::1])

    def test_specialize_cache_same(self):
        # Ensure that the same dispatcher is returned for the same argument
        # types, and that different dispatchers are returned for different
        # argument types.
        @cuda.jit
        def f(x):
            pass

        self.assertEqual(len(f.specializations), 0)

        f_float32 = f.specialize(float32[::1])
        self.assertEqual(len(f.specializations), 1)

        f_float32_2 = f.specialize(float32[::1])
        self.assertEqual(len(f.specializations), 1)
        self.assertIs(f_float32, f_float32_2)

        f_int32 = f.specialize(int32[::1])
        self.assertEqual(len(f.specializations), 2)
        self.assertIsNot(f_int32, f_float32)

    def test_specialize_cache_same_with_ordering(self):
        # Ensure that the same dispatcher is returned for the same argument
        # types, and that different dispatchers are returned for different
        # argument types, taking into account array ordering and multiple
        # arguments.
        @cuda.jit
        def f(x, y):
            pass

        self.assertEqual(len(f.specializations), 0)

        # 'A' order specialization
        f_f32a_f32a = f.specialize(float32[:], float32[:])
        self.assertEqual(len(f.specializations), 1)

        # 'C' order specialization
        f_f32c_f32c = f.specialize(float32[::1], float32[::1])
        self.assertEqual(len(f.specializations), 2)
        self.assertIsNot(f_f32a_f32a, f_f32c_f32c)

        # Reuse 'C' order specialization
        f_f32c_f32c_2 = f.specialize(float32[::1], float32[::1])
        self.assertEqual(len(f.specializations), 2)
        self.assertIs(f_f32c_f32c, f_f32c_f32c_2)


class TestDispatcher(CUDATestCase):
    """Most tests based on those in numba.tests.test_dispatcher."""

    def test_coerce_input_types(self):
        # Do not allow unsafe conversions if we can still compile other
        # specializations.
        c_add = cuda.jit(add_kernel)

        # Using a complex128 allows us to represent any result produced by the
        # test
        r = np.zeros(1, dtype=np.complex128)

        c_add[1, 1](r, 123, 456)
        self.assertEqual(r[0], add(123, 456))

        c_add[1, 1](r, 12.3, 45.6)
        self.assertEqual(r[0], add(12.3, 45.6))

        c_add[1, 1](r, 12.3, 45.6j)
        self.assertEqual(r[0], add(12.3, 45.6j))

        c_add[1, 1](r, 12300000000, 456)
        self.assertEqual(r[0], add(12300000000, 456))

        # Now force compilation of only a single specialization
        c_add = cuda.jit('(i4[::1], i4, i4)')(add_kernel)
        r = np.zeros(1, dtype=np.int32)

        c_add[1, 1](r, 123, 456)
        self.assertPreciseEqual(r[0], add(123, 456))

    @skip_on_cudasim('Simulator ignores signature')
    @unittest.expectedFailure
    def test_coerce_input_types_unsafe(self):
        # Implicit (unsafe) conversion of float to int, originally from
        # test_coerce_input_types. This test presently fails with the CUDA
        # Dispatcher because argument preparation is done by
        # _Kernel._prepare_args, which is currently inflexible with respect to
        # the types it can accept when preparing.
        #
        # This test is marked as xfail until future changes enable this
        # behavior.
        c_add = cuda.jit('(i4[::1], i4, i4)')(add_kernel)
        r = np.zeros(1, dtype=np.int32)

        c_add[1, 1](r, 12.3, 45.6)
        self.assertPreciseEqual(r[0], add(12, 45))

    @skip_on_cudasim('Simulator ignores signature')
    def test_coerce_input_types_unsafe_complex(self):
        # Implicit conversion of complex to int disallowed
        c_add = cuda.jit('(i4[::1], i4, i4)')(add_kernel)
        r = np.zeros(1, dtype=np.int32)

        with self.assertRaises(TypeError):
            c_add[1, 1](r, 12.3, 45.6j)

    @skip_on_cudasim('Simulator does not track overloads')
    def test_ambiguous_new_version(self):
        """Test compiling new version in an ambiguous case
        """
        c_add = cuda.jit(add_kernel)

        r = np.zeros(1, dtype=np.float64)
        INT = 1
        FLT = 1.5

        c_add[1, 1](r, INT, FLT)
        self.assertAlmostEqual(r[0], INT + FLT)
        self.assertEqual(len(c_add.overloads), 1)

        c_add[1, 1](r, FLT, INT)
        self.assertAlmostEqual(r[0], FLT + INT)
        self.assertEqual(len(c_add.overloads), 2)

        c_add[1, 1](r, FLT, FLT)
        self.assertAlmostEqual(r[0], FLT + FLT)
        self.assertEqual(len(c_add.overloads), 3)

        # The following call is ambiguous because (int, int) can resolve
        # to (float, int) or (int, float) with equal weight.
        c_add[1, 1](r, 1, 1)
        self.assertAlmostEqual(r[0], INT + INT)
        self.assertEqual(len(c_add.overloads), 4, "didn't compile a new "
                                                  "version")

    @skip_on_cudasim("Simulator doesn't support concurrent kernels")
    def test_lock(self):
        """
        Test that (lazy) compiling from several threads at once doesn't
        produce errors (see issue #908).
        """
        errors = []

        @cuda.jit
        def foo(r, x):
            r[0] = x + 1

        def wrapper():
            try:
                r = np.zeros(1, dtype=np.int64)
                foo[1, 1](r, 1)
                self.assertEqual(r[0], 2)
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=wrapper) for i in range(16)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertFalse(errors)

    def _test_explicit_signatures(self, sigs):
        f = cuda.jit(sigs)(add_kernel)

        # Exact signature matches
        r = np.zeros(1, dtype=np.int64)
        f[1, 1](r, 1, 2)
        self.assertPreciseEqual(r[0], 3)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, 1.5, 2.5)
        self.assertPreciseEqual(r[0], 4.0)

        if config.ENABLE_CUDASIM:
            # Pass - we can't check for no conversion on the simulator.
            return

        # No conversion
        with self.assertRaises(TypeError) as cm:
            r = np.zeros(1, dtype=np.complex128)
            f[1, 1](r, 1j, 1j)
        self.assertIn("No matching definition", str(cm.exception))
        self.assertEqual(len(f.overloads), 2, f.overloads)

    def test_explicit_signatures_strings(self):
        # Check with a list of strings for signatures
        sigs = ["(int64[::1], int64, int64)",
                "(float64[::1], float64, float64)"]
        self._test_explicit_signatures(sigs)

    def test_explicit_signatures_tuples(self):
        # Check with a list of tuples of argument types for signatures
        sigs = [(int64[::1], int64, int64), (float64[::1], float64, float64)]
        self._test_explicit_signatures(sigs)

    def test_explicit_signatures_signatures(self):
        # Check with a list of Signature objects for signatures
        sigs = [void(int64[::1], int64, int64),
                void(float64[::1], float64, float64)]
        self._test_explicit_signatures(sigs)

    def test_explicit_signatures_mixed(self):
        # Check when we mix types of signature objects in a list of signatures

        # Tuple and string
        sigs = [(int64[::1], int64, int64),
                "(float64[::1], float64, float64)"]
        self._test_explicit_signatures(sigs)

        # Tuple and Signature object
        sigs = [(int64[::1], int64, int64),
                void(float64[::1], float64, float64)]
        self._test_explicit_signatures(sigs)

        # Signature object and string
        sigs = [void(int64[::1], int64, int64),
                "(float64[::1], float64, float64)"]
        self._test_explicit_signatures(sigs)

    def test_explicit_signatures_same_type_class(self):
        # A more interesting one...
        # (Note that the type of r is deliberately float64 in both cases so
        # that dispatch is differentiated on the types of x and y only, to
        # closely preserve the intent of the original test from
        # numba.tests.test_dispatcher)
        sigs = ["(float64[::1], float32, float32)",
                "(float64[::1], float64, float64)"]
        f = cuda.jit(sigs)(add_kernel)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, np.float32(1), np.float32(2**-25))
        self.assertPreciseEqual(r[0], 1.0)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, 1, 2**-25)
        self.assertPreciseEqual(r[0], 1.0000000298023224)

    @skip_on_cudasim('No overload resolution in the simulator')
    def test_explicit_signatures_ambiguous_resolution(self):
        # Fail to resolve ambiguity between the two best overloads
        # (Also deliberate float64[::1] for the first argument in all cases)
        f = cuda.jit(["(float64[::1], float32, float64)",
                      "(float64[::1], float64, float32)",
                      "(float64[::1], int64, int64)"])(add_kernel)
        with self.assertRaises(TypeError) as cm:
            r = np.zeros(1, dtype=np.float64)
            f[1, 1](r, 1.0, 2.0)

        # The two best matches are output in the error message, as well
        # as the actual argument types.
        self.assertRegex(
            str(cm.exception),
            r"Ambiguous overloading for <function add_kernel [^>]*> "
            r"\(Array\(float64, 1, 'C', False, aligned=True\), float64,"
            r" float64\):\n"
            r"\(Array\(float64, 1, 'C', False, aligned=True\), float32,"
            r" float64\) -> none\n"
            r"\(Array\(float64, 1, 'C', False, aligned=True\), float64,"
            r" float32\) -> none"
        )
        # The integer signature is not part of the best matches
        self.assertNotIn("int64", str(cm.exception))

    @skip_on_cudasim('Simulator does not use _prepare_args')
    @unittest.expectedFailure
    def test_explicit_signatures_unsafe(self):
        # These tests are from test_explicit_signatures, but have to be xfail
        # at present because _prepare_args in the CUDA target cannot handle
        # unsafe conversions of arguments.
        f = cuda.jit("(int64[::1], int64, int64)")(add_kernel)
        r = np.zeros(1, dtype=np.int64)

        # Approximate match (unsafe conversion)
        f[1, 1](r, 1.5, 2.5)
        self.assertPreciseEqual(r[0], 3)
        self.assertEqual(len(f.overloads), 1, f.overloads)

        sigs = ["(int64[::1], int64, int64)",
                "(float64[::1], float64, float64)"]
        f = cuda.jit(sigs)(add_kernel)
        r = np.zeros(1, dtype=np.float64)
        # Approximate match (int32 -> float64 is a safe conversion)
        f[1, 1](r, np.int32(1), 2.5)
        self.assertPreciseEqual(r[0], 3.5)

    def add_device_usecase(self, sigs):
        # Generate a kernel that calls the add device function compiled with a
        # given set of signatures
        add_device = cuda.jit(sigs, device=True)(add)

        @cuda.jit
        def f(r, x, y):
            r[0] = add_device(x, y)

        return f

    def test_explicit_signatures_device(self):
        # Tests similar to test_explicit_signatures, but on a device function
        # instead of a kernel
        sigs = ["(int64, int64)", "(float64, float64)"]
        f = self.add_device_usecase(sigs)

        # Exact signature matches
        r = np.zeros(1, dtype=np.int64)
        f[1, 1](r, 1, 2)
        self.assertPreciseEqual(r[0], 3)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, 1.5, 2.5)
        self.assertPreciseEqual(r[0], 4.0)

        if config.ENABLE_CUDASIM:
            # Pass - we can't check for no conversion on the simulator.
            return

        # No conversion
        with self.assertRaises(TypingError) as cm:
            r = np.zeros(1, dtype=np.complex128)
            f[1, 1](r, 1j, 1j)

        msg = str(cm.exception)
        self.assertIn("Invalid use of type", msg)
        self.assertIn("with parameters (complex128, complex128)", msg)
        self.assertEqual(len(f.overloads), 2, f.overloads)

    def test_explicit_signatures_device_same_type_class(self):
        # A more interesting one...
        # (Note that the type of r is deliberately float64 in both cases so
        # that dispatch is differentiated on the types of x and y only, to
        # closely preserve the intent of the original test from
        # numba.tests.test_dispatcher)
        sigs = ["(float32, float32)", "(float64, float64)"]
        f = self.add_device_usecase(sigs)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, np.float32(1), np.float32(2**-25))
        self.assertPreciseEqual(r[0], 1.0)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, 1, 2**-25)
        self.assertPreciseEqual(r[0], 1.0000000298023224)

    def test_explicit_signatures_device_ambiguous(self):
        # Ambiguity between the two best overloads resolves. This is somewhat
        # surprising given that ambiguity is not permitted for dispatching
        # overloads when launching a kernel, but seems to be the general
        # behaviour of Numba (See Issue #8307:
        # https://github.com/numba/numba/issues/8307).
        sigs = ["(float32, float64)", "(float64, float32)", "(int64, int64)"]
        f = self.add_device_usecase(sigs)

        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, 1.5, 2.5)
        self.assertPreciseEqual(r[0], 4.0)

    @skip_on_cudasim('CUDA Simulator does not force casting')
    def test_explicit_signatures_device_unsafe(self):
        # These tests are from test_explicit_signatures. The device function
        # variant of these tests can succeed on CUDA because the compilation
        # can handle unsafe casting (c.f. test_explicit_signatures_unsafe which
        # has to xfail due to _prepare_args not supporting unsafe casting).
        sigs = ["(int64, int64)"]
        f = self.add_device_usecase(sigs)

        # Approximate match (unsafe conversion)
        r = np.zeros(1, dtype=np.int64)
        f[1, 1](r, 1.5, 2.5)
        self.assertPreciseEqual(r[0], 3)
        self.assertEqual(len(f.overloads), 1, f.overloads)

        sigs = ["(int64, int64)", "(float64, float64)"]
        f = self.add_device_usecase(sigs)

        # Approximate match (int32 -> float64 is a safe conversion)
        r = np.zeros(1, dtype=np.float64)
        f[1, 1](r, np.int32(1), 2.5)
        self.assertPreciseEqual(r[0], 3.5)

    def test_dispatcher_docstring(self):
        # Ensure that CUDA-jitting a function preserves its docstring. See
        # Issue #5902: https://github.com/numba/numba/issues/5902

        @cuda.jit
        def add_kernel(a, b):
            """Add two integers, kernel version"""

        @cuda.jit(device=True)
        def add_device(a, b):
            """Add two integers, device version"""

        self.assertEqual("Add two integers, kernel version", add_kernel.__doc__)
        self.assertEqual("Add two integers, device version", add_device.__doc__)


@skip_on_cudasim("CUDA simulator doesn't implement kernel properties")
class TestDispatcherKernelProperties(CUDATestCase):
    def test_get_regs_per_thread_unspecialized(self):
        # A kernel where the register usage per thread is likely to differ
        # between different specializations
        @cuda.jit
        def pi_sin_array(x, n):
            i = cuda.grid(1)
            if i < n:
                x[i] = 3.14 * math.sin(x[i])

        # Call the kernel with different arguments to create two different
        # definitions within the Dispatcher object
        N = 10
        arr_f32 = np.zeros(N, dtype=np.float32)
        arr_f64 = np.zeros(N, dtype=np.float64)

        pi_sin_array[1, N](arr_f32, N)
        pi_sin_array[1, N](arr_f64, N)

        # Check we get a positive integer for the two different variations
        sig_f32 = void(float32[::1], int64)
        sig_f64 = void(float64[::1], int64)
        regs_per_thread_f32 = pi_sin_array.get_regs_per_thread(sig_f32)
        regs_per_thread_f64 = pi_sin_array.get_regs_per_thread(sig_f64)

        self.assertIsInstance(regs_per_thread_f32, int)
        self.assertIsInstance(regs_per_thread_f64, int)

        self.assertGreater(regs_per_thread_f32, 0)
        self.assertGreater(regs_per_thread_f64, 0)

        # Check that getting the registers per thread for all signatures
        # provides the same values as getting the registers per thread for
        # individual signatures.
        regs_per_thread_all = pi_sin_array.get_regs_per_thread()
        self.assertEqual(regs_per_thread_all[sig_f32.args],
                         regs_per_thread_f32)
        self.assertEqual(regs_per_thread_all[sig_f64.args],
                         regs_per_thread_f64)

        if regs_per_thread_f32 == regs_per_thread_f64:
            # If the register usage is the same for both variants, there may be
            # a bug, but this may also be an artifact of the compiler / driver
            # / device combination, so produce an informational message only.
            print('f32 and f64 variant thread usages are equal.')
            print('This may warrant some investigation. Devices:')
            cuda.detect()

    def test_get_regs_per_thread_specialized(self):
        @cuda.jit(void(float32[::1], int64))
        def pi_sin_array(x, n):
            i = cuda.grid(1)
            if i < n:
                x[i] = 3.14 * math.sin(x[i])

        # Check we get a positive integer for the specialized variation
        regs_per_thread = pi_sin_array.get_regs_per_thread()
        self.assertIsInstance(regs_per_thread, int)
        self.assertGreater(regs_per_thread, 0)

    def test_get_const_mem_unspecialized(self):
        @cuda.jit
        def const_fmt_string(val, to_print):
            # We guard the print with a conditional to prevent noise from the
            # test suite
            if to_print:
                print(val)

        # Call the kernel with different arguments to create two different
        # definitions within the Dispatcher object
        const_fmt_string[1, 1](1, False)
        const_fmt_string[1, 1](1.0, False)

        # Check we get a positive integer for the two different variations
        sig_i64 = void(int64, boolean)
        sig_f64 = void(float64, boolean)
        const_mem_size_i64 = const_fmt_string.get_const_mem_size(sig_i64)
        const_mem_size_f64 = const_fmt_string.get_const_mem_size(sig_f64)

        self.assertIsInstance(const_mem_size_i64, int)
        self.assertIsInstance(const_mem_size_f64, int)

        # 6 bytes for the equivalent of b'%lld\n\0'
        self.assertGreaterEqual(const_mem_size_i64, 6)
        # 4 bytes for the equivalent of b'%f\n\0'
        self.assertGreaterEqual(const_mem_size_f64, 4)

        # Check that getting the const memory size for all signatures
        # provides the same values as getting the const memory size for
        # individual signatures.

        const_mem_size_all = const_fmt_string.get_const_mem_size()
        self.assertEqual(const_mem_size_all[sig_i64.args], const_mem_size_i64)
        self.assertEqual(const_mem_size_all[sig_f64.args], const_mem_size_f64)

    def test_get_const_mem_specialized(self):
        arr = np.arange(32, dtype=np.int64)
        sig = void(int64[::1])

        @cuda.jit(sig)
        def const_array_use(x):
            C = cuda.const.array_like(arr)
            i = cuda.grid(1)
            x[i] = C[i]

        const_mem_size = const_array_use.get_const_mem_size(sig)
        self.assertIsInstance(const_mem_size, int)
        self.assertGreaterEqual(const_mem_size, arr.nbytes)

    def test_get_shared_mem_per_block_unspecialized(self):
        N = 10

        # A kernel where the shared memory per block is likely to differ
        # between different specializations
        @cuda.jit
        def simple_smem(ary):
            sm = cuda.shared.array(N, dtype=ary.dtype)
            for j in range(N):
                sm[j] = j
            for j in range(N):
                ary[j] = sm[j]

        # Call the kernel with different arguments to create two different
        # definitions within the Dispatcher object
        arr_f32 = np.zeros(N, dtype=np.float32)
        arr_f64 = np.zeros(N, dtype=np.float64)

        simple_smem[1, 1](arr_f32)
        simple_smem[1, 1](arr_f64)

        sig_f32 = void(float32[::1])
        sig_f64 = void(float64[::1])

        sh_mem_f32 = simple_smem.get_shared_mem_per_block(sig_f32)
        sh_mem_f64 = simple_smem.get_shared_mem_per_block(sig_f64)

        self.assertIsInstance(sh_mem_f32, int)
        self.assertIsInstance(sh_mem_f64, int)

        self.assertEqual(sh_mem_f32, N * 4)
        self.assertEqual(sh_mem_f64, N * 8)

        # Check that getting the shared memory per block for all signatures
        # provides the same values as getting the shared mem per block for
        # individual signatures.
        sh_mem_f32_all = simple_smem.get_shared_mem_per_block()
        sh_mem_f64_all = simple_smem.get_shared_mem_per_block()
        self.assertEqual(sh_mem_f32_all[sig_f32.args], sh_mem_f32)
        self.assertEqual(sh_mem_f64_all[sig_f64.args], sh_mem_f64)

    def test_get_shared_mem_per_block_specialized(self):
        @cuda.jit(void(float32[::1]))
        def simple_smem(ary):
            sm = cuda.shared.array(100, dtype=float32)
            i = cuda.grid(1)
            if i == 0:
                for j in range(100):
                    sm[j] = j
            cuda.syncthreads()
            ary[i] = sm[i]

        shared_mem_per_block = simple_smem.get_shared_mem_per_block()
        self.assertIsInstance(shared_mem_per_block, int)
        self.assertEqual(shared_mem_per_block, 400)

    def test_get_max_threads_per_block_unspecialized(self):
        N = 10

        @cuda.jit
        def simple_maxthreads(ary):
            i = cuda.grid(1)
            ary[i] = i

        arr_f32 = np.zeros(N, dtype=np.float32)
        simple_maxthreads[1, 1](arr_f32)
        sig_f32 = void(float32[::1])
        max_threads_f32 = simple_maxthreads.get_max_threads_per_block(sig_f32)

        self.assertIsInstance(max_threads_f32, int)
        self.assertGreater(max_threads_f32, 0)

        max_threads_f32_all = simple_maxthreads.get_max_threads_per_block()
        self.assertEqual(max_threads_f32_all[sig_f32.args], max_threads_f32)

    def test_get_local_mem_per_thread_unspecialized(self):
        # NOTE: A large amount of local memory must be allocated
        # otherwise the compiler will optimize out the call to
        # cuda.local.array and use local registers instead
        N = 1000

        @cuda.jit
        def simple_lmem(ary):
            lm = cuda.local.array(N, dtype=ary.dtype)
            for j in range(N):
                lm[j] = j
            for j in range(N):
                ary[j] = lm[j]

        # Call the kernel with different arguments to create two different
        # definitions within the Dispatcher object
        arr_f32 = np.zeros(N, dtype=np.float32)
        arr_f64 = np.zeros(N, dtype=np.float64)

        simple_lmem[1, 1](arr_f32)
        simple_lmem[1, 1](arr_f64)

        sig_f32 = void(float32[::1])
        sig_f64 = void(float64[::1])
        local_mem_f32 = simple_lmem.get_local_mem_per_thread(sig_f32)
        local_mem_f64 = simple_lmem.get_local_mem_per_thread(sig_f64)
        self.assertIsInstance(local_mem_f32, int)
        self.assertIsInstance(local_mem_f64, int)

        self.assertGreaterEqual(local_mem_f32, N * 4)
        self.assertGreaterEqual(local_mem_f64, N * 8)

        # Check that getting the local memory per thread for all signatures
        # provides the same values as getting the shared mem per block for
        # individual signatures.
        local_mem_all = simple_lmem.get_local_mem_per_thread()
        self.assertEqual(local_mem_all[sig_f32.args], local_mem_f32)
        self.assertEqual(local_mem_all[sig_f64.args], local_mem_f64)

    def test_get_local_mem_per_thread_specialized(self):
        # NOTE: A large amount of local memory must be allocated
        # otherwise the compiler will optimize out the call to
        # cuda.local.array and use local registers instead
        N = 1000

        @cuda.jit(void(float32[::1]))
        def simple_lmem(ary):
            lm = cuda.local.array(N, dtype=ary.dtype)
            for j in range(N):
                lm[j] = j
            for j in range(N):
                ary[j] = lm[j]

        local_mem_per_thread = simple_lmem.get_local_mem_per_thread()
        self.assertIsInstance(local_mem_per_thread, int)
        self.assertGreaterEqual(local_mem_per_thread, N * 4)


if __name__ == '__main__':
    unittest.main()
