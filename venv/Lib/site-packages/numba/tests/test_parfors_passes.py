"""
Tests for sub-components of parfors.
These tests are aimed to produce a good-enough coverage of parfor passes
so that refactoring on these passes are easier with faster testing turnaround.
"""
import unittest
from functools import reduce

import numpy as np

from numba import njit, typeof, prange, pndindex
import numba.parfors.parfor
from numba.core import (
    rewrites,
    typed_passes,
    untyped_passes,
    inline_closurecall,
    compiler,
    cpu,
    errors
)
from numba.core.registry import cpu_target
from numba.tests.support import (TestCase, is_parfors_unsupported)


class MyPipeline(object):
    def __init__(self, typingctx, targetctx, args, test_ir):
        self.state = compiler.StateDict()
        self.state.typingctx = typingctx
        self.state.targetctx = targetctx
        self.state.args = args
        self.state.func_ir = test_ir
        self.state.typemap = None
        self.state.return_type = None
        self.state.calltypes = None
        self.state.metadata = {}


class BaseTest(TestCase):
    @classmethod
    def _run_parfor(cls, test_func, args, swap_map=None):
        # TODO: refactor this with get_optimized_numba_ir() where this is
        #       copied from
        typingctx = cpu_target.typing_context
        targetctx = cpu_target.target_context
        test_ir = compiler.run_frontend(test_func)
        options = cpu.ParallelOptions(True)

        tp = MyPipeline(typingctx, targetctx, args, test_ir)

        typingctx.refresh()
        targetctx.refresh()

        inline_pass = inline_closurecall.InlineClosureCallPass(
            tp.state.func_ir, options, typed=True
        )
        inline_pass.run()

        rewrites.rewrite_registry.apply("before-inference", tp.state)

        untyped_passes.ReconstructSSA().run_pass(tp.state)

        (
            tp.state.typemap,
            tp.state.return_type,
            tp.state.calltypes,
            _
        ) = typed_passes.type_inference_stage(
            tp.state.typingctx, tp.state.targetctx, tp.state.func_ir,
            tp.state.args, None
        )

        typed_passes.PreLowerStripPhis().run_pass(tp.state)

        diagnostics = numba.parfors.parfor.ParforDiagnostics()

        preparfor_pass = numba.parfors.parfor.PreParforPass(
            tp.state.func_ir,
            tp.state.typemap,
            tp.state.calltypes,
            tp.state.typingctx,
            tp.state.targetctx,
            options,
            swapped=diagnostics.replaced_fns,
            replace_functions_map=swap_map,
        )
        preparfor_pass.run()

        rewrites.rewrite_registry.apply("after-inference", tp.state)
        return tp, options, diagnostics, preparfor_pass

    @classmethod
    def run_parfor_sub_pass(cls, test_func, args):
        tp, options, diagnostics, _ = cls._run_parfor(test_func, args)

        flags = compiler.Flags()
        parfor_pass = numba.parfors.parfor.ParforPass(
            tp.state.func_ir,
            tp.state.typemap,
            tp.state.calltypes,
            tp.state.return_type,
            tp.state.typingctx,
            tp.state.targetctx,
            options,
            flags,
            tp.state.metadata,
            diagnostics=diagnostics,
        )
        parfor_pass._pre_run()
        # Run subpass
        sub_pass = cls.sub_pass_class(parfor_pass)
        sub_pass.run(parfor_pass.func_ir.blocks)

        return sub_pass

    @classmethod
    def run_parfor_pre_pass(cls, test_func, args, swap_map=None):
        tp, options, diagnostics, preparfor_pass = cls._run_parfor(
            test_func, args, swap_map
        )
        return preparfor_pass

    def _run_parallel(self, func, *args, **kwargs):
        cfunc = njit(parallel=True)(func)
        expect = func(*args, **kwargs)
        got = cfunc(*args, **kwargs)
        return expect, got

    def run_parallel(self, func, *args, **kwargs):
        if is_parfors_unsupported:
            # Skip
            return
        expect, got = self._run_parallel(func, *args, **kwargs)
        self.assertPreciseEqual(expect, got)

    def run_parallel_check_output_array(self, func, *args, **kwargs):
        if is_parfors_unsupported:
            # Skip
            return
        expect, got = self._run_parallel(func, *args, **kwargs)
        # Don't match the value, just the return type. must return array
        self.assertIsInstance(expect, np.ndarray)
        self.assertIsInstance(got, np.ndarray)
        self.assertEqual(expect.shape, got.shape)

    def check_records(self, records):
        for rec in records:
            self.assertIsInstance(rec["new"], numba.parfors.parfor.Parfor)


class TestConvertSetItemPass(BaseTest):
    sub_pass_class = numba.parfors.parfor.ConvertSetItemPass

    def test_setitem_full_slice(self):
        def test_impl():
            n = 10
            a = np.ones(n)
            a[:] = 7
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "slice")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_setitem_slice_stop_bound(self):
        def test_impl():
            n = 10
            a = np.ones(n)
            a[:5] = 7
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "slice")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_setitem_slice_start_bound(self):
        def test_impl():
            n = 10
            a = np.ones(n)
            a[4:] = 7
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "slice")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_setitem_gather_if_scalar(self):
        def test_impl():
            n = 10
            a = np.ones(n)
            b = np.ones_like(a, dtype=np.bool_)
            a[b] = 7
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "masked_assign_broadcast_scalar")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_setitem_gather_if_array(self):
        def test_impl():
            n = 10
            a = np.ones(n)
            b = np.ones_like(a, dtype=np.bool_)
            c = np.ones_like(a)
            a[b] = c[b]
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "masked_assign_array")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)


class TestConvertNumpyPass(BaseTest):
    sub_pass_class = numba.parfors.parfor.ConvertNumpyPass

    def check_numpy_allocators(self, fn):
        def test_impl():
            n = 10
            a = fn(n)
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "numpy_allocator")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def check_numpy_random(self, fn):
        def test_impl():
            n = 10
            a = fn(n)
            return a

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "numpy_allocator")
        self.check_records(sub_pass.rewritten)

        self.run_parallel_check_output_array(test_impl)

    def test_numpy_allocators(self):
        fns = [np.ones, np.zeros]
        for fn in fns:
            with self.subTest(fn.__name__):
                self.check_numpy_allocators(fn)

    def test_numpy_random(self):
        fns = [np.random.random]
        for fn in fns:
            with self.subTest(fn.__name__):
                self.check_numpy_random(fn)

    def test_numpy_arrayexpr(self):
        def test_impl(a, b):
            return a + b

        a = b = np.ones(10)

        args = (a, b)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "arrayexpr")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl, *args)

    def test_numpy_arrayexpr_ufunc(self):
        def test_impl(a, b):
            return np.sin(-a) + np.float64(1) / np.sqrt(b)

        a = b = np.ones(10)

        args = (a, b)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "arrayexpr")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl, *args)

    def test_numpy_arrayexpr_boardcast(self):
        def test_impl(a, b):
            return a + b + np.array(1)

        a = np.ones(10)
        b = np.ones((3, 10))

        args = (a, b)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "arrayexpr")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl, *args)

    def test_numpy_arrayexpr_reshaped(self):
        def test_impl(a, b):
            a = a.reshape(1, a.size)  # shape[0] is now constant
            return a + b

        a = np.ones(10)
        b = np.ones(10)

        args = (a, b)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "arrayexpr")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl, *args)


class TestConvertReducePass(BaseTest):
    sub_pass_class = numba.parfors.parfor.ConvertReducePass

    def test_reduce_max_basic(self):
        def test_impl(arr):
            return reduce(lambda x, y: max(x, y), arr, 0.0)

        x = np.ones(10)
        args = (x,)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "reduce")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl, *args)

    def test_reduce_max_masked(self):
        def test_impl(arr):
            return reduce(lambda x, y: max(x, y), arr[arr > 5], 0.0)

        x = np.ones(10)
        args = (x,)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "reduce")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl, *args)


class TestConvertLoopPass(BaseTest):
    sub_pass_class = numba.parfors.parfor.ConvertLoopPass

    def test_prange_reduce_simple(self):
        def test_impl():
            n = 20
            c = 0
            for i in prange(n):
                c += i
            return c

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_prange_map_simple(self):
        def test_impl():
            n = 20
            arr = np.ones(n)
            for i in prange(n):
                arr[i] += i
            return arr

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_prange_two_args(self):
        def test_impl():
            n = 20
            arr = np.ones(n)
            for i in prange(3, n):
                arr[i] += i
            return arr

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_prange_three_args(self):
        def test_impl():
            n = 20
            arr = np.ones(n)
            for i in prange(3, n, 2):
                arr[i] += i
            return arr

        with self.assertRaises(errors.UnsupportedRewriteError) as raises:
            self.run_parfor_sub_pass(test_impl, ())
        self.assertIn(
            "Only constant step size of 1 is supported for prange",
            str(raises.exception),
        )

    def test_prange_map_inner_loop(self):
        def test_impl():
            n = 20
            arr = np.ones((n, n))
            for i in prange(n):
                for j in range(i):
                    arr[i, j] += i + j * n
            return arr

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)

        self.run_parallel(test_impl)

    def test_prange_map_nested_prange(self):
        def test_impl():
            n = 20
            arr = np.ones((n, n))
            for i in prange(n):
                for j in prange(i):
                    arr[i, j] += i + j * n
            return arr

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 2)
        self.check_records(sub_pass.rewritten)
        for record in sub_pass.rewritten:
            self.assertEqual(record["reason"], "loop")

        self.run_parallel(test_impl)

    def test_prange_map_none_index(self):
        def test_impl():
            n = 20
            arr = np.ones(n)
            for i in prange(n):
                inner = arr[i : i + 1]
                inner[()] += 1
            return arr

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        self.check_records(sub_pass.rewritten)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")

        self.run_parallel(test_impl)

    def test_prange_map_overwrite_index(self):
        def test_impl():
            n = 20
            arr = np.ones(n)
            for i in prange(n):
                i += 1
                arr[i - 1] = i
            return arr

        with self.assertRaises(errors.UnsupportedRewriteError) as raises:
            self.run_parfor_sub_pass(test_impl, ())
        self.assertIn(
            "Overwrite of parallel loop index",
            str(raises.exception),
        )

    def test_init_prange(self):
        def test_impl():
            n = 20
            arr = np.ones(n)
            numba.parfors.parfor.init_prange()
            val = 0
            for i in numba.parfors.parfor.internal_prange(len(arr)):
                val += arr[i]
            return val

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        self.check_records(sub_pass.rewritten)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")

        self.run_parallel(test_impl)

    def test_pndindex(self):
        def test_impl():
            n = 20
            arr = np.ones((n, n))
            val = 0
            for idx in pndindex(arr.shape):
                val += idx[0] * idx[1]
            return val

        sub_pass = self.run_parfor_sub_pass(test_impl, ())
        self.assertEqual(len(sub_pass.rewritten), 1)
        self.check_records(sub_pass.rewritten)

        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")

        self.run_parallel(test_impl)

    def test_numpy_sum(self):
        def test_impl(arr):
            return np.sum(arr)

        shape = 11, 13
        arr = np.arange(np.prod(shape)).reshape(shape)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)
        self.run_parallel(test_impl, *args)

    def test_numpy_sum_bool_array_masked(self):
        def test_impl(arr):
            sliced = arr[:, 0]
            return np.sum(arr[sliced >= 3, 1:2])

        shape = 11, 13
        arr = np.arange(np.prod(shape)).reshape(shape)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)
        self.run_parallel(test_impl, *args)

    def test_numpy_sum_int_array_masked(self):
        def test_impl(arr):
            sel = np.arange(arr.shape[1])
            return np.sum(arr[:, sel])

        shape = 11, 13
        arr = np.arange(np.prod(shape)).reshape(shape)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        # 1 for arange; 1 for sum
        self.assertEqual(len(sub_pass.rewritten), 2)
        for record in sub_pass.rewritten:
            self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)
        self.run_parallel(test_impl, *args)

    def test_numpy_fill_method(self):
        def test_impl(arr):
            arr.fill(3)
            return arr

        shape = 11, 13
        arr = np.arange(np.prod(shape)).reshape(shape)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        sub_pass = self.run_parfor_sub_pass(test_impl, argtypes)
        # 1 for arange; 1 for sum
        self.assertEqual(len(sub_pass.rewritten), 1)
        [record] = sub_pass.rewritten
        self.assertEqual(record["reason"], "loop")
        self.check_records(sub_pass.rewritten)
        self.run_parallel(test_impl, *args)


class TestPreParforPass(BaseTest):
    class sub_pass_class:
        def __init__(self, pass_states):
            pass

        def run(self, blocks):
            pass

    def test_dtype_conversion(self):
        # array.dtype are converted to np.dtype(array) in the PreParforPass
        def test_impl(a):
            b = np.ones(20, dtype=a.dtype)
            return b

        arr = np.arange(10)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        pre_pass = self.run_parfor_pre_pass(test_impl, argtypes)
        self.assertEqual(pre_pass.stats["replaced_func"], 0)
        self.assertEqual(pre_pass.stats["replaced_dtype"], 1)
        self.run_parallel(test_impl, *args)

    def test_sum_replacement(self):
        def test_impl(a):
            return np.sum(a)

        arr = np.arange(10)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        pre_pass = self.run_parfor_pre_pass(test_impl, argtypes)
        self.assertEqual(pre_pass.stats["replaced_func"], 1)
        self.assertEqual(pre_pass.stats["replaced_dtype"], 0)
        self.run_parallel(test_impl, *args)

    def test_replacement_map(self):
        def test_impl(a):
            return np.sum(a)

        arr = np.arange(10)
        args = (arr,)
        argtypes = [typeof(x) for x in args]

        swap_map = numba.parfors.parfor.swap_functions_map.copy()
        swap_map.pop(("sum", "numpy"))
        pre_pass = self.run_parfor_pre_pass(test_impl, argtypes, swap_map)
        self.assertEqual(pre_pass.stats["replaced_func"], 0)
        self.assertEqual(pre_pass.stats["replaced_dtype"], 0)
        self.run_parallel(test_impl, *args)


if __name__ == "__main__":
    unittest.main()
