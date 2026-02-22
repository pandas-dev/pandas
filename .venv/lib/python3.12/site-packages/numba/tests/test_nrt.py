import math
import os
import platform
import sys
import re

import numpy as np

from numba import njit
from numba.core import types
from numba.core.runtime import (
    rtsys,
    nrtopt,
    _nrt_python,
    nrt,
)
from numba.core.extending import intrinsic, include_path
from numba.core.typing import signature
from numba.core.imputils import impl_ret_untracked
from llvmlite import ir
import llvmlite.binding as llvm
from numba.core.unsafe.nrt import NRT_get_api

from numba.tests.support import (EnableNRTStatsMixin, TestCase, temp_directory,
                                 import_dynamic, skip_if_32bit,
                                 skip_unless_cffi, run_in_subprocess)
from numba.core.registry import cpu_target
import unittest


linux_only = unittest.skipIf(not sys.platform.startswith('linux'),
                             'linux only test')
x86_only = unittest.skipIf(platform.machine() not in ('i386', 'x86_64'),
                           'x86 only test')


class Dummy(object):
    alive = 0

    def __init__(self):
        type(self).alive += 1

    def __del__(self):
        type(self).alive -= 1


class TestNrtMemInfoNotInitialized(unittest.TestCase):
    """
    Unit test for checking the use of the NRT fails if the
    initialization sequence has not been run.
    """
    _numba_parallel_test_ = False

    def test_init_fail(self):
        methods = {'library': (),
                   'meminfo_new': ((), ()),
                   'meminfo_alloc': ((),),
                   }

        for meth, args in methods.items():
            try:
                with self.assertRaises(RuntimeError) as raises:
                    rtsys._init = False
                    fn = getattr(rtsys, meth)
                    fn(*args)

                msg = "Runtime must be initialized before use."
                self.assertIn(msg, str(raises.exception))
            finally:
                rtsys._init = True


class TestNrtMemInfo(unittest.TestCase):
    """
    Unit test for core MemInfo functionality
    """

    def setUp(self):
        # Reset the Dummy class
        Dummy.alive = 0
        # initialize the NRT (in case the tests are run in isolation)
        rtsys.initialize(cpu_target.target_context)
        super(TestNrtMemInfo, self).setUp()

    def test_meminfo_refct_1(self):
        d = Dummy()
        self.assertEqual(Dummy.alive, 1)
        addr = 0xdeadcafe  # some made up location

        mi = rtsys.meminfo_new(addr, d)
        self.assertEqual(mi.refcount, 1)
        del d
        self.assertEqual(Dummy.alive, 1)
        mi.acquire()
        self.assertEqual(mi.refcount, 2)
        self.assertEqual(Dummy.alive, 1)
        mi.release()
        self.assertEqual(mi.refcount, 1)
        del mi
        self.assertEqual(Dummy.alive, 0)

    def test_meminfo_refct_2(self):
        d = Dummy()
        self.assertEqual(Dummy.alive, 1)
        addr = 0xdeadcafe  # some made up location

        mi = rtsys.meminfo_new(addr, d)
        self.assertEqual(mi.refcount, 1)
        del d
        self.assertEqual(Dummy.alive, 1)
        for ct in range(100):
            mi.acquire()
        self.assertEqual(mi.refcount, 1 + 100)
        self.assertEqual(Dummy.alive, 1)
        for _ in range(100):
            mi.release()
        self.assertEqual(mi.refcount, 1)
        del mi
        self.assertEqual(Dummy.alive, 0)

    def test_fake_memoryview(self):
        d = Dummy()
        self.assertEqual(Dummy.alive, 1)
        addr = 0xdeadcafe  # some made up location

        mi = rtsys.meminfo_new(addr, d)
        self.assertEqual(mi.refcount, 1)
        mview = memoryview(mi)
        self.assertEqual(mi.refcount, 1)
        self.assertEqual(addr, mi.data)
        self.assertFalse(mview.readonly)
        self.assertIs(mi, mview.obj)
        self.assertTrue(mview.c_contiguous)
        self.assertEqual(mview.itemsize, 1)
        self.assertEqual(mview.ndim, 1)
        del d
        del mi

        self.assertEqual(Dummy.alive, 1)
        del mview
        self.assertEqual(Dummy.alive, 0)

    def test_memoryview(self):
        from ctypes import c_uint32, c_void_p, POINTER, cast

        dtype = np.dtype(np.uint32)
        bytesize = dtype.itemsize * 10
        mi = rtsys.meminfo_alloc(bytesize, safe=True)
        addr = mi.data
        c_arr = cast(c_void_p(mi.data), POINTER(c_uint32 * 10))
        # Check 0xCB-filling
        for i in range(10):
            self.assertEqual(c_arr.contents[i], 0xcbcbcbcb)

        # Init array with ctypes
        for i in range(10):
            c_arr.contents[i] = i + 1
        mview = memoryview(mi)
        self.assertEqual(mview.nbytes, bytesize)
        self.assertFalse(mview.readonly)
        self.assertIs(mi, mview.obj)
        self.assertTrue(mview.c_contiguous)
        self.assertEqual(mview.itemsize, 1)
        self.assertEqual(mview.ndim, 1)
        del mi
        arr = np.ndarray(dtype=dtype, shape=mview.nbytes // dtype.itemsize,
                         buffer=mview)
        del mview
        # Modify array with NumPy
        np.testing.assert_equal(np.arange(arr.size) + 1, arr)

        arr += 1

        # Check value reflected in ctypes
        for i in range(10):
            self.assertEqual(c_arr.contents[i], i + 2)

        self.assertEqual(arr.ctypes.data, addr)
        del arr
        # At this point the memory is zero filled
        # We can't check this deterministically because the memory could be
        # consumed by another thread.

    def test_buffer(self):
        from ctypes import c_uint32, c_void_p, POINTER, cast

        dtype = np.dtype(np.uint32)
        bytesize = dtype.itemsize * 10
        mi = rtsys.meminfo_alloc(bytesize, safe=True)
        self.assertEqual(mi.refcount, 1)
        addr = mi.data
        c_arr = cast(c_void_p(addr), POINTER(c_uint32 * 10))
        # Check 0xCB-filling
        for i in range(10):
            self.assertEqual(c_arr.contents[i], 0xcbcbcbcb)

        # Init array with ctypes
        for i in range(10):
            c_arr.contents[i] = i + 1

        arr = np.ndarray(dtype=dtype, shape=bytesize // dtype.itemsize,
                         buffer=mi)
        self.assertEqual(mi.refcount, 1)
        del mi
        # Modify array with NumPy
        np.testing.assert_equal(np.arange(arr.size) + 1, arr)

        arr += 1

        # Check value reflected in ctypes
        for i in range(10):
            self.assertEqual(c_arr.contents[i], i + 2)

        self.assertEqual(arr.ctypes.data, addr)
        del arr
        # At this point the memory is zero filled
        # We can't check this deterministically because the memory could be
        # consumed by another thread.

    @skip_if_32bit
    def test_allocate_invalid_size(self):
        # Checks that attempting to allocate too big a region fails gracefully.
        size = types.size_t.maxval // 8 // 2
        for pred in (True, False):
            with self.assertRaises(MemoryError) as raises:
                rtsys.meminfo_alloc(size, safe=pred)
            self.assertIn(f"Requested allocation of {size} bytes failed.",
                          str(raises.exception))

    def test_allocate_negative_size(self):
        # Checks that attempting to allocate negative number of bytes fails
        # gracefully.
        size = -10
        for pred in (True, False):
            with self.assertRaises(ValueError) as raises:
                rtsys.meminfo_alloc(size, safe=pred)
            msg = f"Cannot allocate a negative number of bytes: {size}."
            self.assertIn(msg, str(raises.exception))


class TestTracemalloc(unittest.TestCase):
    """
    Test NRT-allocated memory can be tracked by tracemalloc.
    """

    def measure_memory_diff(self, func):
        try:
            import tracemalloc
        except ImportError:
            self.skipTest("tracemalloc not available")
        tracemalloc.start()
        try:
            before = tracemalloc.take_snapshot()
            # Keep the result and only delete it after taking a snapshot
            res = func()
            after = tracemalloc.take_snapshot()
            del res
            return after.compare_to(before, 'lineno')
        finally:
            tracemalloc.stop()

    def test_snapshot(self):
        N = 1000000
        dtype = np.int8

        @njit
        def alloc_nrt_memory():
            """
            Allocate and return a large array.
            """
            return np.empty(N, dtype)

        def keep_memory():
            return alloc_nrt_memory()

        def release_memory():
            alloc_nrt_memory()

        alloc_lineno = keep_memory.__code__.co_firstlineno + 1

        # Warmup JIT
        alloc_nrt_memory()

        # The large NRT-allocated array should appear topmost in the diff
        diff = self.measure_memory_diff(keep_memory)
        stat = diff[0]
        # There is a slight overhead, so the allocated size won't exactly be N
        self.assertGreaterEqual(stat.size, N)
        self.assertLess(stat.size, N * 1.015,
                        msg=("Unexpected allocation overhead encountered. "
                             "May be due to difference in CPython "
                             "builds or running under coverage"))
        frame = stat.traceback[0]
        self.assertEqual(os.path.basename(frame.filename), "test_nrt.py")
        self.assertEqual(frame.lineno, alloc_lineno)

        # If NRT memory is released before taking a snapshot, it shouldn't
        # appear.
        diff = self.measure_memory_diff(release_memory)
        stat = diff[0]
        # Something else appears, but nothing the magnitude of N
        self.assertLess(stat.size, N * 0.01)


class TestNRTIssue(TestCase):
    def test_issue_with_refct_op_pruning(self):
        """
        GitHub Issue #1244 https://github.com/numba/numba/issues/1244
        """
        @njit
        def calculate_2D_vector_mag(vector):
            x, y = vector

            return math.sqrt(x ** 2 + y ** 2)

        @njit
        def normalize_2D_vector(vector):
            normalized_vector = np.empty(2, dtype=np.float64)

            mag = calculate_2D_vector_mag(vector)
            x, y = vector

            normalized_vector[0] = x / mag
            normalized_vector[1] = y / mag

            return normalized_vector

        @njit
        def normalize_vectors(num_vectors, vectors):
            normalized_vectors = np.empty((num_vectors, 2), dtype=np.float64)

            for i in range(num_vectors):
                vector = vectors[i]

                normalized_vector = normalize_2D_vector(vector)

                normalized_vectors[i, 0] = normalized_vector[0]
                normalized_vectors[i, 1] = normalized_vector[1]

            return normalized_vectors

        num_vectors = 10
        test_vectors = np.random.random((num_vectors, 2))
        got = normalize_vectors(num_vectors, test_vectors)
        expected = normalize_vectors.py_func(num_vectors, test_vectors)

        np.testing.assert_almost_equal(expected, got)

    def test_incref_after_cast(self):
        # Issue #1427: when casting a value before returning it, the
        # cast result should be incref'ed, not the original value.
        def f():
            return 0.0, np.zeros(1, dtype=np.int32)

        # Note the return type isn't the same as the tuple type above:
        # the first element is a complex rather than a float.
        cfunc = njit((types.Tuple((types.complex128,
                                   types.Array(types.int32, 1, 'C') )))())(f)
        z, arr = cfunc()
        self.assertPreciseEqual(z, 0j)
        self.assertPreciseEqual(arr, np.zeros(1, dtype=np.int32))

    def test_refct_pruning_issue_1511(self):
        @njit
        def f():
            a = np.ones(10, dtype=np.float64)
            b = np.ones(10, dtype=np.float64)
            return a, b[:]

        a, b = f()
        np.testing.assert_equal(a, b)
        np.testing.assert_equal(a, np.ones(10, dtype=np.float64))

    def test_refct_pruning_issue_1526(self):
        @njit
        def udt(image, x, y):
            next_loc = np.where(image == 1)

            if len(next_loc[0]) == 0:
                y_offset = 1
                x_offset = 1
            else:
                y_offset = next_loc[0][0]
                x_offset = next_loc[1][0]

            next_loc_x = (x - 1) + x_offset
            next_loc_y = (y - 1) + y_offset

            return next_loc_x, next_loc_y

        a = np.array([[1, 0, 1, 0, 1, 0, 0, 1, 0, 0]])
        expect = udt.py_func(a, 1, 6)
        got = udt(a, 1, 6)

        self.assertEqual(expect, got)

    @TestCase.run_test_in_subprocess
    def test_no_nrt_on_njit_decoration(self):
        # Checks that the NRT is not initialized/compiled as a result of
        # decorating a function with `@njit`.
        from numba import njit

        # check the NRT is not initialized.
        self.assertFalse(rtsys._init)

        # decorate
        @njit
        def foo():
            return 123

        # check the NRT is still not initialized
        self.assertFalse(rtsys._init)

        # execute
        self.assertEqual(foo(), foo.py_func())

        # check the NRT is still now initialized as execution has definitely
        # occurred.
        self.assertTrue(rtsys._init)


class TestRefCtPruning(unittest.TestCase):

    sample_llvm_ir = '''
define i32 @"MyFunction"(i8** noalias nocapture %retptr, { i8*, i32 }** noalias nocapture %excinfo, i8* noalias nocapture readnone %env, double %arg.vt.0, double %arg.vt.1, double %arg.vt.2, double %arg.vt.3, double %arg.bounds.0, double %arg.bounds.1, double %arg.bounds.2, double %arg.bounds.3, i8* %arg.xs.0, i8* nocapture readnone %arg.xs.1, i64 %arg.xs.2, i64 %arg.xs.3, double* nocapture readonly %arg.xs.4, i64 %arg.xs.5.0, i64 %arg.xs.6.0, i8* %arg.ys.0, i8* nocapture readnone %arg.ys.1, i64 %arg.ys.2, i64 %arg.ys.3, double* nocapture readonly %arg.ys.4, i64 %arg.ys.5.0, i64 %arg.ys.6.0, i8* %arg.aggs_and_cols.0.0, i8* nocapture readnone %arg.aggs_and_cols.0.1, i64 %arg.aggs_and_cols.0.2, i64 %arg.aggs_and_cols.0.3, i32* nocapture %arg.aggs_and_cols.0.4, i64 %arg.aggs_and_cols.0.5.0, i64 %arg.aggs_and_cols.0.5.1, i64 %arg.aggs_and_cols.0.6.0, i64 %arg.aggs_and_cols.0.6.1) local_unnamed_addr {
entry:
tail call void @NRT_incref(i8* %arg.xs.0)
tail call void @NRT_incref(i8* %arg.ys.0)
tail call void @NRT_incref(i8* %arg.aggs_and_cols.0.0)
%.251 = icmp sgt i64 %arg.xs.5.0, 0
br i1 %.251, label %B42.preheader, label %B160

B42.preheader:                                    ; preds = %entry
%0 = add i64 %arg.xs.5.0, 1
br label %B42

B42:                                              ; preds = %B40.backedge, %B42.preheader
%lsr.iv3 = phi i64 [ %lsr.iv.next, %B40.backedge ], [ %0, %B42.preheader ]
%lsr.iv1 = phi double* [ %scevgep2, %B40.backedge ], [ %arg.xs.4, %B42.preheader ]
%lsr.iv = phi double* [ %scevgep, %B40.backedge ], [ %arg.ys.4, %B42.preheader ]
%.381 = load double, double* %lsr.iv1, align 8
%.420 = load double, double* %lsr.iv, align 8
%.458 = fcmp ole double %.381, %arg.bounds.1
%not..432 = fcmp oge double %.381, %arg.bounds.0
%"$phi82.1.1" = and i1 %.458, %not..432
br i1 %"$phi82.1.1", label %B84, label %B40.backedge

B84:                                              ; preds = %B42
%.513 = fcmp ole double %.420, %arg.bounds.3
%not..487 = fcmp oge double %.420, %arg.bounds.2
%"$phi106.1.1" = and i1 %.513, %not..487
br i1 %"$phi106.1.1", label %B108.endif.endif.endif, label %B40.backedge

B160:                                             ; preds = %B40.backedge, %entry
tail call void @NRT_decref(i8* %arg.ys.0)
tail call void @NRT_decref(i8* %arg.xs.0)
tail call void @NRT_decref(i8* %arg.aggs_and_cols.0.0)
store i8* null, i8** %retptr, align 8
ret i32 0

B108.endif.endif.endif:                           ; preds = %B84
%.575 = fmul double %.381, %arg.vt.0
%.583 = fadd double %.575, %arg.vt.1
%.590 = fptosi double %.583 to i64
%.630 = fmul double %.420, %arg.vt.2
%.638 = fadd double %.630, %arg.vt.3
%.645 = fptosi double %.638 to i64
tail call void @NRT_incref(i8* %arg.aggs_and_cols.0.0)              ; GONE 1
tail call void @NRT_decref(i8* null)                                ; GONE 2
tail call void @NRT_incref(i8* %arg.aggs_and_cols.0.0), !noalias !0 ; GONE 3
%.62.i.i = icmp slt i64 %.645, 0
%.63.i.i = select i1 %.62.i.i, i64 %arg.aggs_and_cols.0.5.0, i64 0
%.64.i.i = add i64 %.63.i.i, %.645
%.65.i.i = icmp slt i64 %.590, 0
%.66.i.i = select i1 %.65.i.i, i64 %arg.aggs_and_cols.0.5.1, i64 0
%.67.i.i = add i64 %.66.i.i, %.590
%.84.i.i = mul i64 %.64.i.i, %arg.aggs_and_cols.0.5.1
%.87.i.i = add i64 %.67.i.i, %.84.i.i
%.88.i.i = getelementptr i32, i32* %arg.aggs_and_cols.0.4, i64 %.87.i.i
%.89.i.i = load i32, i32* %.88.i.i, align 4, !noalias !3
%.99.i.i = add i32 %.89.i.i, 1
store i32 %.99.i.i, i32* %.88.i.i, align 4, !noalias !3
tail call void @NRT_decref(i8* %arg.aggs_and_cols.0.0), !noalias !0 ; GONE 4
tail call void @NRT_decref(i8* %arg.aggs_and_cols.0.0)              ; GONE 5
br label %B40.backedge

B40.backedge:                                     ; preds = %B108.endif.endif.endif, %B84, %B42
%scevgep = getelementptr double, double* %lsr.iv, i64 1
%scevgep2 = getelementptr double, double* %lsr.iv1, i64 1
%lsr.iv.next = add i64 %lsr.iv3, -1
%.294 = icmp sgt i64 %lsr.iv.next, 1
br i1 %.294, label %B42, label %B160
}
    ''' # noqa

    def test_refct_pruning_op_recognize(self):
        input_ir = self.sample_llvm_ir
        input_lines = list(input_ir.splitlines())
        before_increfs = [ln for ln in input_lines if 'NRT_incref' in ln]
        before_decrefs = [ln for ln in input_lines if 'NRT_decref' in ln]

        # prune
        output_ir = nrtopt._remove_redundant_nrt_refct(input_ir)
        output_lines = list(output_ir.splitlines())
        after_increfs = [ln for ln in output_lines if 'NRT_incref' in ln]
        after_decrefs = [ln for ln in output_lines if 'NRT_decref' in ln]

        # check
        self.assertNotEqual(before_increfs, after_increfs)
        self.assertNotEqual(before_decrefs, after_decrefs)

        pruned_increfs = set(before_increfs) - set(after_increfs)
        pruned_decrefs = set(before_decrefs) - set(after_decrefs)

        # the symm difference == or-combined
        combined = pruned_increfs | pruned_decrefs
        self.assertEqual(combined, pruned_increfs ^ pruned_decrefs)
        pruned_lines = '\n'.join(combined)

        # all GONE lines are pruned
        for i in [1, 2, 3, 4, 5]:
            gone = '; GONE {}'.format(i)
            self.assertIn(gone, pruned_lines)
        # no other lines
        self.assertEqual(len(list(pruned_lines.splitlines())), len(combined))

    @unittest.skip("Pass removed as it was buggy. Re-enable when fixed.")
    def test_refct_pruning_with_branches(self):
        '''testcase from #2350'''
        @njit
        def _append_non_na(x, y, agg, field):
            if not np.isnan(field):
                agg[y, x] += 1

        @njit
        def _append(x, y, agg, field):
            if not np.isnan(field):
                if np.isnan(agg[y, x]):
                    agg[y, x] = field
                else:
                    agg[y, x] += field

        @njit
        def append(x, y, agg, field):
            _append_non_na(x, y, agg, field)
            _append(x, y, agg, field)

        # Disable python wrapper to avoid detecting necessary
        # refcount inside it
        @njit(no_cpython_wrapper=True)
        def extend(arr, field):
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    append(j, i, arr, field)

        # Compile
        extend.compile("(f4[:,::1], f4)")

        # Test there are no reference count operations
        llvmir = str(extend.inspect_llvm(extend.signatures[0]))
        refops = list(re.finditer(r'(NRT_incref|NRT_decref)\([^\)]+\)', llvmir))
        self.assertEqual(len(refops), 0)

    @linux_only
    @x86_only
    def test_inline_asm(self):
        """The InlineAsm class from llvmlite.ir has no 'name' attr the refcount
        pruning pass should be tolerant to this"""
        llvm.initialize_native_target()
        llvm.initialize_native_asmprinter()
        llvm.initialize_native_asmparser()

        @intrinsic
        def bar(tyctx, x, y):
            def codegen(cgctx, builder, sig, args):
                (arg_0, arg_1) = args
                fty = ir.FunctionType(ir.IntType(32), [ir.IntType(32),
                                                       ir.IntType(32)])
                mul = builder.asm(fty, "mov $2, $0; imul $1, $0", "=&r,r,r",
                                  (arg_0, arg_1), name="asm_mul",
                                  side_effect=False)
                return impl_ret_untracked(cgctx, builder, sig.return_type, mul)
            return signature(types.int32, types.int32, types.int32), codegen

        @njit(['int32(int32)'])
        def foo(x):
            x += 1
            z = bar(x, 2)
            return z

        self.assertEqual(foo(10), 22) # expect (10 + 1) * 2 = 22


@skip_unless_cffi
class TestNrtExternalCFFI(EnableNRTStatsMixin, TestCase):
    """Testing the use of externally compiled C code that use NRT
    """
    def setUp(self):
        # initialize the NRT (in case the tests are run in isolation)
        cpu_target.target_context
        super(TestNrtExternalCFFI, self).setUp()

    def compile_cffi_module(self, name, source, cdef):
        from cffi import FFI

        ffi = FFI()
        ffi.set_source(name, source, include_dirs=[include_path()])
        ffi.cdef(cdef)
        tmpdir = temp_directory("cffi_test_{}".format(name))
        ffi.compile(tmpdir=tmpdir)
        sys.path.append(tmpdir)
        try:
            mod = import_dynamic(name)
        finally:
            sys.path.remove(tmpdir)

        return ffi, mod

    def get_nrt_api_table(self):
        from cffi import FFI

        ffi = FFI()
        nrt_get_api = ffi.cast("void* (*)()", _nrt_python.c_helpers['get_api'])
        table = nrt_get_api()
        return table

    def test_manage_memory(self):
        name = "{}_test_manage_memory".format(self.__class__.__name__)
        source = r"""
#include <stdio.h>
#include "numba/core/runtime/nrt_external.h"

int status = 0;

void my_dtor(void *ptr) {
    free(ptr);
    status = 0xdead;
}

NRT_MemInfo* test_nrt_api(NRT_api_functions *nrt) {
    void * data = malloc(10);
    NRT_MemInfo *mi = nrt->manage_memory(data, my_dtor);
    nrt->acquire(mi);
    nrt->release(mi);
    status = 0xa110c;
    return mi;
}
        """
        cdef = """
void* test_nrt_api(void *nrt);
extern int status;
        """

        ffi, mod = self.compile_cffi_module(name, source, cdef)
        # Init status is 0
        self.assertEqual(mod.lib.status, 0)
        table = self.get_nrt_api_table()
        out = mod.lib.test_nrt_api(table)
        # status is now 0xa110c
        self.assertEqual(mod.lib.status, 0xa110c)
        mi_addr = int(ffi.cast("size_t", out))
        mi = nrt.MemInfo(mi_addr)
        self.assertEqual(mi.refcount, 1)
        del mi   # force deallocation on mi
        # status is now 0xdead
        self.assertEqual(mod.lib.status, 0xdead)

    def test_allocate(self):
        name = "{}_test_allocate".format(self.__class__.__name__)
        source = r"""
#include <stdio.h>
#include "numba/core/runtime/nrt_external.h"

NRT_MemInfo* test_nrt_api(NRT_api_functions *nrt, size_t n) {
    size_t *data = NULL;
    NRT_MemInfo *mi = nrt->allocate(n);
    data = nrt->get_data(mi);
    data[0] = 0xded;
    data[1] = 0xabc;
    data[2] = 0xdef;
    return mi;
}
        """
        cdef = "void* test_nrt_api(void *nrt, size_t n);"
        ffi, mod = self.compile_cffi_module(name, source, cdef)

        table = self.get_nrt_api_table()

        numbytes = 3 * np.dtype(np.intp).itemsize
        out = mod.lib.test_nrt_api(table, numbytes)

        mi_addr = int(ffi.cast("size_t", out))
        mi = nrt.MemInfo(mi_addr)
        self.assertEqual(mi.refcount, 1)

        buffer = ffi.buffer(ffi.cast("char [{}]".format(numbytes), mi.data))
        arr = np.ndarray(shape=(3,), dtype=np.intp, buffer=buffer)
        np.testing.assert_equal(arr, [0xded, 0xabc, 0xdef])

    def test_get_api(self):
        from cffi import FFI

        @njit
        def test_nrt_api():
            return NRT_get_api()

        ffi = FFI()
        expect = int(ffi.cast('size_t', self.get_nrt_api_table()))
        got = test_nrt_api()
        self.assertEqual(expect, got)


class TestNrtStatistics(TestCase):

    def setUp(self):
        # Store the current stats state
        self.__stats_state = _nrt_python.memsys_stats_enabled()

    def tearDown(self):
        # Set stats state back to whatever it was before the test ran
        if self.__stats_state:
            _nrt_python.memsys_enable_stats()
        else:
            _nrt_python.memsys_disable_stats()

    def test_stats_env_var_explicit_on(self):
        # Checks that explicitly turning the stats on via the env var works.
        src = """if 1:
        from numba import njit
        import numpy as np
        from numba.core.runtime import rtsys, _nrt_python
        from numba.core.registry import cpu_target

        @njit
        def foo():
            return np.arange(10)[0]

        # initialize the NRT before use
        rtsys.initialize(cpu_target.target_context)
        assert _nrt_python.memsys_stats_enabled()
        orig_stats = rtsys.get_allocation_stats()
        foo()
        new_stats = rtsys.get_allocation_stats()
        total_alloc = new_stats.alloc - orig_stats.alloc
        total_free = new_stats.free - orig_stats.free
        total_mi_alloc = new_stats.mi_alloc - orig_stats.mi_alloc
        total_mi_free = new_stats.mi_free - orig_stats.mi_free

        expected = 1
        assert total_alloc == expected
        assert total_free == expected
        assert total_mi_alloc == expected
        assert total_mi_free == expected
        """
        # Check env var explicitly being set works
        env = os.environ.copy()
        env['NUMBA_NRT_STATS'] = "1"
        run_in_subprocess(src, env=env)

    def check_env_var_off(self, env):

        src = """if 1:
        from numba import njit
        import numpy as np
        from numba.core.runtime import rtsys, _nrt_python

        @njit
        def foo():
            return np.arange(10)[0]

        assert _nrt_python.memsys_stats_enabled() == False
        try:
            rtsys.get_allocation_stats()
        except RuntimeError as e:
            assert "NRT stats are disabled." in str(e)
        """
        run_in_subprocess(src, env=env)

    def test_stats_env_var_explicit_off(self):
        # Checks that explicitly turning the stats off via the env var works.
        env = os.environ.copy()
        env['NUMBA_NRT_STATS'] = "0"
        self.check_env_var_off(env)

    def test_stats_env_var_default_off(self):
        # Checks that the env var not being set is the same as "off", i.e.
        # default for Numba is off.
        env = os.environ.copy()
        env.pop('NUMBA_NRT_STATS', None)
        self.check_env_var_off(env)

    def test_stats_status_toggle(self):

        @njit
        def foo():
            tmp = np.ones(3)
            return np.arange(5 * tmp[0])

        # Switch on stats
        _nrt_python.memsys_enable_stats()
        # check the stats are on
        self.assertTrue(_nrt_python.memsys_stats_enabled())

        for i in range(2):
            # capture the stats state
            stats_1 = rtsys.get_allocation_stats()
            # Switch off stats
            _nrt_python.memsys_disable_stats()
            # check the stats are off
            self.assertFalse(_nrt_python.memsys_stats_enabled())
            # run something that would move the counters were they enabled
            foo()
            # Switch on stats
            _nrt_python.memsys_enable_stats()
            # check the stats are on
            self.assertTrue(_nrt_python.memsys_stats_enabled())
            # capture the stats state (should not have changed)
            stats_2 = rtsys.get_allocation_stats()
            # run something that will move the counters
            foo()
            # capture the stats state (should have changed)
            stats_3 = rtsys.get_allocation_stats()
            # check stats_1 == stats_2
            self.assertEqual(stats_1, stats_2)
            # check stats_2 < stats_3
            self.assertLess(stats_2, stats_3)

    def test_rtsys_stats_query_raises_exception_when_disabled(self):
        # Checks that the standard rtsys.get_allocation_stats() query raises
        # when stats counters are turned off.

        _nrt_python.memsys_disable_stats()
        self.assertFalse(_nrt_python.memsys_stats_enabled())

        with self.assertRaises(RuntimeError) as raises:
            rtsys.get_allocation_stats()

        self.assertIn("NRT stats are disabled.", str(raises.exception))

    def test_nrt_explicit_stats_query_raises_exception_when_disabled(self):
        # Checks the various memsys_get_stats functions raise if queried when
        # the stats counters are disabled.
        method_variations = ('alloc', 'free', 'mi_alloc', 'mi_free')
        for meth in method_variations:
            stats_func = getattr(_nrt_python, f'memsys_get_stats_{meth}')
            with self.subTest(stats_func=stats_func):
                # Turn stats off
                _nrt_python.memsys_disable_stats()
                self.assertFalse(_nrt_python.memsys_stats_enabled())
                with self.assertRaises(RuntimeError) as raises:
                    stats_func()
                self.assertIn("NRT stats are disabled.", str(raises.exception))


if __name__ == '__main__':
    unittest.main()
