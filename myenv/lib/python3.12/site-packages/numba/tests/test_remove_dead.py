#
# Copyright (c) 2017 Intel Corporation
# SPDX-License-Identifier: BSD-2-Clause
#

import numba
import numba.parfors.parfor
from numba import njit
from numba.core import ir_utils
from numba.core import types, ir,  compiler
from numba.core.registry import cpu_target
from numba.core.ir_utils import (copy_propagate, apply_copy_propagate,
                            get_name_var_table, remove_dels, remove_dead,
                            remove_call_handlers, alias_func_extensions)
from numba.core.typed_passes import type_inference_stage
from numba.core.compiler_machinery import FunctionPass, register_pass, PassManager
from numba.core.untyped_passes import (ExtractByteCode, TranslateByteCode, FixupArgs,
                             IRProcessing, DeadBranchPrune,
                             RewriteSemanticConstants, GenericRewrites,
                             WithLifting, PreserveIR, InlineClosureLikes)

from numba.core.typed_passes import (NopythonTypeInference, AnnotateTypes,
                           NopythonRewrites, PreParforPass, ParforPass,
                           DumpParforDiagnostics, NativeLowering,
                           IRLegalization, NoPythonBackend, NativeLowering)
import numpy as np
from numba.tests.support import skip_parfors_unsupported, needs_blas
import unittest


def test_will_propagate(b, z, w):
    x1 = 3
    x = x1
    if b > 0:
        y = z + w
    else:
        y = 0
    a = 2 * x
    return a < b

def null_func(a,b,c,d):
    False

@numba.njit
def dummy_aliased_func(A):
    return A

def alias_ext_dummy_func(lhs_name, args, alias_map, arg_aliases):
    ir_utils._add_alias(lhs_name, args[0].name, alias_map, arg_aliases)

def findLhsAssign(func_ir, var):
    for label, block in func_ir.blocks.items():
        for i, inst in enumerate(block.body):
            if isinstance(inst, ir.Assign) and inst.target.name==var:
                return True

    return False

class TestRemoveDead(unittest.TestCase):

    _numba_parallel_test_ = False

    def compile_parallel(self, func, arg_types):
        return njit(arg_types, parallel=True, fastmath=True)(func)

    def test1(self):
        typingctx = cpu_target.typing_context
        targetctx = cpu_target.target_context
        test_ir = compiler.run_frontend(test_will_propagate)

        typingctx.refresh()
        targetctx.refresh()
        args = (types.int64, types.int64, types.int64)
        typemap, _, calltypes, _ = type_inference_stage(typingctx, targetctx, test_ir, args, None)
        remove_dels(test_ir.blocks)
        in_cps, out_cps = copy_propagate(test_ir.blocks, typemap)
        apply_copy_propagate(test_ir.blocks, in_cps, get_name_var_table(test_ir.blocks), typemap, calltypes)

        remove_dead(test_ir.blocks, test_ir.arg_names, test_ir)
        self.assertFalse(findLhsAssign(test_ir, "x"))

    def test2(self):
        def call_np_random_seed():
            np.random.seed(2)

        def seed_call_exists(func_ir):
            for inst in func_ir.blocks[0].body:
                if (isinstance(inst, ir.Assign) and
                    isinstance(inst.value, ir.Expr) and
                    inst.value.op == 'call' and
                    func_ir.get_definition(inst.value.func).attr == 'seed'):
                    return True
            return False

        test_ir = compiler.run_frontend(call_np_random_seed)
        remove_dead(test_ir.blocks, test_ir.arg_names, test_ir)
        self.assertTrue(seed_call_exists(test_ir))

    def run_array_index_test(self, func):
        A1 = np.arange(6).reshape(2,3)
        A2 = A1.copy()
        i = 0
        pfunc = self.compile_parallel(func, (numba.typeof(A1), numba.typeof(i)))

        func(A1, i)
        pfunc(A2, i)
        np.testing.assert_array_equal(A1, A2)

    def test_alias_ravel(self):
        def func(A, i):
            B = A.ravel()
            B[i] = 3

        self.run_array_index_test(func)

    def test_alias_flat(self):
        def func(A, i):
            B = A.flat
            B[i] = 3

        self.run_array_index_test(func)

    def test_alias_transpose1(self):
        def func(A, i):
            B = A.T
            B[i,0] = 3

        self.run_array_index_test(func)

    def test_alias_transpose2(self):
        def func(A, i):
            B = A.transpose()
            B[i,0] = 3

        self.run_array_index_test(func)

    def test_alias_transpose3(self):
        def func(A, i):
            B = np.transpose(A)
            B[i,0] = 3

        self.run_array_index_test(func)

    @skip_parfors_unsupported
    @needs_blas
    def test_alias_ctypes(self):
        # use xxnrm2 to test call a C function with ctypes
        from numba.np.linalg import _BLAS
        xxnrm2 = _BLAS().numba_xxnrm2(types.float64)

        def remove_dead_xxnrm2(rhs, lives, call_list):
            if call_list == [xxnrm2]:
                return rhs.args[4].name not in lives
            return False

        # adding this handler has no-op effect since this function won't match
        # anything else but it's a bit cleaner to save the state and recover
        old_remove_handlers = remove_call_handlers[:]
        remove_call_handlers.append(remove_dead_xxnrm2)

        def func(ret):
            a = np.ones(4)
            xxnrm2(100, 4, a.ctypes, 1, ret.ctypes)

        A1 = np.zeros(1)
        A2 = A1.copy()

        try:
            pfunc = self.compile_parallel(func, (numba.typeof(A1),))
            numba.njit(func)(A1)
            pfunc(A2)
        finally:
            # recover global state
            remove_call_handlers[:] = old_remove_handlers

        self.assertEqual(A1[0], A2[0])

    def test_alias_reshape1(self):
        def func(A, i):
            B = np.reshape(A, (3,2))
            B[i,0] = 3

        self.run_array_index_test(func)

    def test_alias_reshape2(self):
        def func(A, i):
            B = A.reshape(3,2)
            B[i,0] = 3

        self.run_array_index_test(func)

    def test_alias_func_ext(self):
        def func(A, i):
            B = dummy_aliased_func(A)
            B[i, 0] = 3

        # save global state
        old_ext_handlers = alias_func_extensions.copy()
        try:
            alias_func_extensions[('dummy_aliased_func',
                'numba.tests.test_remove_dead')] = alias_ext_dummy_func
            self.run_array_index_test(func)
        finally:
            # recover global state
            ir_utils.alias_func_extensions = old_ext_handlers

    def test_rm_dead_rhs_vars(self):
        """make sure lhs variable of assignment is considered live if used in
        rhs (test for #6715).
        """
        def func():
            for i in range(3):
                a = (lambda j: j)(i)
                a = np.array(a)
            return a

        self.assertEqual(func(), numba.njit(func)())

    @skip_parfors_unsupported
    def test_alias_parfor_extension(self):
        """Make sure aliases are considered in remove dead extension for
        parfors.
        """
        def func():
            n = 11
            numba.parfors.parfor.init_prange()
            A = np.empty(n)
            B = A  # create alias to A
            for i in numba.prange(n):
                A[i] = i

            return B

        @register_pass(analysis_only=False, mutates_CFG=True)
        class LimitedParfor(FunctionPass):
            _name = "limited_parfor"

            def __init__(self):
                FunctionPass.__init__(self)

            def run_pass(self, state):
                parfor_pass = numba.parfors.parfor.ParforPass(
                    state.func_ir,
                    state.typemap,
                    state.calltypes,
                    state.return_type,
                    state.typingctx,
                    state.flags.auto_parallel,
                    state.flags,
                    state.metadata,
                    state.parfor_diagnostics
                )
                remove_dels(state.func_ir.blocks)
                parfor_pass.array_analysis.run(state.func_ir.blocks)
                parfor_pass._convert_loop(state.func_ir.blocks)
                remove_dead(state.func_ir.blocks,
                            state.func_ir.arg_names,
                            state.func_ir,
                            state.typemap)
                numba.parfors.parfor.get_parfor_params(state.func_ir.blocks,
                                                parfor_pass.options.fusion,
                                                parfor_pass.nested_fusion_info)
                return True

        class TestPipeline(compiler.Compiler):
            """Test pipeline that just converts prange() to parfor and calls
            remove_dead(). Copy propagation can replace B in the example code
            which this pipeline avoids.
            """
            def define_pipelines(self):
                name = 'test parfor aliasing'
                pm = PassManager(name)
                pm.add_pass(TranslateByteCode, "analyzing bytecode")
                pm.add_pass(FixupArgs, "fix up args")
                pm.add_pass(IRProcessing, "processing IR")
                pm.add_pass(WithLifting, "Handle with contexts")
                # pre typing
                if not self.state.flags.no_rewrites:
                    pm.add_pass(GenericRewrites, "nopython rewrites")
                    pm.add_pass(RewriteSemanticConstants, "rewrite semantic constants")
                    pm.add_pass(DeadBranchPrune, "dead branch pruning")
                pm.add_pass(InlineClosureLikes,
                            "inline calls to locally defined closures")
                # typing
                pm.add_pass(NopythonTypeInference, "nopython frontend")

                # lower
                pm.add_pass(NativeLowering, "native lowering")
                pm.add_pass(NoPythonBackend, "nopython mode backend")
                pm.finalize()
                return [pm]

        test_res = numba.jit(pipeline_class=TestPipeline)(func)()
        py_res = func()
        np.testing.assert_array_equal(test_res, py_res)


if __name__ == "__main__":
    unittest.main()
