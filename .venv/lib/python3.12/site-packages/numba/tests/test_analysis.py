# Tests numba.analysis functions
import collections
import types as pytypes

import numpy as np
from numba.core.compiler import run_frontend, Flags, StateDict
from numba import jit, njit, literal_unroll
from numba.core import types, errors, ir, rewrites, ir_utils, cpu
from numba.core import postproc
from numba.core.inline_closurecall import InlineClosureCallPass
from numba.tests.support import (TestCase, MemoryLeakMixin, SerialMixin,
                                 IRPreservingTestPipeline)
from numba.core.analysis import dead_branch_prune, rewrite_semantic_constants
from numba.core.untyped_passes import (ReconstructSSA, TranslateByteCode,
                                       IRProcessing, DeadBranchPrune,
                                       PreserveIR)
from numba.core.compiler import DefaultPassBuilder, CompilerBase, PassManager
from numba.core.utils import PYVERSION


_GLOBAL = 123

enable_pyobj_flags = Flags()
enable_pyobj_flags.enable_pyobject = True


def compile_to_ir(func):
    func_ir = run_frontend(func)
    state = StateDict()
    state.func_ir = func_ir
    state.typemap = None
    state.calltypes = None
    # Transform to SSA
    ReconstructSSA().run_pass(state)
    # call this to get print etc rewrites
    rewrites.rewrite_registry.apply('before-inference', state)
    return func_ir


class TestBranchPruneBase(MemoryLeakMixin, TestCase):
    """
    Tests branch pruning
    """
    _DEBUG = False

    # find *all* branches
    def find_branches(self, the_ir):
        branches = []
        for blk in the_ir.blocks.values():
            tmp = [_ for _ in blk.find_insts(cls=ir.Branch)]
            branches.extend(tmp)
        return branches

    def assert_prune(self, func, args_tys, prune, *args, **kwargs):
        # This checks that the expected pruned branches have indeed been pruned.
        # func is a python function to assess
        # args_tys is the numba types arguments tuple
        # prune arg is a list, one entry per branch. The value in the entry is
        # encoded as follows:
        # True: using constant inference only, the True branch will be pruned
        # False: using constant inference only, the False branch will be pruned
        # None: under no circumstances should this branch be pruned
        # *args: the argument instances to pass to the function to check
        #        execution is still valid post transform
        # **kwargs:
        #        - flags: args to pass to `jit` default is `nopython=True`,
        #          e.g. permits use of e.g. object mode.

        func_ir = compile_to_ir(func)
        before = func_ir.copy()
        if self._DEBUG:
            print("=" * 80)
            print("before inline")
            func_ir.dump()

        # run closure inlining to ensure that nonlocals in closures are visible
        inline_pass = InlineClosureCallPass(func_ir,
                                            cpu.ParallelOptions(False),)
        inline_pass.run()

        # Remove all Dels, and re-run postproc
        post_proc = postproc.PostProcessor(func_ir)
        post_proc.run()

        rewrite_semantic_constants(func_ir, args_tys)
        if self._DEBUG:
            print("=" * 80)
            print("before prune")
            func_ir.dump()

        dead_branch_prune(func_ir, args_tys)

        after = func_ir
        if self._DEBUG:
            print("after prune")
            func_ir.dump()

        before_branches = self.find_branches(before)
        self.assertEqual(len(before_branches), len(prune))

        # what is expected to be pruned
        expect_removed = []
        for idx, prune in enumerate(prune):
            branch = before_branches[idx]
            if prune is True:
                expect_removed.append(branch.truebr)
            elif prune is False:
                expect_removed.append(branch.falsebr)
            elif prune is None:
                pass  # nothing should be removed!
            elif prune == 'both':
                expect_removed.append(branch.falsebr)
                expect_removed.append(branch.truebr)
            else:
                assert 0, "unreachable"

        # compare labels
        original_labels = set([_ for _ in before.blocks.keys()])
        new_labels = set([_ for _ in after.blocks.keys()])
        # assert that the new labels are precisely the original less the
        # expected pruned labels
        try:
            self.assertEqual(new_labels, original_labels - set(expect_removed))
        except AssertionError as e:
            print("new_labels", sorted(new_labels))
            print("original_labels", sorted(original_labels))
            print("expect_removed", sorted(expect_removed))
            raise e

        supplied_flags = kwargs.pop('flags', {'nopython': True})
        # NOTE: original testing used `compile_isolated` hence use of `cres`.
        cres = jit(args_tys, **supplied_flags)(func).overloads[args_tys]
        if args is None:
            res = cres.entry_point()
            expected = func()
        else:
            res = cres.entry_point(*args)
            expected = func(*args)
        self.assertEqual(res, expected)


class TestBranchPrune(TestBranchPruneBase, SerialMixin):

    def test_single_if(self):

        def impl(x):
            if 1 == 0:
                return 3.14159

        self.assert_prune(impl, (types.NoneType('none'),), [True], None)

        def impl(x):
            if 1 == 1:
                return 3.14159

        self.assert_prune(impl, (types.NoneType('none'),), [False], None)

        def impl(x):
            if x is None:
                return 3.14159

        self.assert_prune(impl, (types.NoneType('none'),), [False], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True], 10)

        def impl(x):
            if x == 10:
                return 3.14159

        self.assert_prune(impl, (types.NoneType('none'),), [True], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [None], 10)

        def impl(x):
            if x == 10:
                z = 3.14159  # noqa: F841 # no effect

        self.assert_prune(impl, (types.NoneType('none'),), [True], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [None], 10)

        def impl(x):
            z = None
            y = z
            if x == y:
                return 100

        self.assert_prune(impl, (types.NoneType('none'),), [False], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True], 10)

    def test_single_if_else(self):

        def impl(x):
            if x is None:
                return 3.14159
            else:
                return 1.61803

        self.assert_prune(impl, (types.NoneType('none'),), [False], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True], 10)

    def test_single_if_const_val(self):

        def impl(x):
            if x == 100:
                return 3.14159

        self.assert_prune(impl, (types.NoneType('none'),), [True], None)
        self.assert_prune(impl, (types.IntegerLiteral(100),), [None], 100)

        def impl(x):
            # switch the condition order
            if 100 == x:
                return 3.14159

        self.assert_prune(impl, (types.NoneType('none'),), [True], None)
        self.assert_prune(impl, (types.IntegerLiteral(100),), [None], 100)

    def test_single_if_else_two_const_val(self):

        def impl(x, y):
            if x == y:
                return 3.14159
            else:
                return 1.61803

        self.assert_prune(impl, (types.IntegerLiteral(100),) * 2, [None], 100,
                          100)
        self.assert_prune(impl, (types.NoneType('none'),) * 2, [False], None,
                          None)
        self.assert_prune(impl, (types.IntegerLiteral(100),
                                 types.NoneType('none'),), [True], 100, None)
        self.assert_prune(impl, (types.IntegerLiteral(100),
                                 types.IntegerLiteral(1000)), [None], 100, 1000)

    def test_single_if_else_w_following_undetermined(self):

        def impl(x):
            x_is_none_work = False
            if x is None:
                x_is_none_work = True
            else:
                dead = 7  # noqa: F841 # no effect

            if x_is_none_work:
                y = 10
            else:
                y = -3
            return y

        self.assert_prune(impl, (types.NoneType('none'),), [False, None], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True, None], 10)

        def impl(x):
            x_is_none_work = False
            if x is None:
                x_is_none_work = True
            else:
                pass

            if x_is_none_work:
                y = 10
            else:
                y = -3
            return y

        # Python 3.10 creates a block with a NOP in it for the `pass` which
        # means it gets pruned.
        self.assert_prune(impl, (types.NoneType('none'),), [False, None],
                          None)

        self.assert_prune(impl, (types.IntegerLiteral(10),), [True, None], 10)

    def test_double_if_else_rt_const(self):

        def impl(x):
            one_hundred = 100
            x_is_none_work = 4
            if x is None:
                x_is_none_work = 100
            else:
                dead = 7  # noqa: F841 # no effect

            if x_is_none_work == one_hundred:
                y = 10
            else:
                y = -3

            return y, x_is_none_work

        self.assert_prune(impl, (types.NoneType('none'),), [False, None], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True, None], 10)

    def test_double_if_else_non_literal_const(self):

        def impl(x):
            one_hundred = 100
            if x == one_hundred:
                y = 3.14159
            else:
                y = 1.61803
            return y

        # no prune as compilation specialization on literal value not permitted
        self.assert_prune(impl, (types.IntegerLiteral(10),), [None], 10)
        self.assert_prune(impl, (types.IntegerLiteral(100),), [None], 100)

    def test_single_two_branches_same_cond(self):

        def impl(x):
            if x is None:
                y = 10
            else:
                y = 40

            if x is not None:
                z = 100
            else:
                z = 400

            return z, y

        self.assert_prune(impl, (types.NoneType('none'),), [False, True], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True, False], 10)

    def test_cond_is_kwarg_none(self):

        def impl(x=None):
            if x is None:
                y = 10
            else:
                y = 40

            if x is not None:
                z = 100
            else:
                z = 400

            return z, y

        self.assert_prune(impl, (types.Omitted(None),),
                          [False, True], None)
        self.assert_prune(impl, (types.NoneType('none'),), [False, True], None)
        self.assert_prune(impl, (types.IntegerLiteral(10),), [True, False], 10)

    def test_cond_is_kwarg_value(self):

        def impl(x=1000):
            if x == 1000:
                y = 10
            else:
                y = 40

            if x != 1000:
                z = 100
            else:
                z = 400

            return z, y

        self.assert_prune(impl, (types.Omitted(1000),), [None, None], 1000)
        self.assert_prune(impl, (types.IntegerLiteral(1000),), [None, None],
                          1000)
        self.assert_prune(impl, (types.IntegerLiteral(0),), [None, None], 0)
        self.assert_prune(impl, (types.NoneType('none'),), [True, False], None)

    def test_cond_rewrite_is_correct(self):
        # this checks that when a condition is replaced, it is replace by a
        # true/false bit that correctly represents the evaluated condition
        def fn(x):
            if x is None:
                return 10
            return 12

        def check(func, arg_tys, bit_val):
            func_ir = compile_to_ir(func)

            # check there is 1 branch
            before_branches = self.find_branches(func_ir)
            self.assertEqual(len(before_branches), 1)

            # check the condition in the branch is a binop
            pred_var = before_branches[0].cond
            pred_defn = ir_utils.get_definition(func_ir, pred_var)
            self.assertEqual(pred_defn.op, 'call')
            condition_var = pred_defn.args[0]
            condition_op = ir_utils.get_definition(func_ir, condition_var)
            self.assertEqual(condition_op.op, 'binop')

            # do the prune, this should kill the dead branch and rewrite the
            #'condition to a true/false const bit
            if self._DEBUG:
                print("=" * 80)
                print("before prune")
                func_ir.dump()
            dead_branch_prune(func_ir, arg_tys)
            if self._DEBUG:
                print("=" * 80)
                print("after prune")
                func_ir.dump()

            # after mutation, the condition should be a const value `bit_val`
            new_condition_defn = ir_utils.get_definition(func_ir, condition_var)
            self.assertTrue(isinstance(new_condition_defn, ir.Const))
            self.assertEqual(new_condition_defn.value, bit_val)

        check(fn, (types.NoneType('none'),), 1)
        check(fn, (types.IntegerLiteral(10),), 0)

    def test_global_bake_in(self):

        def impl(x):
            if _GLOBAL == 123:
                return x
            else:
                return x + 10

        self.assert_prune(impl, (types.IntegerLiteral(1),), [False], 1)

        global _GLOBAL
        tmp = _GLOBAL

        try:
            _GLOBAL = 5

            def impl(x):
                if _GLOBAL == 123:
                    return x
                else:
                    return x + 10

            self.assert_prune(impl, (types.IntegerLiteral(1),), [True], 1)
        finally:
            _GLOBAL = tmp

    def test_freevar_bake_in(self):

        _FREEVAR = 123

        def impl(x):
            if _FREEVAR == 123:
                return x
            else:
                return x + 10

        self.assert_prune(impl, (types.IntegerLiteral(1),), [False], 1)

        _FREEVAR = 12

        def impl(x):
            if _FREEVAR == 123:
                return x
            else:
                return x + 10

        self.assert_prune(impl, (types.IntegerLiteral(1),), [True], 1)

    def test_redefined_variables_are_not_considered_in_prune(self):
        # see issue #4163, checks that if a variable that is an argument is
        # redefined in the user code it is not considered const

        def impl(array, a=None):
            if a is None:
                a = 0
            if a < 0:
                return 10
            return 30

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.NoneType('none'),),
                          [None, None],
                          np.zeros((2, 3)), None)

    def test_comparison_operators(self):
        # see issue #4163, checks that a variable that is an argument and has
        # value None survives TypeError from invalid comparison which should be
        # dead

        def impl(array, a=None):
            x = 0
            if a is None:
                return 10 # dynamic exec would return here
            # static analysis requires that this is executed with a=None,
            # hence TypeError
            if a < 0:
                return 20
            return x

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.NoneType('none'),),
                          [False, 'both'],
                          np.zeros((2, 3)), None)

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.float64,),
                          [None, None],
                          np.zeros((2, 3)), 12.)

    def test_redefinition_analysis_same_block(self):
        # checks that a redefinition in a block with prunable potential doesn't
        # break

        def impl(array, x, a=None):
            b = 2
            if x < 4:
                b = 12
            if a is None: # known true
                a = 7 # live
            else:
                b = 15 # dead
            if a < 0: # valid as a result of the redefinition of 'a'
                return 10
            return 30 + b + a

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.float64, types.NoneType('none'),),
                          [None, False, None],
                          np.zeros((2, 3)), 1., None)

    def test_redefinition_analysis_different_block_can_exec(self):
        # checks that a redefinition in a block that may be executed prevents
        # pruning

        def impl(array, x, a=None):
            b = 0
            if x > 5:
                a = 11 # a redefined, cannot tell statically if this will exec
            if x < 4:
                b = 12
            if a is None: # cannot prune, cannot determine if re-defn occurred
                b += 5
            else:
                b += 7
                if a < 0:
                    return 10
            return 30 + b

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.float64, types.NoneType('none'),),
                          [None, None, None, None],
                          np.zeros((2, 3)), 1., None)

    def test_redefinition_analysis_different_block_cannot_exec(self):
        # checks that a redefinition in a block guarded by something that
        # has prune potential

        def impl(array, x=None, a=None):
            b = 0
            if x is not None:
                a = 11
            if a is None:
                b += 5
            else:
                b += 7
            return 30 + b

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.NoneType('none'), types.NoneType('none')),
                          [True, None],
                          np.zeros((2, 3)), None, None)

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.NoneType('none'), types.float64),
                          [True, None],
                          np.zeros((2, 3)), None, 1.2)

        self.assert_prune(impl,
                          (types.Array(types.float64, 2, 'C'),
                           types.float64, types.NoneType('none')),
                          [None, None],
                          np.zeros((2, 3)), 1.2, None)

    def test_closure_and_nonlocal_can_prune(self):
        # Closures must be inlined ahead of branch pruning in case nonlocal
        # is used. See issue #6585.
        def impl():
            x = 1000

            def closure():
                nonlocal x
                x = 0

            closure()

            if x == 0:
                return True
            else:
                return False

        self.assert_prune(impl, (), [False,],)

    def test_closure_and_nonlocal_cannot_prune(self):
        # Closures must be inlined ahead of branch pruning in case nonlocal
        # is used. See issue #6585.
        def impl(n):
            x = 1000

            def closure(t):
                nonlocal x
                x = t

            closure(n)

            if x == 0:
                return True
            else:
                return False

        self.assert_prune(impl, (types.int64,), [None,], 1)


class TestBranchPrunePredicates(TestBranchPruneBase, SerialMixin):
    # Really important thing to remember... the branch on predicates end up as
    # POP_JUMP_IF_<bool> and the targets are backwards compared to normal, i.e.
    # the true condition is far jump and the false the near i.e. `if x` would
    # end up in Numba IR as e.g. `branch x 10, 6`.

    _TRUTHY = (1, "String", True, 7.4, 3j)
    _FALSEY = (0, "", False, 0.0, 0j, None)

    def _literal_const_sample_generator(self, pyfunc, consts):
        """
        This takes a python function, pyfunc, and manipulates its co_const
        __code__ member to create a new function with different co_consts as
        supplied in argument consts.

        consts is a dict {index: value} of co_const tuple index to constant
        value used to update a pyfunc clone's co_const.
        """
        pyfunc_code = pyfunc.__code__

        # translate consts spec to update the constants
        co_consts = {k: v for k, v in enumerate(pyfunc_code.co_consts)}
        for k, v in consts.items():
            co_consts[k] = v
        new_consts = tuple([v for _, v in sorted(co_consts.items())])

        # create code object with mutation
        new_code = pyfunc_code.replace(co_consts=new_consts)

        # get function
        return pytypes.FunctionType(new_code, globals())

    def test_literal_const_code_gen(self):
        def impl(x):
            _CONST1 = "PLACEHOLDER1"
            if _CONST1:
                return 3.14159
            else:
                _CONST2 = "PLACEHOLDER2"
            return _CONST2 + 4

        if PYVERSION in ((3, 14), ):
            # The order of the __code__.co_consts changes with 3.14
            new = self._literal_const_sample_generator(impl, {0:0, 2:20})
        elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            new = self._literal_const_sample_generator(impl, {1:0, 3:20})
        else:
            raise NotImplementedError(PYVERSION)

        iconst = impl.__code__.co_consts
        nconst = new.__code__.co_consts

        if PYVERSION in ((3, 14), ):
            self.assertEqual(iconst, ("PLACEHOLDER1", 3.14159,
                                      "PLACEHOLDER2"))
            self.assertEqual(nconst, (0, 3.14159,  20))
        elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            self.assertEqual(iconst, (None, "PLACEHOLDER1", 3.14159,
                                      "PLACEHOLDER2", 4))
            self.assertEqual(nconst, (None, 0, 3.14159,  20, 4))
        else:
            raise NotImplementedError(PYVERSION)

        self.assertEqual(impl(None), 3.14159)
        self.assertEqual(new(None), 24)

    def test_single_if_const(self):

        def impl(x):
            _CONST1 = "PLACEHOLDER1"
            if _CONST1:
                return 3.14159

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                if PYVERSION in ((3, 14), ):
                    # The order of the __code__.co_consts changes with 3.14
                    func = self._literal_const_sample_generator(impl,
                                                                {0: const})
                elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
                    func = self._literal_const_sample_generator(impl,
                                                                {1: const})
                else:
                    raise NotImplementedError(PYVERSION)

                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_negate_const(self):

        def impl(x):
            _CONST1 = "PLACEHOLDER1"
            if not _CONST1:
                return 3.14159

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                if PYVERSION in ((3, 14), ):
                    # The order of the __code__.co_consts changes with 3.14
                    func = self._literal_const_sample_generator(impl,
                                                                {0: const})
                elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
                    func = self._literal_const_sample_generator(impl,
                                                                {1: const})
                else:
                    raise NotImplementedError(PYVERSION)

                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_else_const(self):

        def impl(x):
            _CONST1 = "PLACEHOLDER1"
            if _CONST1:
                return 3.14159
            else:
                return 1.61803

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                if PYVERSION in ((3, 14), ):
                    # The order of the __code__.co_consts changes with 3.14
                    func = self._literal_const_sample_generator(impl,
                                                                {0: const})
                elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
                    func = self._literal_const_sample_generator(impl,
                                                                {1: const})
                else:
                    raise NotImplementedError(PYVERSION)

                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_else_negate_const(self):

        def impl(x):
            _CONST1 = "PLACEHOLDER1"
            if not _CONST1:
                return 3.14159
            else:
                return 1.61803

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                if PYVERSION in ((3, 14), ):
                    # The order of the __code__.co_consts changes with 3.14
                    func = self._literal_const_sample_generator(impl,
                                                                {0: const})
                elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
                    func = self._literal_const_sample_generator(impl,
                                                                {1: const})
                else:
                    raise NotImplementedError(PYVERSION)

                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_freevar(self):
        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                def func(x):
                    if const:
                        return 3.14159, const
                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_negate_freevar(self):
        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                def func(x):
                    if not const:
                        return 3.14159, const
                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_else_freevar(self):
        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                def func(x):
                    if const:
                        return 3.14159, const
                    else:
                        return 1.61803, const
                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_else_negate_freevar(self):
        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for const in c_inp:

                def func(x):
                    if not const:
                        return 3.14159, const
                    else:
                        return 1.61803, const
                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    # globals in this section have absurd names after their test usecase names
    # so as to prevent collisions and permit tests to run in parallel
    def test_single_if_global(self):
        global c_test_single_if_global

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for c in c_inp:
                c_test_single_if_global = c

                def func(x):
                    if c_test_single_if_global:
                        return 3.14159, c_test_single_if_global

                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_negate_global(self):
        global c_test_single_if_negate_global

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for c in c_inp:
                c_test_single_if_negate_global = c

                def func(x):
                    if c_test_single_if_negate_global:
                        return 3.14159, c_test_single_if_negate_global

                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_else_global(self):
        global c_test_single_if_else_global

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for c in c_inp:
                c_test_single_if_else_global = c

                def func(x):
                    if c_test_single_if_else_global:
                        return 3.14159, c_test_single_if_else_global
                    else:
                        return 1.61803, c_test_single_if_else_global
                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_single_if_else_negate_global(self):
        global c_test_single_if_else_negate_global

        for c_inp, prune in (self._TRUTHY, False), (self._FALSEY, True):
            for c in c_inp:
                c_test_single_if_else_negate_global = c

                def func(x):
                    if not c_test_single_if_else_negate_global:
                        return 3.14159, c_test_single_if_else_negate_global
                    else:
                        return 1.61803, c_test_single_if_else_negate_global
                self.assert_prune(func, (types.NoneType('none'),), [prune],
                                  None)

    def test_issue_5618(self):

        @njit
        def foo():
            values = np.zeros(1)
            tmp = 666
            if tmp:
                values[0] = tmp
            return values

        self.assertPreciseEqual(foo.py_func()[0], 666.)
        self.assertPreciseEqual(foo()[0], 666.)


class TestBranchPruneSSA(MemoryLeakMixin, TestCase):
    # Tests SSA rewiring of phi nodes after branch pruning.

    class SSAPrunerCompiler(CompilerBase):
        def define_pipelines(self):
            # This is a simple pipeline that does branch pruning on IR in SSA
            # form, then types and lowers as per the standard nopython pipeline.
            pm = PassManager("testing pm")
            pm.add_pass(TranslateByteCode, "analyzing bytecode")
            pm.add_pass(IRProcessing, "processing IR")
            # SSA early
            pm.add_pass(ReconstructSSA, "ssa")
            pm.add_pass(DeadBranchPrune, "dead branch pruning")
            # type and then lower as usual
            pm.add_pass(PreserveIR, "preserves the IR as metadata")
            dpb = DefaultPassBuilder
            typed_passes = dpb.define_typed_pipeline(self.state)
            pm.passes.extend(typed_passes.passes)
            lowering_passes = dpb.define_nopython_lowering_pipeline(self.state)
            pm.passes.extend(lowering_passes.passes)
            pm.finalize()
            return [pm]

    def test_ssa_update_phi(self):
        # This checks that dead branch pruning is rewiring phi nodes correctly
        # after a block containing an incoming for a phi is removed.

        @njit(pipeline_class=self.SSAPrunerCompiler)
        def impl(p=None, q=None):
            z = 1
            r = False
            if p is None:
                r = True # live

            if r and q is not None:
                z = 20 # dead

            # one of the incoming blocks for z is dead, the phi needs an update
            # were this not done, it would refer to variables that do not exist
            # and result in a lowering error.
            return z, r

        self.assertPreciseEqual(impl(), impl.py_func())

    def test_ssa_replace_phi(self):
        # This checks that when a phi only has one incoming, because the other
        # has been pruned, that a direct assignment is used instead.

        @njit(pipeline_class=self.SSAPrunerCompiler)
        def impl(p=None):
            z = 0
            if p is None:
                z = 10
            else:
                z = 20

            return z

        self.assertPreciseEqual(impl(), impl.py_func())
        func_ir = impl.overloads[impl.signatures[0]].metadata['preserved_ir']

        # check the func_ir, make sure there's no phi nodes
        for blk in func_ir.blocks.values():
            self.assertFalse([*blk.find_exprs('phi')])


class TestBranchPrunePostSemanticConstRewrites(TestBranchPruneBase):
    # Tests that semantic constants rewriting works by virtue of branch pruning

    def test_array_ndim_attr(self):

        def impl(array):
            if array.ndim == 2:
                if array.shape[1] == 2:
                    return 1
            else:
                return 10

        self.assert_prune(impl, (types.Array(types.float64, 2, 'C'),), [False,
                                                                        None],
                          np.zeros((2, 3)))
        self.assert_prune(impl, (types.Array(types.float64, 1, 'C'),), [True,
                                                                        'both'],
                          np.zeros((2,)))

    def test_tuple_len(self):

        def impl(tup):
            if len(tup) == 3:
                if tup[2] == 2:
                    return 1
            else:
                return 0

        self.assert_prune(impl, (types.UniTuple(types.int64, 3),), [False,
                                                                    None],
                          tuple([1, 2, 3]))
        self.assert_prune(impl, (types.UniTuple(types.int64, 2),), [True,
                                                                    'both'],
                          tuple([1, 2]))

    def test_attr_not_len(self):
        # The purpose of this test is to make sure that the conditions guarding
        # the rewrite part do not themselves raise exceptions.
        # This produces an `ir.Expr` call node for `float.as_integer_ratio`,
        # which is a getattr() on `float`.

        @njit
        def test():
            float.as_integer_ratio(1.23)

        # this should raise a TypingError
        with self.assertRaises(errors.TypingError) as e:
            test()

        self.assertIn("Unknown attribute 'as_integer_ratio'", str(e.exception))

    def test_ndim_not_on_array(self):

        FakeArray = collections.namedtuple('FakeArray', ['ndim'])
        fa = FakeArray(ndim=2)

        def impl(fa):
            if fa.ndim == 2:
                return fa.ndim
            else:
                object()

        # check prune works for array ndim
        self.assert_prune(impl, (types.Array(types.float64, 2, 'C'),), [False],
                          np.zeros((2, 3)))

        # check prune fails for something with `ndim` attr that is not array
        FakeArrayType = types.NamedUniTuple(types.int64, 1, FakeArray)
        self.assert_prune(impl, (FakeArrayType,), [None], fa,
                          flags={'nopython':False, 'forceobj':True})

    def test_semantic_const_propagates_before_static_rewrites(self):
        # see issue #5015, the ndim needs writing in as a const before
        # the rewrite passes run to make e.g. getitems static where possible
        @njit
        def impl(a, b):
            return a.shape[:b.ndim]

        args = (np.zeros((5, 4, 3, 2)), np.zeros((1, 1)))

        self.assertPreciseEqual(impl(*args), impl.py_func(*args))

    def test_tuple_const_propagation(self):
        @njit(pipeline_class=IRPreservingTestPipeline)
        def impl(*args):
            s = 0
            for arg in literal_unroll(args):
                s += len(arg)
            return s

        inp = ((), (1, 2, 3), ())
        self.assertPreciseEqual(impl(*inp), impl.py_func(*inp))

        ol = impl.overloads[impl.signatures[0]]
        func_ir = ol.metadata['preserved_ir']
        # make sure one of the inplace binop args is a Const
        binop_consts = set()
        for blk in func_ir.blocks.values():
            for expr in blk.find_exprs('inplace_binop'):
                inst = blk.find_variable_assignment(expr.rhs.name)
                self.assertIsInstance(inst.value, ir.Const)
                binop_consts.add(inst.value.value)
        self.assertEqual(binop_consts, {len(x) for x in inp})
