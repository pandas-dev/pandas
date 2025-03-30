import numba
from numba.tests.support import TestCase, unittest
from numba.core.registry import cpu_target
from numba.core.compiler import CompilerBase, Flags
from numba.core.compiler_machinery import PassManager
from numba.core import types, ir, bytecode, compiler, ir_utils, registry
from numba.core.untyped_passes import (ExtractByteCode, TranslateByteCode,
                                       FixupArgs, IRProcessing,)

from numba.core.typed_passes import (NopythonTypeInference,
                                     type_inference_stage, DeadCodeElimination)
from numba.experimental import jitclass

# global constant for testing find_const
GLOBAL_B = 11


@jitclass([('val', numba.core.types.List(numba.intp))])
class Dummy(object):
    def __init__(self, val):
        self.val = val


class TestIrUtils(TestCase):
    """
    Tests ir handling utility functions like find_callname.
    """

    def test_obj_func_match(self):
        """Test matching of an object method (other than Array see #3449)
        """

        def test_func():
            d = Dummy([1])
            d.val.append(2)

        test_ir = compiler.run_frontend(test_func)
        typingctx = cpu_target.typing_context
        targetctx = cpu_target.target_context
        typing_res = type_inference_stage(
            typingctx, targetctx, test_ir, (), None)
        matched_call = ir_utils.find_callname(
            test_ir, test_ir.blocks[0].body[7].value, typing_res.typemap)
        self.assertTrue(isinstance(matched_call, tuple) and
                        len(matched_call) == 2 and
                        matched_call[0] == 'append')

    def test_dead_code_elimination(self):

        class Tester(CompilerBase):

            @classmethod
            def mk_pipeline(cls, args, return_type=None, flags=None, locals={},
                            library=None, typing_context=None,
                            target_context=None):
                if not flags:
                    flags = Flags()
                flags.nrt = True
                if typing_context is None:
                    typing_context = registry.cpu_target.typing_context
                if target_context is None:
                    target_context = registry.cpu_target.target_context
                return cls(typing_context, target_context, library, args,
                           return_type, flags, locals)

            def compile_to_ir(self, func, DCE=False):
                """
                Compile and return IR
                """
                func_id = bytecode.FunctionIdentity.from_function(func)
                self.state.func_id = func_id
                ExtractByteCode().run_pass(self.state)
                state = self.state

                name = "DCE_testing"
                pm = PassManager(name)
                pm.add_pass(TranslateByteCode, "analyzing bytecode")
                pm.add_pass(FixupArgs, "fix up args")
                pm.add_pass(IRProcessing, "processing IR")
                pm.add_pass(NopythonTypeInference, "nopython frontend")
                if DCE is True:
                    pm.add_pass(DeadCodeElimination, "DCE after typing")
                pm.finalize()
                pm.run(state)
                return state.func_ir

        def check_initial_ir(the_ir):
            # dead stuff:
            # a const int value 0xdead
            # an assign of above into to variable `dead`
            # a const int above 0xdeaddead
            # an assign of said int to variable `deaddead`
            # this is 2 statements to remove

            self.assertEqual(len(the_ir.blocks), 1)
            block = the_ir.blocks[0]
            deads = []
            for x in block.find_insts(ir.Assign):
                if isinstance(getattr(x, 'target', None), ir.Var):
                    if 'dead' in getattr(x.target, 'name', ''):
                        deads.append(x)

            self.assertEqual(len(deads), 2)
            for d in deads:
                # check the ir.Const is the definition and the value is expected
                const_val = the_ir.get_definition(d.value)
                self.assertTrue(int('0x%s' % d.target.name, 16),
                                const_val.value)

            return deads

        def check_dce_ir(the_ir):
            self.assertEqual(len(the_ir.blocks), 1)
            block = the_ir.blocks[0]
            deads = []
            consts = []
            for x in block.find_insts(ir.Assign):
                if isinstance(getattr(x, 'target', None), ir.Var):
                    if 'dead' in getattr(x.target, 'name', ''):
                        deads.append(x)
                if isinstance(getattr(x, 'value', None), ir.Const):
                    consts.append(x)
            self.assertEqual(len(deads), 0)

            # check the consts to make sure there's no reference to 0xdead or
            # 0xdeaddead
            for x in consts:
                self.assertTrue(x.value.value not in [0xdead, 0xdeaddead])

        def foo(x):
            y = x + 1
            dead = 0xdead  # noqa
            z = y + 2
            deaddead = 0xdeaddead  # noqa
            ret = z * z
            return ret

        test_pipeline = Tester.mk_pipeline((types.intp,))
        no_dce = test_pipeline.compile_to_ir(foo)
        removed = check_initial_ir(no_dce)

        test_pipeline = Tester.mk_pipeline((types.intp,))
        w_dce = test_pipeline.compile_to_ir(foo, DCE=True)
        check_dce_ir(w_dce)

        # check that the count of initial - removed = dce
        self.assertEqual(len(no_dce.blocks[0].body) - len(removed),
                         len(w_dce.blocks[0].body))

    def test_find_const_global(self):
        """
        Test find_const() for values in globals (ir.Global) and freevars
        (ir.FreeVar) that are considered constants for compilation.
        """
        FREEVAR_C = 12

        def foo(a):
            b = GLOBAL_B
            c = FREEVAR_C
            return a + b + c

        f_ir = compiler.run_frontend(foo)
        block = f_ir.blocks[0]
        const_b = None
        const_c = None

        for inst in block.body:
            if isinstance(inst, ir.Assign) and inst.target.name == 'b':
                const_b = ir_utils.guard(
                    ir_utils.find_const, f_ir, inst.target)
            if isinstance(inst, ir.Assign) and inst.target.name == 'c':
                const_c = ir_utils.guard(
                    ir_utils.find_const, f_ir, inst.target)

        self.assertEqual(const_b, GLOBAL_B)
        self.assertEqual(const_c, FREEVAR_C)

    def test_flatten_labels(self):
        """ tests flatten_labels """
        def foo(a):
            acc = 0
            if a > 3:
                acc += 1
                if a > 19:
                    return 53
            elif a < 1000:
                if a >= 12:
                    acc += 1
                for x in range(10):
                    acc -= 1
                    if acc < 2:
                        break
                else:
                    acc += 7
            else:
                raise ValueError("some string")
            # prevents inline of return on py310
            py310_defeat1 = 1  # noqa
            py310_defeat2 = 2  # noqa
            py310_defeat3 = 3  # noqa
            py310_defeat4 = 4  # noqa
            return acc

        def bar(a):
            acc = 0
            z = 12
            if a > 3:
                acc += 1
                z += 12
                if a > 19:
                    z += 12
                    return 53
            elif a < 1000:
                if a >= 12:
                    z += 12
                    acc += 1
                for x in range(10):
                    z += 12
                    acc -= 1
                    if acc < 2:
                        break
                else:
                    z += 12
                    acc += 7
            else:
                raise ValueError("some string")
            py310_defeat1 = 1  # noqa
            py310_defeat2 = 2  # noqa
            py310_defeat3 = 3  # noqa
            py310_defeat4 = 4  # noqa
            return acc

        def baz(a):
            acc = 0
            if a > 3:
                acc += 1
                if a > 19:
                    return 53
                else: # extra control flow in comparison to foo
                    return 55
            elif a < 1000:
                if a >= 12:
                    acc += 1
                for x in range(10):
                    acc -= 1
                    if acc < 2:
                        break
                else:
                    acc += 7
            else:
                raise ValueError("some string")
            py310_defeat1 = 1  # noqa
            py310_defeat2 = 2  # noqa
            py310_defeat3 = 3  # noqa
            py310_defeat4 = 4  # noqa
            return acc

        def get_flat_cfg(func):
            func_ir = ir_utils.compile_to_numba_ir(func, dict())
            flat_blocks = ir_utils.flatten_labels(func_ir.blocks)
            self.assertEqual(max(flat_blocks.keys()) + 1, len(func_ir.blocks))
            return ir_utils.compute_cfg_from_blocks(flat_blocks)

        foo_cfg = get_flat_cfg(foo)
        bar_cfg = get_flat_cfg(bar)
        baz_cfg = get_flat_cfg(baz)

        self.assertEqual(foo_cfg, bar_cfg)
        self.assertNotEqual(foo_cfg, baz_cfg)


if __name__ == "__main__":
    unittest.main()
