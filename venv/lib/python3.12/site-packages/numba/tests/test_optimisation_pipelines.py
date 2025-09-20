import unittest

from numba.tests.support import (TestCase, override_config)
from numba import njit
from numba.core import types
import llvmlite.binding as llvm


class TestPassManagerOptimization(TestCase):
    """ Tests that pass manager is not overriding the intended
    optimization level.
    """

    def _get_llvmir(self, fn, sig):
        with override_config('OPT', 0):
            fn.compile(sig)
            return fn.inspect_llvm(sig)

    def test_override_config(self):
        @njit(debug=True, error_model='numpy')
        def foo(a):
            b = a + 1.23
            c = b * 2.34
            d = b / c
            print(d)
            return d

        sig = (types.float64,)
        full_ir = self._get_llvmir(foo, sig=sig)

        module = llvm.parse_assembly(full_ir)

        name = foo.overloads[foo.signatures[0]].fndesc.mangled_name
        funcs = [x for x in module.functions if x.name == name]
        self.assertEqual(len(funcs), 1)
        func = funcs[0]
        blocks = [x for x in func.blocks]
        self.assertGreater(len(blocks), 1)
        block = blocks[0]

        # Find sequence with non-debug instructions
        instrs = [x for x in block.instructions if x.opcode != 'call']
        op_expect = {'fadd', 'fmul', 'fdiv'}
        started = False
        for x in instrs:
            if x.opcode in op_expect:
                op_expect.remove(x.opcode)
                if not started:
                    started = True
            elif op_expect and started:
                break

        self.assertGreater(len(op_expect), 0,
                           "Function was optimized unexpectedly")


if __name__ == '__main__':
    unittest.main()
