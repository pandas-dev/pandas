from numba.cuda.testing import unittest, skip_on_cudasim
import operator
from numba.core import types, typing
from numba.cuda.cudadrv import nvvm


@unittest.skipIf(not nvvm.is_available(), "No libNVVM")
@skip_on_cudasim("Skip on simulator due to use of cuda_target")
class TestFunctionResolution(unittest.TestCase):
    def test_fp16_binary_operators(self):
        from numba.cuda.descriptor import cuda_target
        ops = (operator.add, operator.iadd, operator.sub, operator.isub,
               operator.mul, operator.imul)
        for op in ops:
            fp16 = types.float16
            typingctx = cuda_target.typing_context
            typingctx.refresh()
            fnty = typingctx.resolve_value_type(op)
            out = typingctx.resolve_function_type(fnty, (fp16, fp16), {})
            self.assertEqual(out, typing.signature(fp16, fp16, fp16),
                             msg=str(out))

    def test_fp16_unary_operators(self):
        from numba.cuda.descriptor import cuda_target
        ops = (operator.neg, abs)
        for op in ops:
            fp16 = types.float16
            typingctx = cuda_target.typing_context
            typingctx.refresh()
            fnty = typingctx.resolve_value_type(op)
            out = typingctx.resolve_function_type(fnty, (fp16,), {})
            self.assertEqual(out, typing.signature(fp16, fp16), msg=str(out))


if __name__ == '__main__':
    unittest.main()
