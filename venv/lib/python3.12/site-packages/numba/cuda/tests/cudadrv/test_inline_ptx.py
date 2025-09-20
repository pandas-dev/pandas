from llvmlite import ir

from numba.cuda.cudadrv import nvvm
from numba.cuda.testing import unittest, ContextResettingTestCase
from numba.cuda.testing import skip_on_cudasim


@skip_on_cudasim('Inline PTX cannot be used in the simulator')
class TestCudaInlineAsm(ContextResettingTestCase):
    def test_inline_rsqrt(self):
        mod = ir.Module(__name__)
        mod.triple = 'nvptx64-nvidia-cuda'
        nvvm.add_ir_version(mod)
        fnty = ir.FunctionType(ir.VoidType(), [ir.PointerType(ir.FloatType())])
        fn = ir.Function(mod, fnty, 'cu_rsqrt')
        bldr = ir.IRBuilder(fn.append_basic_block('entry'))

        rsqrt_approx_fnty = ir.FunctionType(ir.FloatType(), [ir.FloatType()])
        inlineasm = ir.InlineAsm(rsqrt_approx_fnty,
                                 'rsqrt.approx.f32 $0, $1;',
                                 '=f,f', side_effect=True)
        val = bldr.load(fn.args[0])
        res = bldr.call(inlineasm, [val])

        bldr.store(res, fn.args[0])
        bldr.ret_void()

        # generate ptx
        mod.data_layout = nvvm.NVVM().data_layout
        nvvm.set_cuda_kernel(fn)
        nvvmir = str(mod)
        ptx = nvvm.compile_ir(nvvmir)
        self.assertTrue('rsqrt.approx.f32' in str(ptx))


if __name__ == '__main__':
    unittest.main()
