from numba.cuda.cudadrv import nvvm
from numba.cuda.testing import skip_on_cudasim
from numba.core import utils

from llvmlite import ir
from llvmlite import binding as llvm

import unittest


original = "call void @llvm.memset.p0i8.i64(" \
           "i8* align 4 %arg.x.41, i8 0, i64 %0, i1 false)"

missing_align = "call void @llvm.memset.p0i8.i64(" \
                "i8* %arg.x.41, i8 0, i64 %0, i1 false)"


@skip_on_cudasim('libNVVM not supported in simulator')
@unittest.skipIf(utils.MACHINE_BITS == 32, "CUDA not support for 32-bit")
@unittest.skipIf(not nvvm.is_available(), "No libNVVM")
class TestNvvmWithoutCuda(unittest.TestCase):
    def test_nvvm_accepts_encoding(self):
        # Test that NVVM will accept a constant containing all possible 8-bit
        # characters. Taken from the test case added in llvmlite PR #53:
        #
        #     https://github.com/numba/llvmlite/pull/53
        #
        # This test case is included in Numba to ensure that the encoding used
        # by llvmlite (e.g. utf-8, latin1, etc.) does not result in an input to
        # NVVM that it cannot parse correctly

        # Create a module with a constant containing all 8-bit characters
        c = ir.Constant(ir.ArrayType(ir.IntType(8), 256),
                        bytearray(range(256)))
        m = ir.Module()
        m.triple = 'nvptx64-nvidia-cuda'
        nvvm.add_ir_version(m)
        gv = ir.GlobalVariable(m, c.type, "myconstant")
        gv.global_constant = True
        gv.initializer = c
        m.data_layout = nvvm.NVVM().data_layout

        # Parse with LLVM then dump the parsed module into NVVM
        parsed = llvm.parse_assembly(str(m))
        ptx = nvvm.compile_ir(str(parsed))

        # Ensure all characters appear in the generated constant array.
        elements = ", ".join([str(i) for i in range(256)])
        myconstant = f"myconstant[256] = {{{elements}}}".encode('utf-8')
        self.assertIn(myconstant, ptx)


if __name__ == '__main__':
    unittest.main()
