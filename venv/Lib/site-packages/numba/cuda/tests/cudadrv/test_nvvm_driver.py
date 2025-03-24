import warnings

from llvmlite import ir
from numba.cuda.cudadrv import nvvm, runtime
from numba.cuda.testing import unittest
from numba.cuda.cudadrv.nvvm import LibDevice, NvvmError, NVVM
from numba.cuda.testing import skip_on_cudasim


@skip_on_cudasim('NVVM Driver unsupported in the simulator')
class TestNvvmDriver(unittest.TestCase):
    def get_nvvmir(self):
        versions = NVVM().get_ir_version()
        data_layout = NVVM().data_layout
        return nvvmir_generic.format(data_layout=data_layout, v=versions)

    def test_nvvm_compile_simple(self):
        nvvmir = self.get_nvvmir()
        ptx = nvvm.compile_ir(nvvmir).decode('utf8')
        self.assertTrue('simple' in ptx)
        self.assertTrue('ave' in ptx)

    def test_nvvm_compile_nullary_option(self):
        # Tests compilation with an option that doesn't take an argument
        # ("-gen-lto") - all other NVVM options are of the form
        # "-<name>=<value>"

        # -gen-lto is not available prior to CUDA 11.5
        if runtime.get_version() < (11, 5):
            self.skipTest("-gen-lto unavailable in this toolkit version")

        nvvmir = self.get_nvvmir()
        ltoir = nvvm.compile_ir(nvvmir, opt=3, gen_lto=None, arch="compute_52")

        # Verify we correctly passed the option by checking if we got LTOIR
        # from NVVM (by looking for the expected magic number for LTOIR)
        self.assertEqual(ltoir[:4], b'\xed\x43\x4e\x7f')

    def test_nvvm_bad_option(self):
        # Ensure that unsupported / non-existent options are reported as such
        # to the user / caller
        msg = "-made-up-option=2 is an unsupported option"
        with self.assertRaisesRegex(NvvmError, msg):
            nvvm.compile_ir("", made_up_option=2)

    def test_nvvm_from_llvm(self):
        m = ir.Module("test_nvvm_from_llvm")
        m.triple = 'nvptx64-nvidia-cuda'
        nvvm.add_ir_version(m)
        fty = ir.FunctionType(ir.VoidType(), [ir.IntType(32)])
        kernel = ir.Function(m, fty, name='mycudakernel')
        bldr = ir.IRBuilder(kernel.append_basic_block('entry'))
        bldr.ret_void()
        nvvm.set_cuda_kernel(kernel)

        m.data_layout = NVVM().data_layout
        ptx = nvvm.compile_ir(str(m)).decode('utf8')
        self.assertTrue('mycudakernel' in ptx)
        self.assertTrue('.address_size 64' in ptx)

    def test_used_list(self):
        # Construct a module
        m = ir.Module("test_used_list")
        m.triple = 'nvptx64-nvidia-cuda'
        m.data_layout = NVVM().data_layout
        nvvm.add_ir_version(m)

        # Add a function and mark it as a kernel
        fty = ir.FunctionType(ir.VoidType(), [ir.IntType(32)])
        kernel = ir.Function(m, fty, name='mycudakernel')
        bldr = ir.IRBuilder(kernel.append_basic_block('entry'))
        bldr.ret_void()
        nvvm.set_cuda_kernel(kernel)

        # Verify that the used list was correctly constructed
        used_lines = [line for line in str(m).splitlines()
                      if 'llvm.used' in line]
        msg = 'Expected exactly one @"llvm.used" array'
        self.assertEqual(len(used_lines), 1, msg)

        used_line = used_lines[0]
        # Kernel should be referenced in the used list
        self.assertIn("mycudakernel", used_line)
        # Check linkage of the used list
        self.assertIn("appending global", used_line)
        # Ensure used list is in the metadata section
        self.assertIn('section "llvm.metadata"', used_line)

    def test_nvvm_ir_verify_fail(self):
        m = ir.Module("test_bad_ir")
        m.triple = "unknown-unknown-unknown"
        m.data_layout = NVVM().data_layout
        nvvm.add_ir_version(m)
        with self.assertRaisesRegex(NvvmError, 'Invalid target triple'):
            nvvm.compile_ir(str(m))

    def _test_nvvm_support(self, arch):
        compute_xx = 'compute_{0}{1}'.format(*arch)
        nvvmir = self.get_nvvmir()
        ptx = nvvm.compile_ir(nvvmir, arch=compute_xx, ftz=1, prec_sqrt=0,
                              prec_div=0).decode('utf8')
        self.assertIn(".target sm_{0}{1}".format(*arch), ptx)
        self.assertIn('simple', ptx)
        self.assertIn('ave', ptx)

    def test_nvvm_support(self):
        """Test supported CC by NVVM
        """
        for arch in nvvm.get_supported_ccs():
            self._test_nvvm_support(arch=arch)

    def test_nvvm_warning(self):
        m = ir.Module("test_nvvm_warning")
        m.triple = 'nvptx64-nvidia-cuda'
        m.data_layout = NVVM().data_layout
        nvvm.add_ir_version(m)

        fty = ir.FunctionType(ir.VoidType(), [])
        kernel = ir.Function(m, fty, name='inlinekernel')
        builder = ir.IRBuilder(kernel.append_basic_block('entry'))
        builder.ret_void()
        nvvm.set_cuda_kernel(kernel)

        # Add the noinline attribute to trigger NVVM to generate a warning
        kernel.attributes.add('noinline')

        with warnings.catch_warnings(record=True) as w:
            nvvm.compile_ir(str(m))

        self.assertEqual(len(w), 1)
        self.assertIn('overriding noinline attribute', str(w[0]))


@skip_on_cudasim('NVVM Driver unsupported in the simulator')
class TestArchOption(unittest.TestCase):
    def test_get_arch_option(self):
        # Test returning the nearest lowest arch.
        self.assertEqual(nvvm.get_arch_option(5, 3), 'compute_53')
        self.assertEqual(nvvm.get_arch_option(7, 5), 'compute_75')
        self.assertEqual(nvvm.get_arch_option(7, 7), 'compute_75')
        # Test known arch.
        supported_cc = nvvm.get_supported_ccs()
        for arch in supported_cc:
            self.assertEqual(nvvm.get_arch_option(*arch), 'compute_%d%d' % arch)
        self.assertEqual(nvvm.get_arch_option(1000, 0),
                         'compute_%d%d' % supported_cc[-1])


@skip_on_cudasim('NVVM Driver unsupported in the simulator')
class TestLibDevice(unittest.TestCase):
    def test_libdevice_load(self):
        # Test that constructing LibDevice gives a bitcode file
        libdevice = LibDevice()
        self.assertEqual(libdevice.bc[:4], b'BC\xc0\xde')


nvvmir_generic = '''\
target triple="nvptx64-nvidia-cuda"
target datalayout = "{data_layout}"

define i32 @ave(i32 %a, i32 %b) {{
entry:
%add = add nsw i32 %a, %b
%div = sdiv i32 %add, 2
ret i32 %div
}}

define void @simple(i32* %data) {{
entry:
%0 = call i32 @llvm.nvvm.read.ptx.sreg.ctaid.x()
%1 = call i32 @llvm.nvvm.read.ptx.sreg.ntid.x()
%mul = mul i32 %0, %1
%2 = call i32 @llvm.nvvm.read.ptx.sreg.tid.x()
%add = add i32 %mul, %2
%call = call i32 @ave(i32 %add, i32 %add)
%idxprom = sext i32 %add to i64
%arrayidx = getelementptr inbounds i32, i32* %data, i64 %idxprom
store i32 %call, i32* %arrayidx, align 4
ret void
}}

declare i32 @llvm.nvvm.read.ptx.sreg.ctaid.x() nounwind readnone

declare i32 @llvm.nvvm.read.ptx.sreg.ntid.x() nounwind readnone

declare i32 @llvm.nvvm.read.ptx.sreg.tid.x() nounwind readnone

!nvvmir.version = !{{!1}}
!1 = !{{i32 {v[0]}, i32 {v[1]}, i32 {v[2]}, i32 {v[3]}}}

!nvvm.annotations = !{{!2}}
!2 = !{{void (i32*)* @simple, !"kernel", i32 1}}

@"llvm.used" = appending global [1 x i8*] [i8* bitcast (void (i32*)* @simple to i8*)], section "llvm.metadata"
'''  # noqa: E501


if __name__ == '__main__':
    unittest.main()
