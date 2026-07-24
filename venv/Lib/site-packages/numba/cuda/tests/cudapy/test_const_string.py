import re
import numpy as np
from numba import cuda
from numba.cuda.testing import unittest, skip_on_cudasim, CUDATestCase
from llvmlite import ir


@skip_on_cudasim("This is testing CUDA backend code generation")
class TestConstStringCodegen(unittest.TestCase):
    def test_const_string(self):
        # These imports are incompatible with CUDASIM
        from numba.cuda.descriptor import cuda_target
        from numba.cuda.cudadrv.nvvm import compile_ir

        targetctx = cuda_target.target_context
        mod = targetctx.create_module("")
        textstring = 'A Little Brown Fox'
        gv0 = targetctx.insert_const_string(mod, textstring)
        # Insert the same const string a second time - the first should be
        # reused.
        targetctx.insert_const_string(mod, textstring)

        res = re.findall(r"@\"__conststring__.*internal.*constant.*\["
                         r"19\s+x\s+i8\]", str(mod))
        # Ensure that the const string was only inserted once
        self.assertEqual(len(res), 1)

        fnty = ir.FunctionType(ir.IntType(8).as_pointer(), [])

        # Using insert_const_string
        fn = ir.Function(mod, fnty, "test_insert_const_string")
        builder = ir.IRBuilder(fn.append_basic_block())
        res = builder.addrspacecast(gv0, ir.PointerType(ir.IntType(8)),
                                    'generic')
        builder.ret(res)

        matches = re.findall(r"@\"__conststring__.*internal.*constant.*\["
                             r"19\s+x\s+i8\]", str(mod))
        self.assertEqual(len(matches), 1)

        # Using insert_string_const_addrspace
        fn = ir.Function(mod, fnty, "test_insert_string_const_addrspace")
        builder = ir.IRBuilder(fn.append_basic_block())
        res = targetctx.insert_string_const_addrspace(builder, textstring)
        builder.ret(res)

        matches = re.findall(r"@\"__conststring__.*internal.*constant.*\["
                             r"19\s+x\s+i8\]", str(mod))
        self.assertEqual(len(matches), 1)

        ptx = compile_ir(str(mod)).decode('ascii')
        matches = list(re.findall(r"\.const.*__conststring__", ptx))

        self.assertEqual(len(matches), 1)


# Inspired by the reproducer from Issue #7041.
class TestConstString(CUDATestCase):
    def test_assign_const_unicode_string(self):
        @cuda.jit
        def str_assign(arr):
            i = cuda.grid(1)
            if i < len(arr):
                arr[i] = "XYZ"

        n_strings = 8
        arr = np.zeros(n_strings + 1, dtype="<U12")
        str_assign[1, n_strings](arr)

        # Expected result, e.g.:
        #     ['XYZ' 'XYZ' 'XYZ' 'XYZ' 'XYZ' 'XYZ' 'XYZ' 'XYZ' '']
        expected = np.zeros_like(arr)
        expected[:-1] = 'XYZ'
        expected[-1] = ''
        np.testing.assert_equal(arr, expected)

    def test_assign_const_byte_string(self):
        @cuda.jit
        def bytes_assign(arr):
            i = cuda.grid(1)
            if i < len(arr):
                arr[i] = b"XYZ"

        n_strings = 8
        arr = np.zeros(n_strings + 1, dtype="S12")
        bytes_assign[1, n_strings](arr)

        # Expected result, e.g.:
        #     [b'XYZ' b'XYZ' b'XYZ' b'XYZ' b'XYZ' b'XYZ' b'XYZ' b'XYZ' b'']
        expected = np.zeros_like(arr)
        expected[:-1] = b'XYZ'
        expected[-1] = b''
        np.testing.assert_equal(arr, expected)

    def test_assign_const_string_in_record(self):
        @cuda.jit
        def f(a):
            a[0]['x'] = 1
            a[0]['y'] = 'ABC'
            a[1]['x'] = 2
            a[1]['y'] = 'XYZ'

        dt = np.dtype([('x', np.int32), ('y', np.dtype('<U12'))])
        a = np.zeros(2, dt)

        f[1, 1](a)

        reference = np.asarray([(1, 'ABC'), (2, 'XYZ')], dtype=dt)
        np.testing.assert_array_equal(reference, a)

    def test_assign_const_bytes_in_record(self):
        @cuda.jit
        def f(a):
            a[0]['x'] = 1
            a[0]['y'] = b'ABC'
            a[1]['x'] = 2
            a[1]['y'] = b'XYZ'

        dt = np.dtype([('x', np.float32), ('y', np.dtype('S12'))])
        a = np.zeros(2, dt)

        f[1, 1](a)

        reference = np.asarray([(1, b'ABC'), (2, b'XYZ')], dtype=dt)
        np.testing.assert_array_equal(reference, a)


if __name__ == '__main__':
    unittest.main()
