from ctypes import byref, c_int, c_void_p, sizeof

from numba.cuda.cudadrv.driver import (host_to_device, device_to_host, driver,
                                       launch_kernel)
from numba.cuda.cudadrv import devices, drvapi, driver as _driver
from numba.cuda.testing import unittest, CUDATestCase
from numba.cuda.testing import skip_on_cudasim


ptx1 = '''
    .version 1.4
    .target sm_10, map_f64_to_f32

    .entry _Z10helloworldPi (
    .param .u64 __cudaparm__Z10helloworldPi_A)
    {
    .reg .u32 %r<3>;
    .reg .u64 %rd<6>;
    .loc	14	4	0
$LDWbegin__Z10helloworldPi:
    .loc	14	6	0
    cvt.s32.u16 	%r1, %tid.x;
    ld.param.u64 	%rd1, [__cudaparm__Z10helloworldPi_A];
    cvt.u64.u16 	%rd2, %tid.x;
    mul.lo.u64 	%rd3, %rd2, 4;
    add.u64 	%rd4, %rd1, %rd3;
    st.global.s32 	[%rd4+0], %r1;
    .loc	14	7	0
    exit;
$LDWend__Z10helloworldPi:
    } // _Z10helloworldPi
'''

ptx2 = '''
.version 3.0
.target sm_20
.address_size 64

    .file	1 "/tmp/tmpxft_000012c7_00000000-9_testcuda.cpp3.i"
    .file	2 "testcuda.cu"

.entry _Z10helloworldPi(
    .param .u64 _Z10helloworldPi_param_0
)
{
    .reg .s32 	%r<3>;
    .reg .s64 	%rl<5>;


    ld.param.u64 	%rl1, [_Z10helloworldPi_param_0];
    cvta.to.global.u64 	%rl2, %rl1;
    .loc 2 6 1
    mov.u32 	%r1, %tid.x;
    mul.wide.u32 	%rl3, %r1, 4;
    add.s64 	%rl4, %rl2, %rl3;
    st.global.u32 	[%rl4], %r1;
    .loc 2 7 2
    ret;
}
'''


@skip_on_cudasim('CUDA Driver API unsupported in the simulator')
class TestCudaDriver(CUDATestCase):
    def setUp(self):
        super().setUp()
        self.assertTrue(len(devices.gpus) > 0)
        self.context = devices.get_context()
        device = self.context.device
        ccmajor, _ = device.compute_capability
        if ccmajor >= 2:
            self.ptx = ptx2
        else:
            self.ptx = ptx1

    def tearDown(self):
        super().tearDown()
        del self.context

    def test_cuda_driver_basic(self):
        module = self.context.create_module_ptx(self.ptx)
        function = module.get_function('_Z10helloworldPi')

        array = (c_int * 100)()

        memory = self.context.memalloc(sizeof(array))
        host_to_device(memory, array, sizeof(array))

        ptr = memory.device_ctypes_pointer
        stream = 0

        if _driver.USE_NV_BINDING:
            ptr = c_void_p(int(ptr))
            stream = _driver.binding.CUstream(stream)

        launch_kernel(function.handle,  # Kernel
                      1,   1, 1,        # gx, gy, gz
                      100, 1, 1,        # bx, by, bz
                      0,                # dynamic shared mem
                      stream,           # stream
                      [ptr])            # arguments

        device_to_host(array, memory, sizeof(array))
        for i, v in enumerate(array):
            self.assertEqual(i, v)

        module.unload()

    def test_cuda_driver_stream_operations(self):
        module = self.context.create_module_ptx(self.ptx)
        function = module.get_function('_Z10helloworldPi')

        array = (c_int * 100)()

        stream = self.context.create_stream()

        with stream.auto_synchronize():
            memory = self.context.memalloc(sizeof(array))
            host_to_device(memory, array, sizeof(array), stream=stream)

            ptr = memory.device_ctypes_pointer
            if _driver.USE_NV_BINDING:
                ptr = c_void_p(int(ptr))

            launch_kernel(function.handle,  # Kernel
                          1,   1, 1,        # gx, gy, gz
                          100, 1, 1,        # bx, by, bz
                          0,                # dynamic shared mem
                          stream.handle,    # stream
                          [ptr])            # arguments

        device_to_host(array, memory, sizeof(array), stream=stream)

        for i, v in enumerate(array):
            self.assertEqual(i, v)

    def test_cuda_driver_default_stream(self):
        # Test properties of the default stream
        ds = self.context.get_default_stream()
        self.assertIn("Default CUDA stream", repr(ds))
        self.assertEqual(0, int(ds))
        # bool(stream) is the check that is done in memcpy to decide if async
        # version should be used. So the default (0) stream should be true-ish
        # even though 0 is usually false-ish in Python.
        self.assertTrue(ds)
        self.assertFalse(ds.external)

    def test_cuda_driver_legacy_default_stream(self):
        # Test properties of the legacy default stream
        ds = self.context.get_legacy_default_stream()
        self.assertIn("Legacy default CUDA stream", repr(ds))
        self.assertEqual(1, int(ds))
        self.assertTrue(ds)
        self.assertFalse(ds.external)

    def test_cuda_driver_per_thread_default_stream(self):
        # Test properties of the per-thread default stream
        ds = self.context.get_per_thread_default_stream()
        self.assertIn("Per-thread default CUDA stream", repr(ds))
        self.assertEqual(2, int(ds))
        self.assertTrue(ds)
        self.assertFalse(ds.external)

    def test_cuda_driver_stream(self):
        # Test properties of non-default streams
        s = self.context.create_stream()
        self.assertIn("CUDA stream", repr(s))
        self.assertNotIn("Default", repr(s))
        self.assertNotIn("External", repr(s))
        self.assertNotEqual(0, int(s))
        self.assertTrue(s)
        self.assertFalse(s.external)

    def test_cuda_driver_external_stream(self):
        # Test properties of a stream created from an external stream object.
        # We use the driver API directly to create a stream, to emulate an
        # external library creating a stream
        if _driver.USE_NV_BINDING:
            handle = driver.cuStreamCreate(0)
            ptr = int(handle)
        else:
            handle = drvapi.cu_stream()
            driver.cuStreamCreate(byref(handle), 0)
            ptr = handle.value
        s = self.context.create_external_stream(ptr)

        self.assertIn("External CUDA stream", repr(s))
        # Ensure neither "Default" nor "default"
        self.assertNotIn("efault", repr(s))
        self.assertEqual(ptr, int(s))
        self.assertTrue(s)
        self.assertTrue(s.external)

    def test_cuda_driver_occupancy(self):
        module = self.context.create_module_ptx(self.ptx)
        function = module.get_function('_Z10helloworldPi')

        value = self.context.get_active_blocks_per_multiprocessor(function,
                                                                  128, 128)
        self.assertTrue(value > 0)

        def b2d(bs):
            return bs

        grid, block = self.context.get_max_potential_block_size(function, b2d,
                                                                128, 128)
        self.assertTrue(grid > 0)
        self.assertTrue(block > 0)


class TestDevice(CUDATestCase):
    def test_device_get_uuid(self):
        # A device UUID looks like:
        #
        #     GPU-e6489c45-5b68-3b03-bab7-0e7c8e809643
        #
        # To test, we construct an RE that matches this form and verify that
        # the returned UUID matches.
        #
        # Device UUIDs may not conform to parts of the UUID specification (RFC
        # 4122) pertaining to versions and variants, so we do not extract and
        # validate the values of these bits.

        h = '[0-9a-f]{%d}'
        h4 = h % 4
        h8 = h % 8
        h12 = h % 12
        uuid_format = f'^GPU-{h8}-{h4}-{h4}-{h4}-{h12}$'

        dev = devices.get_context().device
        self.assertRegex(dev.uuid, uuid_format)


if __name__ == '__main__':
    unittest.main()
