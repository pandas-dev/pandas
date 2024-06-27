import ctypes
import numpy as np
import weakref

from numba import cuda
from numba.core import config
from numba.cuda.testing import unittest, CUDATestCase, skip_on_cudasim
from numba.tests.support import linux_only

if not config.ENABLE_CUDASIM:
    class DeviceOnlyEMMPlugin(cuda.HostOnlyCUDAMemoryManager):
        """
        Dummy EMM Plugin implementation for testing. It memorises which plugin
        API methods have been called so that the tests can check that Numba
        called into the plugin as expected.
        """

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            # For tracking our dummy allocations
            self.allocations = {}
            self.count = 0

            # For tracking which methods have been called
            self.initialized = False
            self.memalloc_called = False
            self.reset_called = False
            self.get_memory_info_called = False
            self.get_ipc_handle_called = False

        def memalloc(self, size):
            # We maintain a list of allocations and keep track of them, so that
            # we can test that the finalizers of objects returned by memalloc
            # get called.

            # Numba should have initialized the memory manager when preparing
            # the context for use, prior to any memalloc call.
            if not self.initialized:
                raise RuntimeError("memalloc called before initialize")
            self.memalloc_called = True

            # Create an allocation and record it
            self.count += 1
            alloc_count = self.count
            self.allocations[alloc_count] = size

            # The finalizer deletes the record from our internal dict of
            # allocations.
            finalizer_allocs = self.allocations

            def finalizer():
                del finalizer_allocs[alloc_count]

            # We use an AutoFreePointer so that the finalizer will be run when
            # the reference count drops to zero.
            ctx = weakref.proxy(self.context)
            ptr = ctypes.c_void_p(alloc_count)
            return cuda.cudadrv.driver.AutoFreePointer(ctx, ptr, size,
                                                       finalizer=finalizer)

        def initialize(self):
            # No special initialization needed.
            self.initialized = True

        def reset(self):
            # We remove all allocations on reset, just as a real EMM Plugin
            # would do. Note that our finalizers in memalloc don't check
            # whether the allocations are still alive, so running them after
            # reset will detect any allocations that are floating around at
            # exit time; however, the atexit finalizer for weakref will only
            # print a traceback, not terminate the interpreter abnormally.
            self.reset_called = True

        def get_memory_info(self):
            # Return some dummy memory information
            self.get_memory_info_called = True
            return cuda.MemoryInfo(free=32, total=64)

        def get_ipc_handle(self, memory):
            # The dummy IPC handle is only a string, so it is important that
            # the tests don't try to do too much with it (e.g. open / close
            # it).
            self.get_ipc_handle_called = True
            return "Dummy IPC handle for alloc %s" % memory.device_pointer.value

        @property
        def interface_version(self):
            # The expected version for an EMM Plugin.
            return 1

    class BadVersionEMMPlugin(DeviceOnlyEMMPlugin):
        """A plugin that claims to implement a different interface version"""

        @property
        def interface_version(self):
            return 2


@skip_on_cudasim('EMM Plugins not supported on CUDA simulator')
class TestDeviceOnlyEMMPlugin(CUDATestCase):
    """
    Tests that the API of an EMM Plugin that implements device allocations
    only is used correctly by Numba.
    """

    def setUp(self):
        super().setUp()
        # Always start afresh with a new context and memory manager
        cuda.close()
        cuda.set_memory_manager(DeviceOnlyEMMPlugin)

    def tearDown(self):
        super().tearDown()
        # Unset the memory manager for subsequent tests
        cuda.close()
        cuda.cudadrv.driver._memory_manager = None

    def test_memalloc(self):
        mgr = cuda.current_context().memory_manager

        # Allocate an array and check that memalloc was called with the correct
        # size.
        arr_1 = np.arange(10)
        d_arr_1 = cuda.device_array_like(arr_1)
        self.assertTrue(mgr.memalloc_called)
        self.assertEqual(mgr.count, 1)
        self.assertEqual(mgr.allocations[1], arr_1.nbytes)

        # Allocate again, with a different size, and check that it is also
        # correct.
        arr_2 = np.arange(5)
        d_arr_2 = cuda.device_array_like(arr_2)
        self.assertEqual(mgr.count, 2)
        self.assertEqual(mgr.allocations[2], arr_2.nbytes)

        # Remove the first array, and check that our finalizer was called for
        # the first array only.
        del d_arr_1
        self.assertNotIn(1, mgr.allocations)
        self.assertIn(2, mgr.allocations)

        # Remove the second array and check that its finalizer was also
        # called.
        del d_arr_2
        self.assertNotIn(2, mgr.allocations)

    def test_initialized_in_context(self):
        # If we have a CUDA context, it should already have initialized its
        # memory manager.
        self.assertTrue(cuda.current_context().memory_manager.initialized)

    def test_reset(self):
        ctx = cuda.current_context()
        ctx.reset()
        self.assertTrue(ctx.memory_manager.reset_called)

    def test_get_memory_info(self):
        ctx = cuda.current_context()
        meminfo = ctx.get_memory_info()
        self.assertTrue(ctx.memory_manager.get_memory_info_called)
        self.assertEqual(meminfo.free, 32)
        self.assertEqual(meminfo.total, 64)

    @linux_only
    def test_get_ipc_handle(self):
        # We don't attempt to close the IPC handle in this test because Numba
        # will be expecting a real IpcHandle object to have been returned from
        # get_ipc_handle, and it would cause problems to do so.
        arr = np.arange(2)
        d_arr = cuda.device_array_like(arr)
        ipch = d_arr.get_ipc_handle()
        ctx = cuda.current_context()
        self.assertTrue(ctx.memory_manager.get_ipc_handle_called)
        self.assertIn("Dummy IPC handle for alloc 1", ipch._ipc_handle)


@skip_on_cudasim('EMM Plugins not supported on CUDA simulator')
class TestBadEMMPluginVersion(CUDATestCase):
    """
    Ensure that Numba rejects EMM Plugins with incompatible version
    numbers.
    """

    def test_bad_plugin_version(self):
        with self.assertRaises(RuntimeError) as raises:
            cuda.set_memory_manager(BadVersionEMMPlugin)
        self.assertIn('version 1 required', str(raises.exception))


if __name__ == '__main__':
    unittest.main()
