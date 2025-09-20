import multiprocessing as mp
import itertools
import traceback
import pickle

import numpy as np

from numba import cuda
from numba.cuda.cudadrv import driver
from numba.cuda.testing import (skip_on_arm, skip_on_cudasim,
                                skip_under_cuda_memcheck,
                                ContextResettingTestCase, ForeignArray)
from numba.tests.support import linux_only, windows_only
import unittest


def core_ipc_handle_test(the_work, result_queue):
    try:
        arr = the_work()
    # Catch anything going wrong in the worker function
    except:  # noqa: E722
        # FAILED. propagate the exception as a string
        succ = False
        out = traceback.format_exc()
    else:
        # OK. send the ndarray back
        succ = True
        out = arr
    result_queue.put((succ, out))


def base_ipc_handle_test(handle, size, result_queue):
    def the_work():
        dtype = np.dtype(np.intp)
        with cuda.open_ipc_array(handle, shape=size // dtype.itemsize,
                                 dtype=dtype) as darr:
            # copy the data to host
            return darr.copy_to_host()

    core_ipc_handle_test(the_work, result_queue)


def serialize_ipc_handle_test(handle, result_queue):
    def the_work():
        dtype = np.dtype(np.intp)
        darr = handle.open_array(cuda.current_context(),
                                 shape=handle.size // dtype.itemsize,
                                 dtype=dtype)
        # copy the data to host
        arr = darr.copy_to_host()
        handle.close()
        return arr

    core_ipc_handle_test(the_work, result_queue)


def ipc_array_test(ipcarr, result_queue):
    try:
        with ipcarr as darr:
            arr = darr.copy_to_host()
            try:
                # should fail to reopen
                with ipcarr:
                    pass
            except ValueError as e:
                if str(e) != 'IpcHandle is already opened':
                    raise AssertionError('invalid exception message')
            else:
                raise AssertionError('did not raise on reopen')
    # Catch any exception so we can propagate it
    except:  # noqa: E722
        # FAILED. propagate the exception as a string
        succ = False
        out = traceback.format_exc()
    else:
        # OK. send the ndarray back
        succ = True
        out = arr
    result_queue.put((succ, out))


@linux_only
@skip_under_cuda_memcheck('Hangs cuda-memcheck')
@skip_on_cudasim('Ipc not available in CUDASIM')
@skip_on_arm('CUDA IPC not supported on ARM in Numba')
class TestIpcMemory(ContextResettingTestCase):

    def test_ipc_handle(self):
        # prepare data for IPC
        arr = np.arange(10, dtype=np.intp)
        devarr = cuda.to_device(arr)

        # create IPC handle
        ctx = cuda.current_context()
        ipch = ctx.get_ipc_handle(devarr.gpu_data)

        # manually prepare for serialization as bytes
        if driver.USE_NV_BINDING:
            handle_bytes = ipch.handle.reserved
        else:
            handle_bytes = bytes(ipch.handle)
        size = ipch.size

        # spawn new process for testing
        ctx = mp.get_context('spawn')
        result_queue = ctx.Queue()
        args = (handle_bytes, size, result_queue)
        proc = ctx.Process(target=base_ipc_handle_test, args=args)
        proc.start()
        succ, out = result_queue.get()
        if not succ:
            self.fail(out)
        else:
            np.testing.assert_equal(arr, out)
        proc.join(3)

    def variants(self):
        # Test with no slicing and various different slices
        indices = (None, slice(3, None), slice(3, 8), slice(None, 8))
        # Test with a Numba DeviceNDArray, or an array from elsewhere through
        # the CUDA Array Interface
        foreigns = (False, True)
        return itertools.product(indices, foreigns)

    def check_ipc_handle_serialization(self, index_arg=None, foreign=False):
        # prepare data for IPC
        arr = np.arange(10, dtype=np.intp)
        devarr = cuda.to_device(arr)
        if index_arg is not None:
            devarr = devarr[index_arg]
        if foreign:
            devarr = cuda.as_cuda_array(ForeignArray(devarr))
        expect = devarr.copy_to_host()

        # create IPC handle
        ctx = cuda.current_context()
        ipch = ctx.get_ipc_handle(devarr.gpu_data)

        # pickle
        buf = pickle.dumps(ipch)
        ipch_recon = pickle.loads(buf)
        self.assertIs(ipch_recon.base, None)
        self.assertEqual(ipch_recon.size, ipch.size)

        if driver.USE_NV_BINDING:
            self.assertEqual(ipch_recon.handle.reserved, ipch.handle.reserved)
        else:
            self.assertEqual(tuple(ipch_recon.handle), tuple(ipch.handle))

        # spawn new process for testing
        ctx = mp.get_context('spawn')
        result_queue = ctx.Queue()
        args = (ipch, result_queue)
        proc = ctx.Process(target=serialize_ipc_handle_test, args=args)
        proc.start()
        succ, out = result_queue.get()
        if not succ:
            self.fail(out)
        else:
            np.testing.assert_equal(expect, out)
        proc.join(3)

    def test_ipc_handle_serialization(self):
        for index, foreign, in self.variants():
            with self.subTest(index=index, foreign=foreign):
                self.check_ipc_handle_serialization(index, foreign)

    def check_ipc_array(self, index_arg=None, foreign=False):
        # prepare data for IPC
        arr = np.arange(10, dtype=np.intp)
        devarr = cuda.to_device(arr)
        # Slice
        if index_arg is not None:
            devarr = devarr[index_arg]
        if foreign:
            devarr = cuda.as_cuda_array(ForeignArray(devarr))
        expect = devarr.copy_to_host()
        ipch = devarr.get_ipc_handle()

        # spawn new process for testing
        ctx = mp.get_context('spawn')
        result_queue = ctx.Queue()
        args = (ipch, result_queue)
        proc = ctx.Process(target=ipc_array_test, args=args)
        proc.start()
        succ, out = result_queue.get()
        if not succ:
            self.fail(out)
        else:
            np.testing.assert_equal(expect, out)
        proc.join(3)

    def test_ipc_array(self):
        for index, foreign, in self.variants():
            with self.subTest(index=index, foreign=foreign):
                self.check_ipc_array(index, foreign)


def staged_ipc_handle_test(handle, device_num, result_queue):
    def the_work():
        with cuda.gpus[device_num]:
            this_ctx = cuda.devices.get_context()
            deviceptr = handle.open_staged(this_ctx)
            arrsize = handle.size // np.dtype(np.intp).itemsize
            hostarray = np.zeros(arrsize, dtype=np.intp)
            cuda.driver.device_to_host(
                hostarray, deviceptr, size=handle.size,
            )
            handle.close()
        return hostarray

    core_ipc_handle_test(the_work, result_queue)


def staged_ipc_array_test(ipcarr, device_num, result_queue):
    try:
        with cuda.gpus[device_num]:
            with ipcarr as darr:
                arr = darr.copy_to_host()
                try:
                    # should fail to reopen
                    with ipcarr:
                        pass
                except ValueError as e:
                    if str(e) != 'IpcHandle is already opened':
                        raise AssertionError('invalid exception message')
                else:
                    raise AssertionError('did not raise on reopen')
    # Catch any exception so we can propagate it
    except:  # noqa: E722
        # FAILED. propagate the exception as a string
        succ = False
        out = traceback.format_exc()
    else:
        # OK. send the ndarray back
        succ = True
        out = arr
    result_queue.put((succ, out))


@linux_only
@skip_under_cuda_memcheck('Hangs cuda-memcheck')
@skip_on_cudasim('Ipc not available in CUDASIM')
@skip_on_arm('CUDA IPC not supported on ARM in Numba')
class TestIpcStaged(ContextResettingTestCase):
    def test_staged(self):
        # prepare data for IPC
        arr = np.arange(10, dtype=np.intp)
        devarr = cuda.to_device(arr)

        # spawn new process for testing
        mpctx = mp.get_context('spawn')
        result_queue = mpctx.Queue()

        # create IPC handle
        ctx = cuda.current_context()
        ipch = ctx.get_ipc_handle(devarr.gpu_data)
        # pickle
        buf = pickle.dumps(ipch)
        ipch_recon = pickle.loads(buf)
        self.assertIs(ipch_recon.base, None)
        if driver.USE_NV_BINDING:
            self.assertEqual(ipch_recon.handle.reserved, ipch.handle.reserved)
        else:
            self.assertEqual(tuple(ipch_recon.handle), tuple(ipch.handle))
        self.assertEqual(ipch_recon.size, ipch.size)

        # Test on every CUDA devices
        for device_num in range(len(cuda.gpus)):
            args = (ipch, device_num, result_queue)
            proc = mpctx.Process(target=staged_ipc_handle_test, args=args)
            proc.start()
            succ, out = result_queue.get()
            proc.join(3)
            if not succ:
                self.fail(out)
            else:
                np.testing.assert_equal(arr, out)

    def test_ipc_array(self):
        for device_num in range(len(cuda.gpus)):
            # prepare data for IPC
            arr = np.random.random(10)
            devarr = cuda.to_device(arr)
            ipch = devarr.get_ipc_handle()

            # spawn new process for testing
            ctx = mp.get_context('spawn')
            result_queue = ctx.Queue()
            args = (ipch, device_num, result_queue)
            proc = ctx.Process(target=staged_ipc_array_test, args=args)
            proc.start()
            succ, out = result_queue.get()
            proc.join(3)
            if not succ:
                self.fail(out)
            else:
                np.testing.assert_equal(arr, out)


@windows_only
@skip_on_cudasim('Ipc not available in CUDASIM')
class TestIpcNotSupported(ContextResettingTestCase):
    def test_unsupported(self):
        arr = np.arange(10, dtype=np.intp)
        devarr = cuda.to_device(arr)
        with self.assertRaises(OSError) as raises:
            devarr.get_ipc_handle()
        errmsg = str(raises.exception)
        self.assertIn('OS does not support CUDA IPC', errmsg)


if __name__ == '__main__':
    unittest.main()
