from contextlib import contextmanager
import functools
import sys
import threading

import numpy as np

from .cudadrv.devicearray import FakeCUDAArray, FakeWithinKernelCUDAArray
from .kernelapi import Dim3, FakeCUDAModule, swapped_cuda_module
from ..errors import normalize_kernel_dimensions
from ..args import wrap_arg, ArgHint


"""
Global variable to keep track of the current "kernel context", i.e the
FakeCUDAModule.  We only support one kernel launch at a time.
No support for concurrent kernel launch.
"""
_kernel_context = None


@contextmanager
def _push_kernel_context(mod):
    """
    Push the current kernel context.
    """
    global _kernel_context
    assert _kernel_context is None, "concurrent simulated kernel not supported"
    _kernel_context = mod
    try:
        yield
    finally:
        _kernel_context = None


def _get_kernel_context():
    """
    Get the current kernel context. This is usually done by a device function.
    """
    return _kernel_context


class FakeOverload:
    '''
    Used only to provide the max_cooperative_grid_blocks method
    '''
    def max_cooperative_grid_blocks(self, blockdim):
        # We can only run one block in a cooperative grid because we have no
        # mechanism for synchronization between different blocks
        return 1


class FakeOverloadDict(dict):
    def __getitem__(self, key):
        # Always return a fake overload for any signature, as we don't keep
        # track of overloads in the simulator.
        return FakeOverload()


class FakeCUDAKernel(object):
    '''
    Wraps a @cuda.jit-ed function.
    '''

    def __init__(
        self, fn, device, fastmath=False, extensions=None, debug=False
    ):
        if extensions is None:
            extensions = []
        self.fn = fn
        self._device = device
        self._fastmath = fastmath
        self._debug = debug
        self.extensions = list(extensions) # defensive copy
        # Initial configuration: grid unconfigured, stream 0, no dynamic shared
        # memory.
        self.grid_dim = None
        self.block_dim = None
        self.stream = 0
        self.dynshared_size = 0
        functools.update_wrapper(self, fn)

    def __call__(self, *args):
        if self._device:
            with swapped_cuda_module(self.fn, _get_kernel_context()):
                return self.fn(*args)

        # Ensure we've been given a valid grid configuration
        grid_dim, block_dim = normalize_kernel_dimensions(self.grid_dim,
                                                          self.block_dim)

        fake_cuda_module = FakeCUDAModule(grid_dim, block_dim,
                                          self.dynshared_size)
        with _push_kernel_context(fake_cuda_module):
            # fake_args substitutes all numpy arrays for FakeCUDAArrays
            # because they implement some semantics differently
            retr = []

            def fake_arg(arg):
                # map the arguments using any extension you've registered
                _, arg = functools.reduce(
                    lambda ty_val, extension: extension.prepare_args(
                        *ty_val,
                        stream=0,
                        retr=retr),
                    self.extensions,
                    (None, arg)
                )

                if isinstance(arg, np.ndarray) and arg.ndim > 0:
                    ret = wrap_arg(arg).to_device(retr)
                elif isinstance(arg, ArgHint):
                    ret = arg.to_device(retr)
                elif isinstance(arg, np.void):
                    ret = FakeCUDAArray(arg)  # In case a np record comes in.
                else:
                    ret = arg
                if isinstance(ret, FakeCUDAArray):
                    return FakeWithinKernelCUDAArray(ret)
                return ret

            fake_args = [fake_arg(arg) for arg in args]
            with swapped_cuda_module(self.fn, fake_cuda_module):
                # Execute one block at a time
                for grid_point in np.ndindex(*grid_dim):
                    bm = BlockManager(self.fn, grid_dim, block_dim, self._debug)
                    bm.run(grid_point, *fake_args)

            for wb in retr:
                wb()

    def __getitem__(self, configuration):
        self.grid_dim, self.block_dim = \
            normalize_kernel_dimensions(*configuration[:2])

        if len(configuration) == 4:
            self.dynshared_size = configuration[3]

        return self

    def bind(self):
        pass

    def specialize(self, *args):
        return self

    def forall(self, ntasks, tpb=0, stream=0, sharedmem=0):
        if ntasks < 0:
            raise ValueError("Can't create ForAll with negative task count: %s"
                             % ntasks)
        return self[ntasks, 1, stream, sharedmem]

    @property
    def overloads(self):
        return FakeOverloadDict()

    @property
    def py_func(self):
        return self.fn


# Thread emulation

class BlockThread(threading.Thread):
    '''
    Manages the execution of a function for a single CUDA thread.
    '''
    def __init__(self, f, manager, blockIdx, threadIdx, debug):
        if debug:
            def debug_wrapper(*args, **kwargs):
                np.seterr(divide='raise')
                f(*args, **kwargs)
            target = debug_wrapper
        else:
            target = f

        super(BlockThread, self).__init__(target=target)
        self.syncthreads_event = threading.Event()
        self.syncthreads_blocked = False
        self._manager = manager
        self.blockIdx = Dim3(*blockIdx)
        self.threadIdx = Dim3(*threadIdx)
        self.exception = None
        self.daemon = True
        self.abort = False
        self.debug = debug
        blockDim = Dim3(*self._manager._block_dim)
        self.thread_id = self.threadIdx.x + (blockDim.x * (self.threadIdx.y +
                                                           blockDim.y *
                                                           self.threadIdx.z))

    def run(self):
        try:
            super(BlockThread, self).run()
        except Exception as e:
            tid = 'tid=%s' % list(self.threadIdx)
            ctaid = 'ctaid=%s' % list(self.blockIdx)
            if str(e) == '':
                msg = '%s %s' % (tid, ctaid)
            else:
                msg = '%s %s: %s' % (tid, ctaid, e)
            tb = sys.exc_info()[2]
            # Using `with_traceback` here would cause it to be mutated by
            # future raise statements, which may or may not matter.
            self.exception = (type(e)(msg), tb)

    def syncthreads(self):

        if self.abort:
            raise RuntimeError("abort flag set on syncthreads call")

        self.syncthreads_blocked = True
        self.syncthreads_event.wait()
        self.syncthreads_event.clear()

        if self.abort:
            raise RuntimeError("abort flag set on syncthreads clear")

    def syncthreads_count(self, value):
        idx = self.threadIdx.x, self.threadIdx.y, self.threadIdx.z
        self._manager.block_state[idx] = value
        self.syncthreads()
        count = np.count_nonzero(self._manager.block_state)
        self.syncthreads()
        return count

    def syncthreads_and(self, value):
        idx = self.threadIdx.x, self.threadIdx.y, self.threadIdx.z
        self._manager.block_state[idx] = value
        self.syncthreads()
        test = np.all(self._manager.block_state)
        self.syncthreads()
        return 1 if test else 0

    def syncthreads_or(self, value):
        idx = self.threadIdx.x, self.threadIdx.y, self.threadIdx.z
        self._manager.block_state[idx] = value
        self.syncthreads()
        test = np.any(self._manager.block_state)
        self.syncthreads()
        return 1 if test else 0

    def __str__(self):
        return 'Thread <<<%s, %s>>>' % (self.blockIdx, self.threadIdx)


class BlockManager(object):
    '''
    Manages the execution of a thread block.

    When run() is called, all threads are started. Each thread executes until it
    hits syncthreads(), at which point it sets its own syncthreads_blocked to
    True so that the BlockManager knows it is blocked. It then waits on its
    syncthreads_event.

    The BlockManager polls threads to determine if they are blocked in
    syncthreads(). If it finds a blocked thread, it adds it to the set of
    blocked threads. When all threads are blocked, it unblocks all the threads.
    The thread are unblocked by setting their syncthreads_blocked back to False
    and setting their syncthreads_event.

    The polling continues until no threads are alive, when execution is
    complete.
    '''
    def __init__(self, f, grid_dim, block_dim, debug):
        self._grid_dim = grid_dim
        self._block_dim = block_dim
        self._f = f
        self._debug = debug
        self.block_state = np.zeros(block_dim, dtype=np.bool_)

    def run(self, grid_point, *args):
        # Create all threads
        threads = set()
        livethreads = set()
        blockedthreads = set()
        for block_point in np.ndindex(*self._block_dim):
            def target():
                self._f(*args)
            t = BlockThread(target, self, grid_point, block_point, self._debug)
            t.start()
            threads.add(t)
            livethreads.add(t)

        # Potential optimisations:
        # 1. Continue the while loop immediately after finding a blocked thread
        # 2. Don't poll already-blocked threads
        while livethreads:
            for t in livethreads:
                if t.syncthreads_blocked:
                    blockedthreads.add(t)
                elif t.exception:

                    # Abort all other simulator threads on exception,
                    # do *not* join immediately to facilitate debugging.
                    for t_other in threads:
                        t_other.abort = True
                        t_other.syncthreads_blocked = False
                        t_other.syncthreads_event.set()

                    raise t.exception[0].with_traceback(t.exception[1])
            if livethreads == blockedthreads:
                for t in blockedthreads:
                    t.syncthreads_blocked = False
                    t.syncthreads_event.set()
                blockedthreads = set()
            livethreads = set([ t for t in livethreads if t.is_alive() ])
        # Final check for exceptions in case any were set prior to thread
        # finishing, before we could check it
        for t in threads:
            if t.exception:
                raise t.exception[0].with_traceback(t.exception[1])
