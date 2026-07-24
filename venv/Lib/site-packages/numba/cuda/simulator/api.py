'''
Contains CUDA API functions
'''

# Imports here bring together parts of the API from other modules, so some of
# them appear unused.
from contextlib import contextmanager

from .cudadrv.devices import require_context, reset, gpus  # noqa: F401
from .kernel import FakeCUDAKernel
from numba.core.sigutils import is_signature
from warnings import warn
from ..args import In, Out, InOut  # noqa: F401


def select_device(dev=0):
    assert dev == 0, 'Only a single device supported by the simulator'


def is_float16_supported():
    return True


class stream(object):
    '''
    The stream API is supported in the simulator - however, all execution
    occurs synchronously, so synchronization requires no operation.
    '''
    @contextmanager
    def auto_synchronize(self):
        yield

    def synchronize(self):
        pass


def synchronize():
    pass


def close():
    gpus.closed = True


def declare_device(*args, **kwargs):
    pass


def detect():
    print('Found 1 CUDA devices')
    print('id %d    %20s %40s' % (0, 'SIMULATOR', '[SUPPORTED]'))
    print('%40s: 5.0' % 'compute capability')


def list_devices():
    return gpus


# Events

class Event(object):
    '''
    The simulator supports the event API, but they do not record timing info,
    and all simulation is synchronous. Execution time is not recorded.
    '''
    def record(self, stream=0):
        pass

    def wait(self, stream=0):
        pass

    def synchronize(self):
        pass

    def elapsed_time(self, event):
        warn('Simulator timings are bogus')
        return 0.0


event = Event


def jit(func_or_sig=None, device=False, debug=False, argtypes=None,
        inline=False, restype=None, fastmath=False, link=None,
        boundscheck=None, opt=True, cache=None
        ):
    # Here for API compatibility
    if boundscheck:
        raise NotImplementedError("bounds checking is not supported for CUDA")

    if link is not None:
        raise NotImplementedError('Cannot link PTX in the simulator')

    # Check for first argument specifying types - in that case the
    # decorator is not being passed a function
    if (func_or_sig is None or is_signature(func_or_sig)
            or isinstance(func_or_sig, list)):
        def jitwrapper(fn):
            return FakeCUDAKernel(fn,
                                  device=device,
                                  fastmath=fastmath,
                                  debug=debug)
        return jitwrapper
    return FakeCUDAKernel(func_or_sig, device=device, debug=debug)


@contextmanager
def defer_cleanup():
    # No effect for simulator
    yield
