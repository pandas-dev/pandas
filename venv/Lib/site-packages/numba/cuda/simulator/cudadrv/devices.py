import numpy as np
from collections import namedtuple

_MemoryInfo = namedtuple("_MemoryInfo", "free,total")

_SIMULATOR_CC = (5, 2)


class FakeCUDADevice:
    def __init__(self):
        self.uuid = 'GPU-00000000-0000-0000-0000-000000000000'

    @property
    def compute_capability(self):
        return _SIMULATOR_CC


class FakeCUDAContext:
    '''
    This stub implements functionality only for simulating a single GPU
    at the moment.
    '''
    def __init__(self, device_id):
        self._device_id = device_id
        self._device = FakeCUDADevice()

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __str__(self):
        return "<Managed Device {self.id}>".format(self=self)

    @property
    def id(self):
        return self._device_id

    @property
    def device(self):
        return self._device

    @property
    def compute_capability(self):
        return _SIMULATOR_CC

    def reset(self):
        pass

    def get_memory_info(self):
        """
        Cross-platform free / total host memory is hard without external
        dependencies, e.g. `psutil` - so return infinite memory to maintain API
        type compatibility
        """
        return _MemoryInfo(float('inf'), float('inf'))

    def memalloc(self, sz):
        """
        Allocates memory on the simulated device
        At present, there is no division between simulated
        host memory and simulated device memory.
        """
        return np.ndarray(sz, dtype='u1')

    def memhostalloc(self, sz, mapped=False, portable=False, wc=False):
        '''Allocates memory on the host'''
        return self.memalloc(sz)


class FakeDeviceList:
    '''
    This stub implements a device list containing a single GPU. It also
    keeps track of the GPU status, i.e. whether the context is closed or not,
    which may have been set by the user calling reset()
    '''
    def __init__(self):
        self.lst = (FakeCUDAContext(0),)
        self.closed = False

    def __getitem__(self, devnum):
        self.closed = False
        return self.lst[devnum]

    def __str__(self):
        return ', '.join([str(d) for d in self.lst])

    def __iter__(self):
        return iter(self.lst)

    def __len__(self):
        return len(self.lst)

    @property
    def current(self):
        if self.closed:
            return None
        return self.lst[0]


gpus = FakeDeviceList()


def reset():
    gpus[0].closed = True


def get_context(devnum=0):
    return FakeCUDAContext(devnum)


def require_context(func):
    '''
    In the simulator, a context is always "available", so this is a no-op.
    '''
    return func
