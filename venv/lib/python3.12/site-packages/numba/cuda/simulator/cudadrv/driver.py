'''
Most of the driver API is unsupported in the simulator, but some stubs are
provided to allow tests to import correctly.
'''


def device_memset(dst, val, size, stream=0):
    dst.view('u1')[:size].fill(bytes([val])[0])


def host_to_device(dst, src, size, stream=0):
    dst.view('u1')[:size] = src.view('u1')[:size]


def device_to_host(dst, src, size, stream=0):
    host_to_device(dst, src, size)


def device_memory_size(obj):
    return obj.itemsize * obj.size


def device_to_device(dst, src, size, stream=0):
    host_to_device(dst, src, size)


class FakeDriver(object):
    def get_device_count(self):
        return 1


driver = FakeDriver()


class Linker:
    @classmethod
    def new(cls, max_registers=0, lineinfo=False, cc=None):
        return Linker()

    @property
    def lto(self):
        return False


class LinkerError(RuntimeError):
    pass


class NvrtcError(RuntimeError):
    pass


class CudaAPIError(RuntimeError):
    pass


def launch_kernel(*args, **kwargs):
    msg = 'Launching kernels directly is not supported in the simulator'
    raise RuntimeError(msg)


USE_NV_BINDING = False
