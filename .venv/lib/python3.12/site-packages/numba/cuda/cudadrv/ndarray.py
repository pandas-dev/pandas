from numba.cuda.cudadrv import devices, driver
from numba.core.registry import cpu_target


def _calc_array_sizeof(ndim):
    """
    Use the ABI size in the CPU target
    """
    ctx = cpu_target.target_context
    return ctx.calc_array_sizeof(ndim)


def ndarray_device_allocate_data(ary):
    """
    Allocate gpu data buffer
    """
    datasize = driver.host_memory_size(ary)
    # allocate
    gpu_data = devices.get_context().memalloc(datasize)
    return gpu_data
