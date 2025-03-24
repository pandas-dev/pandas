import numbers
from numba.core.errors import LoweringError


class KernelRuntimeError(RuntimeError):
    def __init__(self, msg, tid=None, ctaid=None):
        self.tid = tid
        self.ctaid = ctaid
        self.msg = msg
        t = ("An exception was raised in thread=%s block=%s\n"
             "\t%s")
        msg = t % (self.tid, self.ctaid, self.msg)
        super(KernelRuntimeError, self).__init__(msg)


class CudaLoweringError(LoweringError):
    pass


_launch_help_url = ("https://numba.readthedocs.io/en/stable/cuda/"
                    "kernels.html#kernel-invocation")
missing_launch_config_msg = """
Kernel launch configuration was not specified. Use the syntax:

kernel_function[blockspergrid, threadsperblock](arg0, arg1, ..., argn)

See {} for help.

""".format(_launch_help_url)


def normalize_kernel_dimensions(griddim, blockdim):
    """
    Normalize and validate the user-supplied kernel dimensions.
    """

    def check_dim(dim, name):
        if not isinstance(dim, (tuple, list)):
            dim = [dim]
        else:
            dim = list(dim)
        if len(dim) > 3:
            raise ValueError('%s must be a sequence of 1, 2 or 3 integers, '
                             'got %r' % (name, dim))
        for v in dim:
            if not isinstance(v, numbers.Integral):
                raise TypeError('%s must be a sequence of integers, got %r'
                                % (name, dim))
        while len(dim) < 3:
            dim.append(1)
        return tuple(dim)

    if None in (griddim, blockdim):
        raise ValueError(missing_launch_config_msg)

    griddim = check_dim(griddim, 'griddim')
    blockdim = check_dim(blockdim, 'blockdim')

    return griddim, blockdim
