from numba.core import types


class Dim3(types.Type):
    """
    A 3-tuple (x, y, z) representing the position of a block or thread.
    """
    def __init__(self):
        super().__init__(name='Dim3')


class GridGroup(types.Type):
    """
    The grid of all threads in a cooperative kernel launch.
    """
    def __init__(self):
        super().__init__(name='GridGroup')


dim3 = Dim3()
grid_group = GridGroup()


class CUDADispatcher(types.Dispatcher):
    """The type of CUDA dispatchers"""
    # This type exists (instead of using types.Dispatcher as the type of CUDA
    # dispatchers) so that we can have an alternative lowering for them to the
    # lowering of CPU dispatchers - the CPU target lowers all dispatchers as a
    # constant address, but we need to lower to a dummy value because it's not
    # generally valid to use the address of CUDA kernels and functions.
    #
    # Notes: it may be a bug in the CPU target that it lowers all dispatchers to
    # a constant address - it should perhaps only lower dispatchers acting as
    # first-class functions to a constant address. Even if that bug is fixed, it
    # is still probably a good idea to have a separate type for CUDA
    # dispatchers, and this type might get other differentiation from the CPU
    # dispatcher type in future.
