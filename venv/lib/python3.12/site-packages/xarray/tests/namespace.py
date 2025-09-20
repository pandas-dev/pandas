from xarray.core import duck_array_ops


def reshape(array, shape, **kwargs):
    return type(array)(duck_array_ops.reshape(array.array, shape=shape, **kwargs))
