"""CUDA Driver

- Driver API binding
- NVVM API binding
- Device array implementation

"""
from numba.core import config
assert not config.ENABLE_CUDASIM, 'Cannot use real driver API with simulator'
