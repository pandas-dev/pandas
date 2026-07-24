"""
Added for symmetry with the core API
"""

from numba.core.extending import intrinsic as _intrinsic

intrinsic = _intrinsic(target='cuda')
