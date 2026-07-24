from .decorators import jit
import numba


@jit(device=True)
def all_sync(mask, predicate):
    """
    If for all threads in the masked warp the predicate is true, then
    a non-zero value is returned, otherwise 0 is returned.
    """
    return numba.cuda.vote_sync_intrinsic(mask, 0, predicate)[1]


@jit(device=True)
def any_sync(mask, predicate):
    """
    If for any thread in the masked warp the predicate is true, then
    a non-zero value is returned, otherwise 0 is returned.
    """
    return numba.cuda.vote_sync_intrinsic(mask, 1, predicate)[1]


@jit(device=True)
def eq_sync(mask, predicate):
    """
    If for all threads in the masked warp the boolean predicate is the same,
    then a non-zero value is returned, otherwise 0 is returned.
    """
    return numba.cuda.vote_sync_intrinsic(mask, 2, predicate)[1]


@jit(device=True)
def ballot_sync(mask, predicate):
    """
    Returns a mask of all threads in the warp whose predicate is true,
    and are within the given mask.
    """
    return numba.cuda.vote_sync_intrinsic(mask, 3, predicate)[0]


@jit(device=True)
def shfl_sync(mask, value, src_lane):
    """
    Shuffles value across the masked warp and returns the value
    from src_lane. If this is outside the warp, then the
    given value is returned.
    """
    return numba.cuda.shfl_sync_intrinsic(mask, 0, value, src_lane, 0x1f)[0]


@jit(device=True)
def shfl_up_sync(mask, value, delta):
    """
    Shuffles value across the masked warp and returns the value
    from (laneid - delta). If this is outside the warp, then the
    given value is returned.
    """
    return numba.cuda.shfl_sync_intrinsic(mask, 1, value, delta, 0)[0]


@jit(device=True)
def shfl_down_sync(mask, value, delta):
    """
    Shuffles value across the masked warp and returns the value
    from (laneid + delta). If this is outside the warp, then the
    given value is returned.
    """
    return numba.cuda.shfl_sync_intrinsic(mask, 2, value, delta, 0x1f)[0]


@jit(device=True)
def shfl_xor_sync(mask, value, lane_mask):
    """
    Shuffles value across the masked warp and returns the value
    from (laneid ^ lane_mask).
    """
    return numba.cuda.shfl_sync_intrinsic(mask, 3, value, lane_mask, 0x1f)[0]
