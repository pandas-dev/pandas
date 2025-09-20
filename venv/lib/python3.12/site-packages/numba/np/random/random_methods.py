import numpy as np

from numba import uint64, uint32, uint16, uint8
from numba.core.extending import register_jitable

from numba.np.random._constants import (UINT32_MAX, UINT64_MAX,
                                        UINT16_MAX, UINT8_MAX)
from numba.np.random.generator_core import next_uint32, next_uint64

# All following implementations are direct translations from:
# https://github.com/numpy/numpy/blob/7cfef93c77599bd387ecc6a15d186c5a46024dac/numpy/random/src/distributions/distributions.c


@register_jitable
def gen_mask(max):
    mask = uint64(max)
    mask |= mask >> 1
    mask |= mask >> 2
    mask |= mask >> 4
    mask |= mask >> 8
    mask |= mask >> 16
    mask |= mask >> 32
    return mask


@register_jitable
def buffered_bounded_bool(bitgen, off, rng, bcnt, buf):
    if (rng == 0):
        return off, bcnt, buf
    if not bcnt:
        buf = next_uint32(bitgen)
        bcnt = 31
    else:
        buf >>= 1
        bcnt -= 1

    return ((buf & 1) != 0), bcnt, buf


@register_jitable
def buffered_uint8(bitgen, bcnt, buf):
    if not bcnt:
        buf = next_uint32(bitgen)
        bcnt = 3
    else:
        buf >>= 8
        bcnt -= 1

    return uint8(buf), bcnt, buf


@register_jitable
def buffered_uint16(bitgen, bcnt, buf):
    if not bcnt:
        buf = next_uint32(bitgen)
        bcnt = 1
    else:
        buf >>= 16
        bcnt -= 1

    return uint16(buf), bcnt, buf


# The following implementations use Lemire's algorithm:
# https://arxiv.org/abs/1805.10941
@register_jitable
def buffered_bounded_lemire_uint8(bitgen, rng, bcnt, buf):
    """
    Generates a random unsigned 8 bit integer bounded
    within a given interval using Lemire's rejection.

    The buffer acts as storage for a 32 bit integer
    drawn from the associated BitGenerator so that
    multiple integers of smaller bitsize can be generated
    from a single draw of the BitGenerator.
    """
    # Note: `rng` should not be 0xFF. When this happens `rng_excl` becomes
    # zero.
    rng_excl = uint8(rng) + uint8(1)

    assert (rng != 0xFF)

    # Generate a scaled random number.
    n, bcnt, buf = buffered_uint8(bitgen, bcnt, buf)
    m = uint16(n * rng_excl)

    # Rejection sampling to remove any bias
    leftover = m & 0xFF

    if (leftover < rng_excl):
        # `rng_excl` is a simple upper bound for `threshold`.
        threshold = ((uint8(UINT8_MAX) - rng) % rng_excl)

        while (leftover < threshold):
            n, bcnt, buf = buffered_uint8(bitgen, bcnt, buf)
            m = uint16(n * rng_excl)
            leftover = m & 0xFF

    return m >> 8, bcnt, buf


@register_jitable
def buffered_bounded_lemire_uint16(bitgen, rng, bcnt, buf):
    """
    Generates a random unsigned 16 bit integer bounded
    within a given interval using Lemire's rejection.

    The buffer acts as storage for a 32 bit integer
    drawn from the associated BitGenerator so that
    multiple integers of smaller bitsize can be generated
    from a single draw of the BitGenerator.
    """
    # Note: `rng` should not be 0xFFFF. When this happens `rng_excl` becomes
    # zero.
    rng_excl = uint16(rng) + uint16(1)

    assert (rng != 0xFFFF)

    # Generate a scaled random number.
    n, bcnt, buf = buffered_uint16(bitgen, bcnt, buf)
    m = uint32(n * rng_excl)

    # Rejection sampling to remove any bias
    leftover = m & 0xFFFF

    if (leftover < rng_excl):
        # `rng_excl` is a simple upper bound for `threshold`.
        threshold = ((uint16(UINT16_MAX) - rng) % rng_excl)

        while (leftover < threshold):
            n, bcnt, buf = buffered_uint16(bitgen, bcnt, buf)
            m = uint32(n * rng_excl)
            leftover = m & 0xFFFF

    return m >> 16, bcnt, buf


@register_jitable
def buffered_bounded_lemire_uint32(bitgen, rng):
    """
    Generates a random unsigned 32 bit integer bounded
    within a given interval using Lemire's rejection.
    """
    rng_excl = uint32(rng) + uint32(1)

    assert (rng != 0xFFFFFFFF)

    # Generate a scaled random number.
    m = uint64(next_uint32(bitgen)) * uint64(rng_excl)

    # Rejection sampling to remove any bias
    leftover = m & 0xFFFFFFFF

    if (leftover < rng_excl):
        # `rng_excl` is a simple upper bound for `threshold`.
        threshold = (UINT32_MAX - rng) % rng_excl

        while (leftover < threshold):
            m = uint64(next_uint32(bitgen)) * uint64(rng_excl)
            leftover = m & 0xFFFFFFFF

    return (m >> 32)


@register_jitable
def bounded_lemire_uint64(bitgen, rng):
    """
    Generates a random unsigned 64 bit integer bounded
    within a given interval using Lemire's rejection.
    """
    rng_excl = uint64(rng) + uint64(1)

    assert (rng != 0xFFFFFFFFFFFFFFFF)

    x = next_uint64(bitgen)

    leftover = uint64(x) * uint64(rng_excl)

    if (leftover < rng_excl):
        threshold = (UINT64_MAX - rng) % rng_excl

        while (leftover < threshold):
            x = next_uint64(bitgen)
            leftover = uint64(x) * uint64(rng_excl)

    x0 = x & uint64(0xFFFFFFFF)
    x1 = x >> 32
    rng_excl0 = rng_excl & uint64(0xFFFFFFFF)
    rng_excl1 = rng_excl >> 32
    w0 = x0 * rng_excl0
    t = x1 * rng_excl0 + (w0 >> 32)
    w1 = t & uint64(0xFFFFFFFF)
    w2 = t >> 32
    w1 += x0 * rng_excl1
    m1 = x1 * rng_excl1 + w2 + (w1 >> 32)

    return m1


@register_jitable
def random_bounded_uint64_fill(bitgen, low, rng, size, dtype):
    """
    Returns a new array of given size with 64 bit integers
    bounded by given interval.
    """
    out = np.empty(size, dtype=dtype)
    if rng == 0:
        for i in np.ndindex(size):
            out[i] = low
    elif rng <= 0xFFFFFFFF:
        if (rng == 0xFFFFFFFF):
            for i in np.ndindex(size):
                out[i] = low + next_uint32(bitgen)
        else:
            for i in np.ndindex(size):
                out[i] = low + buffered_bounded_lemire_uint32(bitgen, rng)

    elif (rng == 0xFFFFFFFFFFFFFFFF):
        for i in np.ndindex(size):
            out[i] = low + next_uint64(bitgen)
    else:
        for i in np.ndindex(size):
            out[i] = low + bounded_lemire_uint64(bitgen, rng)

    return out


@register_jitable
def random_bounded_uint32_fill(bitgen, low, rng, size, dtype):
    """
    Returns a new array of given size with 32 bit integers
    bounded by given interval.
    """
    out = np.empty(size, dtype=dtype)
    if rng == 0:
        for i in np.ndindex(size):
            out[i] = low
    elif rng == 0xFFFFFFFF:
        # Lemire32 doesn't support rng = 0xFFFFFFFF.
        for i in np.ndindex(size):
            out[i] = low + next_uint32(bitgen)
    else:
        for i in np.ndindex(size):
            out[i] = low + buffered_bounded_lemire_uint32(bitgen, rng)
    return out


@register_jitable
def random_bounded_uint16_fill(bitgen, low, rng, size, dtype):
    """
    Returns a new array of given size with 16 bit integers
    bounded by given interval.
    """
    buf = 0
    bcnt = 0

    out = np.empty(size, dtype=dtype)
    if rng == 0:
        for i in np.ndindex(size):
            out[i] = low
    elif rng == 0xFFFF:
        # Lemire16 doesn't support rng = 0xFFFF.
        for i in np.ndindex(size):
            val, bcnt, buf = buffered_uint16(bitgen, bcnt, buf)
            out[i] = low + val

    else:
        for i in np.ndindex(size):
            val, bcnt, buf = \
                buffered_bounded_lemire_uint16(bitgen, rng,
                                               bcnt, buf)
            out[i] = low + val
    return out


@register_jitable
def random_bounded_uint8_fill(bitgen, low, rng, size, dtype):
    """
    Returns a new array of given size with 8 bit integers
    bounded by given interval.
    """
    buf = 0
    bcnt = 0

    out = np.empty(size, dtype=dtype)
    if rng == 0:
        for i in np.ndindex(size):
            out[i] = low
    elif rng == 0xFF:
        # Lemire8 doesn't support rng = 0xFF.
        for i in np.ndindex(size):
            val, bcnt, buf = buffered_uint8(bitgen, bcnt, buf)
            out[i] = low + val
    else:
        for i in np.ndindex(size):
            val, bcnt, buf = \
                buffered_bounded_lemire_uint8(bitgen, rng,
                                              bcnt, buf)
            out[i] = low + val
    return out


@register_jitable
def random_bounded_bool_fill(bitgen, low, rng, size, dtype):
    """
    Returns a new array of given size with boolean values.
    """
    buf = 0
    bcnt = 0
    out = np.empty(size, dtype=dtype)
    for i in np.ndindex(size):
        val, bcnt, buf = buffered_bounded_bool(bitgen, low, rng, bcnt, buf)
        out[i] = low + val
    return out


@register_jitable
def _randint_arg_check(low, high, endpoint, lower_bound, upper_bound):
    """
    Check that low and high are within the bounds
    for the given datatype.
    """

    if low < lower_bound:
        raise ValueError("low is out of bounds")

    # This is being done to avoid high being accidentally
    # casted to int64/32 while subtracting 1 before
    # checking bounds, avoids overflow.
    if high > 0:
        high = uint64(high)
        if not endpoint:
            high -= uint64(1)
        upper_bound = uint64(upper_bound)
        if low > 0:
            low = uint64(low)
        if high > upper_bound:
            raise ValueError("high is out of bounds")
        if low > high:  # -1 already subtracted, closed interval
            raise ValueError("low is greater than high in given interval")
    else:
        if high > upper_bound:
            raise ValueError("high is out of bounds")
        if low > high:  # -1 already subtracted, closed interval
            raise ValueError("low is greater than high in given interval")


@register_jitable
def random_interval(bitgen, max_val):
    if (max_val == 0):
        return 0

    max_val = uint64(max_val)
    mask = uint64(gen_mask(max_val))

    if (max_val <= 0xffffffff):
        value = uint64(next_uint32(bitgen)) & mask
        while value > max_val:
            value = uint64(next_uint32(bitgen)) & mask
    else:
        value = next_uint64(bitgen) & mask
        while value > max_val:
            value = next_uint64(bitgen) & mask

    return uint64(value)
