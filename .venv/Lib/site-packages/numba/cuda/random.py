import math

from numba import (config, cuda, float32, float64, uint32, int64, uint64,
                   from_dtype, jit)

import numpy as np

# This implementation is based upon the xoroshiro128+ and splitmix64 algorithms
# described at:
#
#     http://xoroshiro.di.unimi.it/
#
# and originally implemented by David Blackman and Sebastiano Vigna.
#
# The implementations below are based on the C source code:
#
#  * http://xoroshiro.di.unimi.it/xoroshiro128plus.c
#  * http://xoroshiro.di.unimi.it/splitmix64.c
#
# Splitmix64 is used to generate the initial state of the xoroshiro128+
# generator to ensure that small seeds don't result in predictable output.

# **WARNING**: There is a lot of verbose casting in this file to ensure that
# NumPy casting conventions (which cast uint64 [op] int32 to float64) don't
# turn integers into floats when using these functions in the CUDA simulator.
#
# There are also no function type signatures to ensure that compilation is
# deferred so that import is quick, and Sphinx autodoc works.  We are also
# using the CPU @jit decorator everywhere to create functions that work as
# both CPU and CUDA device functions.

xoroshiro128p_dtype = np.dtype([('s0', np.uint64), ('s1', np.uint64)],
                               align=True)
xoroshiro128p_type = from_dtype(xoroshiro128p_dtype)

# When cudasim is enabled, Fake CUDA arrays are passed to some of the
# @jit-decorated functions. This required fallback to object mode. With
# Numba 0.59.0 object mode must be explicitly enabled.
# https://numba.readthedocs.io/en/stable/reference/deprecation.html#deprecation-of-object-mode-fall-back-behaviour-when-using-jit
# In order to avoid the warning / future error, we explicitly specify that
# object mode with loop lifting is acceptable when using the simulator.
_forceobj = _looplift = config.ENABLE_CUDASIM
_nopython = not config.ENABLE_CUDASIM


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def init_xoroshiro128p_state(states, index, seed):
    '''Use SplitMix64 to generate an xoroshiro128p state from 64-bit seed.

    This ensures that manually set small seeds don't result in a predictable
    initial sequence from the random number generator.

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: uint64
    :param index: offset in states to update
    :type seed: int64
    :param seed: seed value to use when initializing state
    '''
    index = int64(index)
    seed = uint64(seed)

    z = seed + uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> uint32(30))) * uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> uint32(27))) * uint64(0x94D049BB133111EB)
    z = z ^ (z >> uint32(31))

    states[index]['s0'] = z
    states[index]['s1'] = z


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def rotl(x, k):
    '''Left rotate x by k bits.'''
    x = uint64(x)
    k = uint32(k)
    return (x << k) | (x >> uint32(64 - k))


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def xoroshiro128p_next(states, index):
    '''Return the next random uint64 and advance the RNG in states[index].

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: int64
    :param index: offset in states to update
    :rtype: uint64
    '''
    index = int64(index)
    s0 = states[index]['s0']
    s1 = states[index]['s1']
    result = s0 + s1

    s1 ^= s0
    states[index]['s0'] = uint64(rotl(s0, uint32(55))) ^ s1 ^ (s1 << uint32(14))
    states[index]['s1'] = uint64(rotl(s1, uint32(36)))

    return result


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def xoroshiro128p_jump(states, index):
    '''Advance the RNG in ``states[index]`` by 2**64 steps.

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: int64
    :param index: offset in states to update
    '''
    index = int64(index)

    jump = (uint64(0xbeac0467eba5facb), uint64(0xd86b048b86aa9922))

    s0 = uint64(0)
    s1 = uint64(0)

    for i in range(2):
        for b in range(64):
            if jump[i] & (uint64(1) << uint32(b)):
                s0 ^= states[index]['s0']
                s1 ^= states[index]['s1']
            xoroshiro128p_next(states, index)

    states[index]['s0'] = s0
    states[index]['s1'] = s1


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def uint64_to_unit_float64(x):
    '''Convert uint64 to float64 value in the range [0.0, 1.0)'''
    x = uint64(x)
    return (x >> uint32(11)) * (float64(1) / (uint64(1) << uint32(53)))


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def uint64_to_unit_float32(x):
    '''Convert uint64 to float32 value in the range [0.0, 1.0)'''
    x = uint64(x)
    return float32(uint64_to_unit_float64(x))


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def xoroshiro128p_uniform_float32(states, index):
    '''Return a float32 in range [0.0, 1.0) and advance ``states[index]``.

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: int64
    :param index: offset in states to update
    :rtype: float32
    '''
    index = int64(index)
    return uint64_to_unit_float32(xoroshiro128p_next(states, index))


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def xoroshiro128p_uniform_float64(states, index):
    '''Return a float64 in range [0.0, 1.0) and advance ``states[index]``.

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: int64
    :param index: offset in states to update
    :rtype: float64
    '''
    index = int64(index)
    return uint64_to_unit_float64(xoroshiro128p_next(states, index))


TWO_PI_FLOAT32 = np.float32(2 * math.pi)
TWO_PI_FLOAT64 = np.float64(2 * math.pi)


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def xoroshiro128p_normal_float32(states, index):
    '''Return a normally distributed float32 and advance ``states[index]``.

    The return value is drawn from a Gaussian of mean=0 and sigma=1 using the
    Box-Muller transform.  This advances the RNG sequence by two steps.

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: int64
    :param index: offset in states to update
    :rtype: float32
    '''
    index = int64(index)

    u1 = xoroshiro128p_uniform_float32(states, index)
    u2 = xoroshiro128p_uniform_float32(states, index)

    z0 = math.sqrt(-float32(2.0) * math.log(u1)) * math.cos(TWO_PI_FLOAT32 * u2)
    # discarding second normal value
    # z1 = math.sqrt(-float32(2.0) * math.log(u1))
    #                * math.sin(TWO_PI_FLOAT32 * u2)
    return z0


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def xoroshiro128p_normal_float64(states, index):
    '''Return a normally distributed float32 and advance ``states[index]``.

    The return value is drawn from a Gaussian of mean=0 and sigma=1 using the
    Box-Muller transform.  This advances the RNG sequence by two steps.

    :type states: 1D array, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type index: int64
    :param index: offset in states to update
    :rtype: float64
    '''
    index = int64(index)

    u1 = xoroshiro128p_uniform_float32(states, index)
    u2 = xoroshiro128p_uniform_float32(states, index)

    z0 = math.sqrt(-float64(2.0) * math.log(u1)) * math.cos(TWO_PI_FLOAT64 * u2)
    # discarding second normal value
    # z1 = math.sqrt(-float64(2.0) * math.log(u1))
    #                * math.sin(TWO_PI_FLOAT64 * u2)
    return z0


@jit(forceobj=_forceobj, looplift=_looplift, nopython=_nopython)
def init_xoroshiro128p_states_cpu(states, seed, subsequence_start):
    n = states.shape[0]
    seed = uint64(seed)
    subsequence_start = uint64(subsequence_start)

    if n >= 1:
        init_xoroshiro128p_state(states, 0, seed)

        # advance to starting subsequence number
        for _ in range(subsequence_start):
            xoroshiro128p_jump(states, 0)

        # populate the rest of the array
        for i in range(1, n):
            states[i] = states[i - 1]  # take state of previous generator
            xoroshiro128p_jump(states, i)  # and jump forward 2**64 steps


def init_xoroshiro128p_states(states, seed, subsequence_start=0, stream=0):
    '''Initialize RNG states on the GPU for parallel generators.

    This initializes the RNG states so that each state in the array corresponds
    subsequences in the separated by 2**64 steps from each other in the main
    sequence.  Therefore, as long no CUDA thread requests more than 2**64
    random numbers, all of the RNG states produced by this function are
    guaranteed to be independent.

    The subsequence_start parameter can be used to advance the first RNG state
    by a multiple of 2**64 steps.

    :type states: 1D DeviceNDArray, dtype=xoroshiro128p_dtype
    :param states: array of RNG states
    :type seed: uint64
    :param seed: starting seed for list of generators
    '''

    # Initialization on CPU is much faster than the GPU
    states_cpu = np.empty(shape=states.shape, dtype=xoroshiro128p_dtype)
    init_xoroshiro128p_states_cpu(states_cpu, seed, subsequence_start)

    states.copy_to_device(states_cpu, stream=stream)


def create_xoroshiro128p_states(n, seed, subsequence_start=0, stream=0):
    '''Returns a new device array initialized for n random number generators.

    This initializes the RNG states so that each state in the array corresponds
    subsequences in the separated by 2**64 steps from each other in the main
    sequence.  Therefore, as long no CUDA thread requests more than 2**64
    random numbers, all of the RNG states produced by this function are
    guaranteed to be independent.

    The subsequence_start parameter can be used to advance the first RNG state
    by a multiple of 2**64 steps.

    :type n: int
    :param n: number of RNG states to create
    :type seed: uint64
    :param seed: starting seed for list of generators
    :type subsequence_start: uint64
    :param subsequence_start:
    :type stream: CUDA stream
    :param stream: stream to run initialization kernel on
    '''
    states = cuda.device_array(n, dtype=xoroshiro128p_dtype, stream=stream)
    init_xoroshiro128p_states(states, seed, subsequence_start, stream)
    return states
