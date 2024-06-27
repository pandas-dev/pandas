"""
Implementation of method overloads for Generator objects.
"""

import numpy as np
from numba.core import types
from numba.core.extending import overload_method, register_jitable
from numba.np.numpy_support import as_dtype, from_dtype
from numba.np.random.generator_core import next_float, next_double
from numba.np.numpy_support import is_nonelike
from numba.core.errors import TypingError
from numba.core.types.containers import Tuple, UniTuple
from numba.np.random.distributions import \
    (random_standard_exponential_inv_f, random_standard_exponential_inv,
     random_standard_exponential, random_standard_normal_f,
     random_standard_gamma, random_standard_normal, random_uniform,
     random_standard_exponential_f, random_standard_gamma_f, random_normal,
     random_exponential, random_gamma, random_beta, random_power,
     random_f,random_chisquare,random_standard_cauchy,random_pareto,
     random_weibull, random_laplace, random_logistic,
     random_lognormal, random_rayleigh, random_standard_t, random_wald,
     random_geometric, random_zipf, random_triangular,
     random_poisson, random_negative_binomial, random_logseries,
     random_noncentral_chisquare, random_noncentral_f, random_binomial)
from numba.np.random import random_methods


def _get_proper_func(func_32, func_64, dtype, dist_name="the given"):
    """
        Most of the standard NumPy distributions that accept dtype argument
        only support either np.float32 or np.float64 as dtypes.

        This is a helper function that helps Numba select the proper underlying
        implementation according to provided dtype.
    """
    if isinstance(dtype, types.Omitted):
        dtype = dtype.value

    np_dt = dtype
    if isinstance(dtype, type):
        nb_dt = from_dtype(np.dtype(dtype))
    elif isinstance(dtype, types.NumberClass):
        nb_dt = dtype
        np_dt = as_dtype(nb_dt)

    if np_dt not in [np.float32, np.float64]:
        raise TypingError("Argument dtype is not one of the" +
                          " expected type(s): " +
                          " np.float32 or np.float64")

    if np_dt == np.float32:
        next_func = func_32
    else:
        next_func = func_64

    return next_func, nb_dt


def check_size(size):
    if not any([isinstance(size, UniTuple) and
                isinstance(size.dtype, types.Integer),
                isinstance(size, Tuple) and size.count == 0,
                isinstance(size, types.Integer)]):
        raise TypingError("Argument size is not one of the" +
                          " expected type(s): " +
                          " an integer, an empty tuple or a tuple of integers")


def check_types(obj, type_list, arg_name):
    """
    Check if given object is one of the provided types.
    If not raises an TypeError
    """
    if isinstance(obj, types.Omitted):
        obj = obj.value

    if not isinstance(type_list, (list, tuple)):
        type_list = [type_list]

    if not any([isinstance(obj, _type) for _type in type_list]):
        raise TypingError(f"Argument {arg_name} is not one of the" +
                          f" expected type(s): {type_list}")


# Overload the Generator().integers()
@overload_method(types.NumPyRandomGeneratorType, 'integers')
def NumPyRandomGeneratorType_integers(inst, low, high, size=None,
                                      dtype=np.int64, endpoint=False):
    check_types(low, [types.Integer,
                      types.Boolean, bool, int], 'low')
    check_types(high, [types.Integer, types.Boolean,
                       bool, int], 'high')
    check_types(endpoint, [types.Boolean, bool], 'endpoint')

    if isinstance(size, types.Omitted):
        size = size.value

    if isinstance(dtype, types.Omitted):
        dtype = dtype.value

    if isinstance(dtype, type):
        nb_dt = from_dtype(np.dtype(dtype))
        _dtype = dtype
    elif isinstance(dtype, types.NumberClass):
        nb_dt = dtype
        _dtype = as_dtype(nb_dt)
    else:
        raise TypingError("Argument dtype is not one of the" +
                          " expected type(s): " +
                          "np.int32, np.int64, np.int16, np.int8, "
                          "np.uint32, np.uint64, np.uint16, np.uint8, "
                          "np.bool_")

    if _dtype == np.bool_:
        int_func = random_methods.random_bounded_bool_fill
        lower_bound = -1
        upper_bound = 2
    else:
        try:
            i_info = np.iinfo(_dtype)
        except ValueError:
            raise TypingError("Argument dtype is not one of the" +
                              " expected type(s): " +
                              "np.int32, np.int64, np.int16, np.int8, "
                              "np.uint32, np.uint64, np.uint16, np.uint8, "
                              "np.bool_")
        int_func = getattr(random_methods,
                           f'random_bounded_uint{i_info.bits}_fill')
        lower_bound = i_info.min
        upper_bound = i_info.max

    if is_nonelike(size):
        def impl(inst, low, high, size=None,
                 dtype=np.int64, endpoint=False):
            random_methods._randint_arg_check(low, high, endpoint,
                                              lower_bound, upper_bound)
            if not endpoint:
                high -= dtype(1)
                low = dtype(low)
                high = dtype(high)
                rng = high - low
                return int_func(inst.bit_generator, low, rng, 1, dtype)[0]
            else:
                low = dtype(low)
                high = dtype(high)
                rng = high - low
                return int_func(inst.bit_generator, low, rng, 1, dtype)[0]
        return impl
    else:
        check_size(size)

        def impl(inst, low, high, size=None,
                 dtype=np.int64, endpoint=False):
            random_methods._randint_arg_check(low, high, endpoint,
                                              lower_bound, upper_bound)
            if not endpoint:
                high -= dtype(1)
                low = dtype(low)
                high = dtype(high)
                rng = high - low
                return int_func(inst.bit_generator, low, rng, size, dtype)
            else:
                low = dtype(low)
                high = dtype(high)
                rng = high - low
                return int_func(inst.bit_generator, low, rng, size, dtype)
        return impl


# The following `shuffle` implementation is a direct translation from:
# https://github.com/numpy/numpy/blob/95e3e7f445407e4f355b23d6a9991d8774f0eb0c/numpy/random/_generator.pyx#L4578

# Overload the Generator().shuffle()
@overload_method(types.NumPyRandomGeneratorType, 'shuffle')
def NumPyRandomGeneratorType_shuffle(inst, x, axis=0):
    check_types(x, [types.Array], 'x')
    check_types(axis, [int, types.Integer], 'axis')

    def impl(inst, x, axis=0):
        if axis < 0:
            axis = axis + x.ndim
        if axis > x.ndim - 1 or axis < 0:
            raise IndexError("Axis is out of bounds for the given array")

        z = np.swapaxes(x, 0, axis)
        buf = np.empty_like(z[0, ...])

        for i in range(len(z) - 1, 0, -1):
            j = types.intp(random_methods.random_interval(inst.bit_generator,
                                                          i))
            if i == j:
                continue
            buf[...] = z[j, ...]
            z[j, ...] = z[i, ...]
            z[i, ...] = buf

    return impl


# The following `permutation` implementation is a direct translation from:
# https://github.com/numpy/numpy/blob/95e3e7f445407e4f355b23d6a9991d8774f0eb0c/numpy/random/_generator.pyx#L4710
# Overload the Generator().permutation()
@overload_method(types.NumPyRandomGeneratorType, 'permutation')
def NumPyRandomGeneratorType_permutation(inst, x, axis=0):
    check_types(x, [types.Array, types.Integer], 'x')
    check_types(axis, [int, types.Integer], 'axis')

    IS_INT = isinstance(x, types.Integer)

    def impl(inst, x, axis=0):
        if IS_INT:
            new_arr = np.arange(x)
            # NumPy ignores the axis argument when x is an integer
            inst.shuffle(new_arr)
        else:
            new_arr = x.copy()
            inst.shuffle(new_arr, axis=axis)
        return new_arr

    return impl


# Overload the Generator().random()
@overload_method(types.NumPyRandomGeneratorType, 'random')
def NumPyRandomGeneratorType_random(inst, size=None, dtype=np.float64):
    dist_func, nb_dt = _get_proper_func(next_float, next_double,
                                        dtype, "random")
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, size=None, dtype=np.float64):
            return nb_dt(dist_func(inst.bit_generator))
        return impl
    else:
        check_size(size)

        def impl(inst, size=None, dtype=np.float64):
            out = np.empty(size, dtype=dtype)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = dist_func(inst.bit_generator)
            return out
        return impl


# Overload the Generator().standard_exponential() method
@overload_method(types.NumPyRandomGeneratorType, 'standard_exponential')
def NumPyRandomGeneratorType_standard_exponential(inst, size=None,
                                                  dtype=np.float64,
                                                  method='zig'):
    check_types(method, [types.UnicodeType, str], 'method')
    dist_func_inv, nb_dt = _get_proper_func(
        random_standard_exponential_inv_f,
        random_standard_exponential_inv,
        dtype
    )

    dist_func, nb_dt = _get_proper_func(random_standard_exponential_f,
                                        random_standard_exponential,
                                        dtype)

    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, size=None, dtype=np.float64, method='zig'):
            if method == 'zig':
                return nb_dt(dist_func(inst.bit_generator))
            elif method == 'inv':
                return nb_dt(dist_func_inv(inst.bit_generator))
            else:
                raise ValueError("Method must be either 'zig' or 'inv'")
        return impl
    else:
        check_size(size)

        def impl(inst, size=None, dtype=np.float64, method='zig'):
            out = np.empty(size, dtype=dtype)
            out_f = out.flat
            if method == 'zig':
                for i in range(out.size):
                    out_f[i] = dist_func(inst.bit_generator)
            elif method == 'inv':
                for i in range(out.size):
                    out_f[i] = dist_func_inv(inst.bit_generator)
            else:
                raise ValueError("Method must be either 'zig' or 'inv'")
            return out
        return impl


# Overload the Generator().standard_normal() method
@overload_method(types.NumPyRandomGeneratorType, 'standard_normal')
def NumPyRandomGeneratorType_standard_normal(inst, size=None, dtype=np.float64):
    dist_func, nb_dt = _get_proper_func(random_standard_normal_f,
                                        random_standard_normal,
                                        dtype)
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, size=None, dtype=np.float64):
            return nb_dt(dist_func(inst.bit_generator))
        return impl
    else:
        check_size(size)

        def impl(inst, size=None, dtype=np.float64):
            out = np.empty(size, dtype=dtype)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = dist_func(inst.bit_generator)
            return out
        return impl


# Overload the Generator().standard_gamma() method
@overload_method(types.NumPyRandomGeneratorType, 'standard_gamma')
def NumPyRandomGeneratorType_standard_gamma(inst, shape, size=None,
                                            dtype=np.float64):
    check_types(shape, [types.Float, types.Integer, int, float], 'shape')
    dist_func, nb_dt = _get_proper_func(random_standard_gamma_f,
                                        random_standard_gamma,
                                        dtype)
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, shape, size=None, dtype=np.float64):
            return nb_dt(dist_func(inst.bit_generator, shape))
        return impl
    else:
        check_size(size)

        def impl(inst, shape, size=None, dtype=np.float64):
            out = np.empty(size, dtype=dtype)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = dist_func(inst.bit_generator, shape)
            return out
        return impl


# Overload the Generator().normal() method
@overload_method(types.NumPyRandomGeneratorType, 'normal')
def NumPyRandomGeneratorType_normal(inst, loc=0.0, scale=1.0,
                                    size=None):
    check_types(loc, [types.Float, types.Integer, int, float], 'loc')
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, loc=0.0, scale=1.0, size=None):
            return random_normal(inst.bit_generator, loc, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, loc=0.0, scale=1.0, size=None):
            out = np.empty(size, dtype=np.float64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_normal(inst.bit_generator, loc, scale)
            return out
        return impl


# Overload the Generator().uniform() method
@overload_method(types.NumPyRandomGeneratorType, 'uniform')
def NumPyRandomGeneratorType_uniform(inst, low=0.0, high=1.0,
                                     size=None):
    check_types(low, [types.Float, types.Integer, int, float], 'low')
    check_types(high, [types.Float, types.Integer, int, float], 'high')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, low=0.0, high=1.0, size=None):
            return random_uniform(inst.bit_generator, low, high - low)
        return impl
    else:
        check_size(size)

        def impl(inst, low=0.0, high=1.0, size=None):
            out = np.empty(size, dtype=np.float64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_uniform(inst.bit_generator, low, high - low)
            return out
        return impl


# Overload the Generator().exponential() method
@overload_method(types.NumPyRandomGeneratorType, 'exponential')
def NumPyRandomGeneratorType_exponential(inst, scale=1.0, size=None):
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, scale=1.0, size=None):
            return random_exponential(inst.bit_generator, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, scale=1.0, size=None):
            out = np.empty(size, dtype=np.float64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_exponential(inst.bit_generator, scale)
            return out
        return impl


# Overload the Generator().gamma() method
@overload_method(types.NumPyRandomGeneratorType, 'gamma')
def NumPyRandomGeneratorType_gamma(inst, shape, scale=1.0, size=None):
    check_types(shape, [types.Float, types.Integer, int, float], 'shape')
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, shape, scale=1.0, size=None):
            return random_gamma(inst.bit_generator, shape, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, shape, scale=1.0, size=None):
            out = np.empty(size, dtype=np.float64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_gamma(inst.bit_generator, shape, scale)
            return out
        return impl


# Overload the Generator().beta() method
@overload_method(types.NumPyRandomGeneratorType, 'beta')
def NumPyRandomGeneratorType_beta(inst, a, b, size=None):
    check_types(a, [types.Float, types.Integer, int, float], 'a')
    check_types(b, [types.Float, types.Integer, int, float], 'b')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, a, b, size=None):
            return random_beta(inst.bit_generator, a, b)
        return impl
    else:
        check_size(size)

        def impl(inst, a, b, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_beta(inst.bit_generator, a, b)
            return out
        return impl


# Overload the Generator().f() method
@overload_method(types.NumPyRandomGeneratorType, 'f')
def NumPyRandomGeneratorType_f(inst, dfnum, dfden, size=None):
    check_types(dfnum, [types.Float, types.Integer, int, float], 'dfnum')
    check_types(dfden, [types.Float, types.Integer, int, float], 'dfden')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, dfnum, dfden, size=None):
            return random_f(inst.bit_generator, dfnum, dfden)
        return impl
    else:
        check_size(size)

        def impl(inst, dfnum, dfden, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_f(inst.bit_generator, dfnum, dfden)
            return out
        return impl


# Overload the Generator().chisquare() method
@overload_method(types.NumPyRandomGeneratorType, 'chisquare')
def NumPyRandomGeneratorType_chisquare(inst, df, size=None):
    check_types(df, [types.Float, types.Integer, int, float], 'df')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, df, size=None):
            return random_chisquare(inst.bit_generator, df)
        return impl
    else:
        check_size(size)

        def impl(inst, df, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_chisquare(inst.bit_generator, df)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'standard_cauchy')
def NumPyRandomGeneratorType_standard_cauchy(inst, size=None):

    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, size=None):
            return random_standard_cauchy(inst.bit_generator)
        return impl
    else:
        check_size(size)

        def impl(inst, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_standard_cauchy(inst.bit_generator)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'pareto')
def NumPyRandomGeneratorType_pareto(inst, a, size=None):
    check_types(a, [types.Float, types.Integer, int, float], 'a')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, a, size=None):
            return random_pareto(inst.bit_generator, a)
        return impl
    else:
        check_size(size)

        def impl(inst, a, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_pareto(inst.bit_generator, a)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'weibull')
def NumPyRandomGeneratorType_weibull(inst, a, size=None):
    check_types(a, [types.Float, types.Integer, int, float], 'a')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, a, size=None):
            return random_weibull(inst.bit_generator, a)
        return impl
    else:
        check_size(size)

        def impl(inst, a, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_weibull(inst.bit_generator, a)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'power')
def NumPyRandomGeneratorType_power(inst, a, size=None):
    check_types(a, [types.Float, types.Integer, int, float], 'a')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, a, size=None):
            return random_power(inst.bit_generator, a)
        return impl
    else:
        check_size(size)

        def impl(inst, a, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_power(inst.bit_generator, a)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'laplace')
def NumPyRandomGeneratorType_laplace(inst, loc=0.0, scale=1.0, size=None):
    check_types(loc, [types.Float, types.Integer, int, float], 'loc')
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, loc=0.0, scale=1.0, size=None):
            return random_laplace(inst.bit_generator, loc, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, loc=0.0, scale=1.0, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_laplace(inst.bit_generator, loc, scale)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'logistic')
def NumPyRandomGeneratorType_logistic(inst, loc=0.0, scale=1.0, size=None):
    check_types(loc, [types.Float, types.Integer, int, float], 'loc')
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, loc=0.0, scale=1.0, size=None):
            return random_logistic(inst.bit_generator, loc, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, loc=0.0, scale=1.0, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_logistic(inst.bit_generator, loc, scale)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'lognormal')
def NumPyRandomGeneratorType_lognormal(inst, mean=0.0, sigma=1.0, size=None):
    check_types(mean, [types.Float, types.Integer, int, float], 'mean')
    check_types(sigma, [types.Float, types.Integer, int, float], 'sigma')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, mean=0.0, sigma=1.0, size=None):
            return random_lognormal(inst.bit_generator, mean, sigma)
        return impl
    else:
        check_size(size)

        def impl(inst, mean=0.0, sigma=1.0, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_lognormal(inst.bit_generator, mean, sigma)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'rayleigh')
def NumPyRandomGeneratorType_rayleigh(inst, scale=1.0, size=None):
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, scale=1.0, size=None):
            return random_rayleigh(inst.bit_generator, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, scale=1.0, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_rayleigh(inst.bit_generator, scale)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'standard_t')
def NumPyRandomGeneratorType_standard_t(inst, df, size=None):
    check_types(df, [types.Float, types.Integer, int, float], 'df')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, df, size=None):
            return random_standard_t(inst.bit_generator, df)
        return impl
    else:
        check_size(size)

        def impl(inst, df, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_standard_t(inst.bit_generator, df)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'wald')
def NumPyRandomGeneratorType_wald(inst, mean, scale, size=None):
    check_types(mean, [types.Float, types.Integer, int, float], 'mean')
    check_types(scale, [types.Float, types.Integer, int, float], 'scale')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, mean, scale, size=None):
            return random_wald(inst.bit_generator, mean, scale)
        return impl
    else:
        check_size(size)

        def impl(inst, mean, scale, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_wald(inst.bit_generator, mean, scale)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'geometric')
def NumPyRandomGeneratorType_geometric(inst, p, size=None):
    check_types(p, [types.Float, types.Integer, int, float], 'p')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, p, size=None):
            return np.int64(random_geometric(inst.bit_generator, p))
        return impl
    else:
        check_size(size)

        def impl(inst, p, size=None):
            out = np.empty(size, dtype=np.int64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_geometric(inst.bit_generator, p)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'zipf')
def NumPyRandomGeneratorType_zipf(inst, a, size=None):
    check_types(a, [types.Float, types.Integer, int, float], 'a')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, a, size=None):
            return np.int64(random_zipf(inst.bit_generator, a))
        return impl
    else:
        check_size(size)

        def impl(inst, a, size=None):
            out = np.empty(size, dtype=np.int64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_zipf(inst.bit_generator, a)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'triangular')
def NumPyRandomGeneratorType_triangular(inst, left, mode, right, size=None):
    check_types(left, [types.Float, types.Integer, int, float], 'left')
    check_types(mode, [types.Float, types.Integer, int, float], 'mode')
    check_types(right, [types.Float, types.Integer, int, float], 'right')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, left, mode, right, size=None):
            return random_triangular(inst.bit_generator, left, mode, right)
        return impl
    else:
        check_size(size)

        def impl(inst, left, mode, right, size=None):
            out = np.empty(size)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_triangular(inst.bit_generator,
                                             left, mode, right)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'poisson')
def NumPyRandomGeneratorType_poisson(inst, lam , size=None):
    check_types(lam, [types.Float, types.Integer, int, float], 'lam')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, lam , size=None):
            return np.int64(random_poisson(inst.bit_generator, lam))
        return impl
    else:
        check_size(size)

        def impl(inst, lam , size=None):
            out = np.empty(size, dtype=np.int64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_poisson(inst.bit_generator, lam)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'negative_binomial')
def NumPyRandomGeneratorType_negative_binomial(inst, n, p, size=None):
    check_types(n, [types.Float, types.Integer, int, float], 'n')
    check_types(p, [types.Float, types.Integer, int, float], 'p')
    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst,  n, p , size=None):
            return np.int64(random_negative_binomial(inst.bit_generator, n, p))
        return impl
    else:
        check_size(size)

        def impl(inst, n, p , size=None):
            out = np.empty(size, dtype=np.int64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_negative_binomial(inst.bit_generator, n, p)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'noncentral_chisquare')
def NumPyRandomGeneratorType_noncentral_chisquare(inst, df, nonc, size=None):
    check_types(df, [types.Float, types.Integer, int, float], 'df')
    check_types(nonc, [types.Float, types.Integer, int, float], 'nonc')
    if isinstance(size, types.Omitted):
        size = size.value

    @register_jitable
    def check_arg_bounds(df, nonc):
        if df <= 0:
            raise ValueError("df <= 0")
        if nonc < 0:
            raise ValueError("nonc < 0")

    if is_nonelike(size):
        def impl(inst, df, nonc, size=None):
            check_arg_bounds(df, nonc)
            return np.float64(random_noncentral_chisquare(inst.bit_generator,
                                                          df, nonc))
        return impl
    else:
        check_size(size)

        def impl(inst, df, nonc, size=None):
            check_arg_bounds(df, nonc)
            out = np.empty(size, dtype=np.float64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_noncentral_chisquare(inst.bit_generator,
                                                       df, nonc)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'noncentral_f')
def NumPyRandomGeneratorType_noncentral_f(inst, dfnum, dfden, nonc, size=None):
    check_types(dfnum, [types.Float, types.Integer, int, float], 'dfnum')
    check_types(dfden, [types.Float, types.Integer, int, float], 'dfden')
    check_types(nonc, [types.Float, types.Integer, int, float], 'nonc')
    if isinstance(size, types.Omitted):
        size = size.value

    @register_jitable
    def check_arg_bounds(dfnum, dfden, nonc):
        if dfnum <= 0:
            raise ValueError("dfnum <= 0")
        if dfden <= 0:
            raise ValueError("dfden <= 0")
        if nonc < 0:
            raise ValueError("nonc < 0")

    if is_nonelike(size):
        def impl(inst,  dfnum, dfden, nonc, size=None):
            check_arg_bounds(dfnum, dfden, nonc)
            return np.float64(random_noncentral_f(inst.bit_generator,
                                                  dfnum, dfden, nonc))
        return impl
    else:
        check_size(size)

        def impl(inst, dfnum, dfden, nonc, size=None):
            check_arg_bounds(dfnum, dfden, nonc)
            out = np.empty(size, dtype=np.float64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_noncentral_f(inst.bit_generator,
                                               dfnum, dfden, nonc)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'logseries')
def NumPyRandomGeneratorType_logseries(inst, p, size=None):
    check_types(p, [types.Float, types.Integer, int, float], 'p')
    if isinstance(size, types.Omitted):
        size = size.value

    @register_jitable
    def check_arg_bounds(p):
        if p < 0 or p >= 1 or np.isnan(p):
            raise ValueError("p < 0, p >= 1 or p is NaN")

    if is_nonelike(size):
        def impl(inst, p, size=None):
            check_arg_bounds(p)
            return np.int64(random_logseries(inst.bit_generator, p))
        return impl
    else:
        check_size(size)

        def impl(inst, p, size=None):
            check_arg_bounds(p)
            out = np.empty(size, dtype=np.int64)
            out_f = out.flat
            for i in range(out.size):
                out_f[i] = random_logseries(inst.bit_generator, p)
            return out
        return impl


@overload_method(types.NumPyRandomGeneratorType, 'binomial')
def NumPyRandomGeneratorType_binomial(inst, n, p, size=None):
    check_types(n, [types.Float, types.Integer, int, float], 'n')
    check_types(p, [types.Float, types.Integer, int, float], 'p')

    if isinstance(size, types.Omitted):
        size = size.value

    if is_nonelike(size):
        def impl(inst, n, p, size=None):
            return np.int64(random_binomial(inst.bit_generator, n, p))
        return impl
    else:
        check_size(size)

        def impl(inst, n, p, size=None):
            out = np.empty(size, dtype=np.int64)
            for i in np.ndindex(size):
                out[i] = random_binomial(inst.bit_generator, n, p)
            return out
        return impl
