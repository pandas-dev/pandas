import pytest
import numpy as np

from scipy import stats
from scipy.stats._quantile import (_xp_searchsorted, _estimated_cdf_methods,
    _estimated_cdf_discontinuous_methods, _estimated_cdf_continuous_methods)
from scipy._lib._array_api import (
    xp_default_dtype,
    is_numpy,
    is_torch,
    is_jax,
    is_cupy,
    is_array_api_strict,
    make_xp_test_case,
    SCIPY_ARRAY_API,
    xp_size,
    xp_copy,
)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal
from scipy._lib._util import _apply_over_batch
import scipy._external.array_api_extra as xpx
from scipy.stats._axis_nan_policy import _broadcast_arrays

skip_xp_backends = pytest.mark.skip_xp_backends

lazy_xp_modules = [stats]

@_apply_over_batch(('x', 1), ('p', 1))
def quantile_reference_last_axis(x, p, nan_policy, method):
    if nan_policy == 'omit':
        x = x[~np.isnan(x)]
    p_mask = np.isnan(p)
    p = p.copy()
    p[p_mask] = 0.5
    if method == 'harrell-davis':
        # hdquantiles returns masked element if length along axis is 1 (bug)
        res = (np.full_like(p, x[0]) if x.size == 1
               else stats.mstats.hdquantiles(x, p).data)
    elif method.startswith('round'):
        res = winsor_reference_1d(np.sort(x), p, method)
    else:
        res = np.quantile(x, p, method=method)

    res = np.asarray(res)
    if nan_policy == 'propagate' and np.any(np.isnan(x)):
        res[:] = np.nan

    res[p_mask] = np.nan
    return res


@np.vectorize(excluded={0, 2})  # type: ignore[call-arg]
def winsor_reference_1d(y, p, method):
    # Adapted directly from the documentation
    # Note: `y` is the sorted data array
    n = len(y)
    if method == 'round_nearest':
        j = int(np.round(p * n) if p < 0.5 else np.round(n * p - 1))
    elif method == 'round_outward':
        j = int(np.floor(p * n) if p < 0.5 else np.ceil(n * p - 1))
    elif method == 'round_inward':
        j = int(np.ceil(p * n) if p < 0.5 else np.floor(n * p - 1))
    return y[j]


def quantile_reference(x, p, *, axis, nan_policy, keepdims, method):
    x, p = np.moveaxis(x, axis, -1), np.moveaxis(np.atleast_1d(p), axis, -1)
    res = quantile_reference_last_axis(x, p, nan_policy, method)
    res = np.moveaxis(res, -1, axis)
    if not keepdims:
        res = np.squeeze(res, axis=axis)
    return res


@make_xp_test_case(stats.quantile)
class TestQuantile:

    def test_input_validation(self, xp):
        x = xp.asarray([1, 2, 3])
        p = xp.asarray(0.5)

        message = "`x` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(xp.asarray([True, False]), p)
        with pytest.raises(ValueError):
            stats.quantile(xp.asarray([1+1j, 2]), p)

        message = "`p` must have real floating dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, xp.asarray([0, 1]))

        message = "`weights` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, weights=xp.astype(x, xp.complex64))

        message = "`axis` must be an integer or None."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=0.5)
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=(0, -1))

        message = "`axis` is not compatible with the shapes of the inputs."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=2)

        if not is_jax(xp):  # no data-dependent input validation for lazy arrays
            message = "The input contains nan values"
            with pytest.raises(ValueError, match=message):
                stats.quantile(xp.asarray([xp.nan, 1, 2]), p, nan_policy='raise')

        message = "`method` must be one of..."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, method='a duck')

        message = "`method='harrell-davis'` does not support `weights`."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, weights=x, method='harrell-davis')

        message = "`method='round_nearest'` does not support `weights`."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, weights=x, method='round_nearest')

        message = "If specified, `keepdims` must be True or False."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, keepdims=42)

        message = "`keepdims` may be False only if the length of `p` along `axis` is 1."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, xp.asarray([0.5, 0.6]), keepdims=False)


    def _get_weights_x_rep(self, x, axis, rng):
        x = np.swapaxes(x, axis, -1)
        ndim = x.ndim
        x = np.atleast_2d(x)
        counts = rng.integers(10, size=x.shape[-1], dtype=np.int32)
        x_rep = []
        weights = []
        for x_ in x:
            counts_ = rng.permuted(counts)
            x_rep.append(np.repeat(x_, counts_))
            weights.append(counts_)
        x_rep, weights = np.stack(x_rep), np.stack(weights)
        if ndim < 2:
            x_rep, weights = np.squeeze(x_rep, axis=0), np.squeeze(weights, axis=0)
        x_rep, weights = np.swapaxes(x_rep, -1, axis), np.swapaxes(weights, -1, axis)
        weights = np.asarray(weights, dtype=x.dtype)
        return weights, x_rep


    @skip_xp_backends(cpu_only=True, reason="PyTorch doesn't have `betainc`.",
                      exceptions=['cupy'])
    @pytest.mark.parametrize('method',
         ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
          'hazen', 'interpolated_inverted_cdf', 'linear',
          'median_unbiased', 'normal_unbiased', 'weibull',
          'harrell-davis', 'round_nearest', 'round_outward', 'round_inward',
          '_lower', '_higher', '_midpoint', '_nearest'])
    @pytest.mark.parametrize('shape_x, shape_p, axis',
        [(10, None, -1), (10, 10, -1), (10, (2, 3), -1), ((10, 2), None, 0)])
    @pytest.mark.parametrize('weights', [False, True])
    def test_against_reference(self, method, shape_x, shape_p, axis, weights, xp):
        # Test all methods with various data shapes
        if weights and (method.startswith('_') or method.startswith('round')
                        or method=='harrell-davis'):
            pytest.skip('`weights` not supported by private (legacy) methods.')
        dtype = xp_default_dtype(xp)
        rng = np.random.default_rng(23458924568734956)
        x = rng.random(size=shape_x)
        p = rng.random(size=shape_p)

        if weights:
            weights, x_rep = self._get_weights_x_rep(x, axis, rng)
        else:
            weights, x_rep = None, x

        ref = quantile_reference(
            x_rep, p, method=method[1:] if method.startswith('_') else method,
            axis=axis, nan_policy='propagate', keepdims=shape_p is not None)

        x, p = xp.asarray(x, dtype=dtype), xp.asarray(p, dtype=dtype)
        weights = weights if weights is None else xp.asarray(weights, dtype=dtype)
        res = stats.quantile(x, p, method=method, weights=weights, axis=axis)

        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @pytest.mark.filterwarnings("ignore:torch.searchsorted:UserWarning")
    @skip_xp_backends(cpu_only=True, reason="PyTorch doesn't have `betainc`.",
                      exceptions=['cupy', 'jax.numpy'])
    @pytest.mark.parametrize('axis', [0, 1])
    @pytest.mark.parametrize('keepdims', [False, True])
    @pytest.mark.parametrize('nan_policy', ['omit', 'propagate', 'marray'])
    @pytest.mark.parametrize('dtype', ['float32', 'float64'])
    @pytest.mark.parametrize('method', ['linear', 'harrell-davis', 'round_nearest'])
    @pytest.mark.parametrize('weights', [False, True])
    def test_against_reference_2(self, axis, keepdims, nan_policy,
                                 dtype, method, weights, xp):
        # Test some methods with various combinations of arguments
        if is_jax(xp) and nan_policy == 'marray':  # mdhaber/marray#146
            pytest.skip("`marray` currently incompatible with JAX")
        if weights and method in {'harrell-davis', 'round_nearest'}:
            pytest.skip("These methods don't yet support weights")
        rng = np.random.default_rng(23458924568734956)
        shape = (5, 6)
        x = rng.random(size=shape).astype(dtype)
        p = rng.random(size=shape).astype(dtype)
        mask = rng.random(size=shape) > 0.8
        assert np.any(mask)
        x[mask] = np.nan
        if not keepdims:
            p = np.mean(p, axis=axis, keepdims=True)

        # inject p = 0 and p = 1 to test edge cases
        # Currently would fail with CuPy/JAX (cupy/cupy#8934, jax-ml/jax#21900);
        # remove the `if` when those are resolved.
        if is_numpy(xp):
            p0 = p.ravel()
            p0[1] = 0.
            p0[-2] = 1.

        dtype = getattr(xp, dtype)

        if weights:
            weights, x_rep = self._get_weights_x_rep(x, axis, rng)
            weights = weights if weights is None else xp.asarray(weights)
        else:
            weights, x_rep = None, x

        if nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            if weights is not None:
                pytest.skip("MArray is not yet compatible with weights")
            marray = pytest.importorskip('marray')
            kwargs = dict(axis=axis, keepdims=keepdims, method=method)
            mxp = marray._get_namespace(xp)
            x_mp = mxp.asarray(x, mask=mask)
            weights = weights if weights is None else mxp.asarray(weights)
            res = stats.quantile(x_mp, mxp.asarray(p), weights=weights, **kwargs)
            ref = quantile_reference(x_rep, p, nan_policy='omit', **kwargs)
            xp_assert_close(res.data, xp.asarray(ref, dtype=dtype))
            return

        kwargs = dict(axis=axis, keepdims=keepdims,
                      nan_policy=nan_policy, method=method)
        res = stats.quantile(xp.asarray(x), xp.asarray(p), weights=weights, **kwargs)
        ref = quantile_reference(x_rep, p, **kwargs)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    def test_integer_input_output_dtype(self, xp):
        res = stats.quantile(xp.arange(10, dtype=xp.int64), 0.5)
        assert res.dtype == xp_default_dtype(xp)

    @pytest.mark.parametrize('x, p, ref, kwargs',
        [([], 0.5, np.nan, {}),
         ([1, 2, 3], [-1, 0, 1, 1.5, np.nan], [np.nan, 1, 3, np.nan, np.nan], {}),
         ([1, 2, 3], [], [], {}),
         ([[np.nan, 2]], 0.5, [np.nan, 2], {'nan_policy': 'omit'}),
         ([[], []], 0.5, np.full(2, np.nan), {'axis': -1}),
         ([[], []], 0.5, np.zeros((0,)), {'axis': 0, 'keepdims': False}),
         ([[], []], 0.5, np.zeros((1, 0)), {'axis': 0, 'keepdims': True}),
         ([], [0.5, 0.6], np.full(2, np.nan), {}),
         (np.arange(1, 28).reshape((3, 3, 3)), 0.5, [[[14.]]],
          {'axis': None, 'keepdims': True}),
         ([[1, 2], [3, 4]], [0.25, 0.5, 0.75], [[1.75, 2.5, 3.25]],
          {'axis': None, 'keepdims': True}),
         # Known issue:
         # ([1, 2, 3], 0.5, 2., {'weights': [0, 0, 0]})
         # See https://github.com/scipy/scipy/pull/23941#issuecomment-3503554361
         ])
    def test_edge_cases(self, x, p, ref, kwargs, xp):
        default_dtype = xp_default_dtype(xp)
        x, p, ref = xp.asarray(x), xp.asarray(p), xp.asarray(ref, dtype=default_dtype)
        res = stats.quantile(x, p, **kwargs)
        xp_assert_equal(res, ref)

    @pytest.mark.parametrize('axis', [0, 1, 2])
    @pytest.mark.parametrize('keepdims', [False, True])
    def test_size_0(self, axis, keepdims, xp):
        shape = [3, 4, 0]
        out_shape = shape.copy()
        if keepdims:
            out_shape[axis] = 1
        else:
            out_shape.pop(axis)
        res = stats.quantile(xp.zeros(tuple(shape)), 0.5, axis=axis, keepdims=keepdims)
        assert res.shape == tuple(out_shape)

    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
         '_lower', '_higher', '_midpoint', '_nearest'])
    def test_transition(self, method, xp):
        # test that values of discontinuous estimators are correct when
        # p*n + m - 1 is integral.
        if method == 'closest_observation' and np.__version__ < '2.0.1':
            pytest.skip('Bug in np.quantile (numpy/numpy#26656) fixed in 2.0.1')
        x = np.arange(8., dtype=np.float64)
        p = np.arange(0, 1.03125, 0.03125)
        res = stats.quantile(xp.asarray(x), xp.asarray(p), method=method)
        ref = np.quantile(x, p, method=method[1:] if method.startswith('_') else method)
        xp_assert_equal(res, xp.asarray(ref, dtype=xp.float64))

    @pytest.mark.parametrize('zero_weights', [False, True])
    def test_weights_against_numpy(self, zero_weights, xp):
        if is_numpy(xp) and xp.__version__ < "2.1.3" and zero_weights:
            pytest.skip('`Bug in np.quantile (numpy/numpy#27563) fixed in 2.1.3')
        dtype = xp_default_dtype(xp)
        rng = np.random.default_rng(85468924398205602)
        method = 'inverted_cdf'
        x = rng.random(size=100)
        weights = rng.random(size=100)
        if zero_weights:
            weights[weights < 0.5] = 0
        p = np.linspace(0., 1., 300)
        res = stats.quantile(xp.asarray(x, dtype=dtype), xp.asarray(p, dtype=dtype),
                             method=method, weights=xp.asarray(weights, dtype=dtype))
        ref = np.quantile(x, p, method=method, weights=weights)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation', 'hazen',
         'interpolated_inverted_cdf', 'linear','median_unbiased', 'normal_unbiased',
         'weibull'])
    def test_zero_weights(self, method, xp):
        rng = np.random.default_rng(85468924398205602)

        # test 1-D versus eliminating zero-weighted values
        n = 100
        x = xp.asarray(rng.random(size=n))
        x0 = xp_copy(x)
        p = xp.asarray(rng.random(size=n))
        i_zero = xp.asarray(rng.random(size=n) < 0.1)
        weights = xp.asarray(rng.random(size=n))
        weights = xp.where(i_zero, 0., weights)
        res = stats.quantile(x, p, weights=weights, method=method)
        ref = stats.quantile(x[~i_zero], p, weights=weights[~i_zero], method=method)
        xp_assert_close(res, ref)
        xp_assert_equal(x, x0)  # no input mutation

        # test multi-D versus `nan_policy='omit'`
        shape = (5, 100)
        x = xp.asarray(rng.random(size=shape))
        x0 = xp_copy(x)
        p = xp.asarray(rng.random(size=shape))
        i_zero = xp.asarray(rng.random(size=shape) < 0.1)
        weights = xp.asarray(rng.random(size=shape))
        x_nanned = xp.where(i_zero, xp.nan, x)
        weights_zeroed = xp.where(i_zero, 0., weights)
        res = stats.quantile(x, p, weights=weights_zeroed, method=method, axis=-1)
        ref = stats.quantile(x_nanned, p, weights=weights,
                             nan_policy='omit', method=method, axis=-1)
        xp_assert_close(res, ref)
        xp_assert_equal(x, x0)  # no input mutation

    @pytest.mark.filterwarnings("ignore:torch.searchsorted:UserWarning")
    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation', 'hazen',
         'interpolated_inverted_cdf', 'linear','median_unbiased', 'normal_unbiased',
         'weibull'])
    @pytest.mark.parametrize('shape', [50, (50, 3)])
    def test_unity_weights(self, method, shape, xp):
        # Check that result is unchanged if all weights are `1.0`
        rng = np.random.default_rng(28546892439820560)
        x = xp.asarray(rng.random(size=shape))
        p = xp.asarray(rng.random(size=shape))
        weights = xp.ones_like(x)
        res = stats.quantile(x, p, weights=weights, method=method)
        ref = stats.quantile(x, p, method=method)
        xp_assert_close(res, ref)

    @skip_xp_backends(cpu_only=True, reason="PyTorch doesn't have `betainc`.",
                      exceptions=['cupy', 'jax.numpy'])
    def test_all_nan_harrell_davis_gh24707(self, xp):
        # While working on gh-24707, there was a case in which if *all* elements of one
        # slice were NaN, only some elements of another slice were NaN, and
        # `nan_policy='omit'`, then the NaNs would not be ignored in the other slice.
        # The same test case could fail with `nan_policy='propagate'` for a different
        # reason, if the fix were not made carefully. Check that both these cases are
        # resolved.
        kwargs = dict(method='harrell-davis', axis=-1)
        x = xp.asarray([[xp.nan, xp.nan, xp.nan], [xp.nan, 2, 3]])

        res = stats.quantile(x, 0.5, **kwargs, nan_policy='omit')
        xp_assert_close(res, xp.asarray([xp.nan, 2.5]))

        res = stats.quantile(x, 0.5, **kwargs, nan_policy='propagate')
        xp_assert_close(res, xp.asarray([xp.nan, xp.nan]))


@_apply_over_batch(('a', 1), ('v', 1))
def np_searchsorted(a, v, side):
    return np.searchsorted(a, v, side=side)


@make_xp_test_case(_xp_searchsorted)
class Test_XPSearchsorted:
    @pytest.mark.parametrize('side', ['left', 'right'])
    @pytest.mark.parametrize('ties', [False, True])
    @pytest.mark.parametrize('shape', [0, 1, 2, 10, 11, 1000, 10001,
                                       (2, 0), (0, 2), (2, 10), (2, 3, 11)])
    @pytest.mark.parametrize('nans_x', [False, True])
    @pytest.mark.parametrize('infs_x', [False, True])
    def test_nd(self, side, ties, shape, nans_x, infs_x, xp):
        if nans_x and is_torch(xp):
            pytest.skip('torch sorts NaNs differently')
        rng = np.random.default_rng(945298725498274853)
        if ties:
            x = rng.integers(5, size=shape)
        else:
            x = rng.random(shape)
        # float32 is to accommodate JAX - nextafter with `float64` is too small?
        x = np.asarray(x, dtype=np.float32)
        xr = np.nextafter(x, np.inf)
        xl = np.nextafter(x, -np.inf)
        x_ = np.asarray([-np.inf, np.inf, np.nan])
        x_ = np.broadcast_to(x_, x.shape[:-1] + (3,))
        y = rng.permuted(np.concatenate((xl, x, xr, x_), axis=-1), axis=-1)
        if nans_x:
            mask = rng.random(shape) < 0.1
            x[mask] = np.nan
        if infs_x:
            mask = rng.random(shape) < 0.1
            x[mask] = -np.inf
            mask = rng.random(shape) > 0.9
            x[mask] = np.inf
        x = np.sort(x, axis=-1)
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        xp_default_int = xp.asarray(1).dtype
        if xp_size(x) == 0 and x.ndim > 0 and x.shape[-1] != 0:
            ref = xp.empty(x.shape[:-1] + (y.shape[-1],), dtype=xp_default_int)
        else:
            ref = xp.asarray(np_searchsorted(x, y, side=side), dtype=xp_default_int)
        x, y = xp.asarray(x), xp.asarray(y)
        res = _xp_searchsorted(x, y, side=side)
        xp_assert_equal(res, ref)


@_apply_over_batch(('x', 1), ('y', 1))
def estimated_cdf_reference_last_axis(x, y, nan_policy, method):
    i_nan = np.isnan(x)
    if nan_policy == 'propagate' and np.any(i_nan):
        return np.full_like(y, np.nan)
    elif nan_policy == 'omit':
        x = x[~i_nan]
    return stats.estimated_cdf(x, y, keepdims=True, method=method)


def estimated_cdf_reference(x, y, *, axis=0, nan_policy='propagate',
                        keepdims=None, method='linear'):
    x, y = _broadcast_arrays((x, y), axis=axis)
    x, y = np.moveaxis(x, axis, -1), np.moveaxis(y, axis, -1)
    res = estimated_cdf_reference_last_axis(x, y, nan_policy, method)
    res = np.moveaxis(res, -1, axis)
    if not keepdims:
        res = np.squeeze(res, axis=axis)
    return res


# avoid variable collection order issues
_estimated_cdf_methods_list = sorted(list(_estimated_cdf_methods))


@make_xp_test_case(stats.estimated_cdf)
class TestEstimatedCDF:
    def test_input_validation(self, xp):
        x = xp.asarray([1, 2, 3])
        y = xp.asarray(2)

        message = "`x` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(xp.asarray([True, False]), y)
        with pytest.raises(ValueError):
            stats.estimated_cdf(xp.asarray([1+1j, 2]), y)

        message = "`y` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, xp.asarray([0+1j, 1]))

        message = "`axis` must be an integer or None."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, y, axis=0.5)
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, y, axis=(0, -1))

        message = "`axis` is not compatible with the shapes of the inputs."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, y, axis=2)

        if not is_jax(xp):  # no data-dependent input validation for lazy arrays
            message = "The input contains nan values"
            with pytest.raises(ValueError, match=message):
                stats.estimated_cdf(xp.asarray([xp.nan, 1, 2]), y, nan_policy='raise')

        message = "method` must be one of..."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, y, method='a duck')

        message = "If specified, `keepdims` must be True or False."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, y, keepdims=42)

        message = "`keepdims` may be False only if the length of `y` along `axis` is 1."
        with pytest.raises(ValueError, match=message):
            stats.estimated_cdf(x, xp.asarray([0.5, 0.6]), keepdims=False)

    @pytest.mark.parametrize('method', _estimated_cdf_methods_list)
    @pytest.mark.parametrize('x_shape', [2, 10, 11, 100, 1001, (2, 10), (2, 3, 11)])
    @pytest.mark.parametrize('y_shape', [None, 25])
    @pytest.mark.parametrize('ties', [False, True])
    def test_against_quantile(self, method, x_shape, y_shape, ties, xp):
        discontinuous = method in _estimated_cdf_discontinuous_methods
        dtype = xp_default_dtype(xp)  # removed parameterization to speed up tests
        rng = np.random.default_rng(394529872549827485)
        y_shape = x_shape if y_shape is None else y_shape

        if ties:
            x = xp.asarray(rng.integers(9, size=x_shape), dtype=dtype)
        else:
            x = xp.asarray(rng.standard_normal(size=x_shape), dtype=dtype)

        p = xp.asarray(rng.random(size=y_shape), dtype=dtype)
        y = stats.quantile(x, p, method=method, axis=-1)
        res = stats.estimated_cdf(x, y, method=method, axis=-1)
        ref = xp.broadcast_to(p, (*x.shape[:-1], y.shape[-1]))

        # check that `quantile` is the inverse of `estimated_cdf`
        # note that for discontinuous methods, res is right on the cusp of a transition,
        # and there can be a tiny bit of error to the right or left. We shift it left
        # to ensure we're on the correct side of the transition, producing the same `y2`
        # as if the probability calculation were exact.
        res = res - 1e-6 if discontinuous else res
        y2 = stats.quantile(x, res, method=method, axis=-1)
        atol = 1e-5 if dtype == xp.float32 else 1e-12
        xp_assert_close(y2, y, atol=atol)

        # if there are ties or method is discontinuous, `quantile` is not invertible
        if ties or discontinuous:
            return

        # `quantile` is not invertible outside this domain
        a, b = _estimated_cdf_continuous_methods[method]
        n = x.shape[-1]
        p_min = (1 - a) / (n + 1 - a - b)
        p_max = (n - a) / (n + 1 - a - b)
        i_very_low = y < xp.min(x, axis=-1, keepdims=True)
        i_very_high = y > xp.max(x, axis=-1, keepdims=True)
        i_low = (ref <= p_min) & ~i_very_low
        i_high = (ref >= p_max) & ~i_very_high
        i_ok = ~(i_low | i_high | i_very_low | i_very_high)

        # check for correct inversion within the domain
        xp_assert_close(res[i_ok], ref[i_ok])

        # check that all other values get mapped to bottom or top of range
        kwargs = dict(check_shape=False, check_dtype=False, check_0d=True)
        xp_assert_close(res[i_low], xp.asarray(p_min), **kwargs)
        xp_assert_close(res[i_high], xp.asarray(p_max), **kwargs)
        xp_assert_close(res[i_very_low], xp.asarray(0.0), **kwargs)
        xp_assert_close(res[i_very_high], xp.asarray(1.0), **kwargs)

    @pytest.mark.filterwarnings("ignore:torch.searchsorted:UserWarning")
    @pytest.mark.parametrize('axis', [0, 1])
    @pytest.mark.parametrize('keepdims', [False, True])
    @pytest.mark.parametrize('nan_policy', ['propagate', 'omit', 'marray'])
    @pytest.mark.parametrize('dtype', ['float32', 'float64'])
    @pytest.mark.parametrize('nans', [False, True])
    @pytest.mark.parametrize('meth', ['linear', 'inverted_cdf'])
    def test_against_reference(self, axis, keepdims, nan_policy, dtype, nans, meth, xp):
        if is_jax(xp) and nan_policy == 'marray':  # mdhaber/marray#146
            pytest.skip("`marray` currently incompatible with JAX")
        rng = np.random.default_rng(23458924568734956)
        shape = (5, 6)
        x = rng.standard_normal(size=shape).astype(dtype)
        y = rng.standard_normal(size=shape).astype(dtype)

        mask = None
        if nans:
            mask = rng.random(size=shape) > 0.8
            assert np.any(mask)
            x[mask] = np.nan

        if not keepdims:
            y = np.mean(y, axis=axis, keepdims=True)

        dtype = getattr(xp, dtype)

        if nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            marray = pytest.importorskip('marray')
            kwargs = dict(axis=axis, keepdims=keepdims, method=meth)
            mxp = marray._get_namespace(xp)
            x_mp = mxp.asarray(x, mask=mask)
            res = stats.estimated_cdf(x_mp, mxp.asarray(y), **kwargs)
            ref = estimated_cdf_reference(x, y, nan_policy='omit', **kwargs)
            xp_assert_close(res.data, xp.asarray(ref, dtype=dtype))
            return

        kwargs = dict(axis=axis, keepdims=keepdims,
                      nan_policy=nan_policy, method=meth)
        res = stats.estimated_cdf(xp.asarray(x), xp.asarray(y), **kwargs)
        ref = estimated_cdf_reference(x, y, **kwargs)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @pytest.mark.skip_xp_backends('torch', reason='issues with sorting NaNs')
    @pytest.mark.parametrize('n', [50, 500])
    @pytest.mark.parametrize('method, ab', _estimated_cdf_continuous_methods.items())
    def test_plotting_positions(self, n, method, ab, xp):
        a, b = ab
        rng = np.random.default_rng(539452987254982748)
        x = rng.standard_normal(n)

        mask = rng.random(n) < 0.1
        x[mask] = np.nan
        mask = xp.asarray(mask)

        x_masked = np.ma.masked_invalid(x)
        ref = stats.mstats.plotting_positions(x_masked, a, b)
        ref = xp.asarray(ref.data)

        x = xp.asarray(x)
        res = stats.estimated_cdf(x, x, nan_policy='omit', method=method)

        xp_assert_close(res[~mask], ref[~mask])
        assert xp.all(xp.isnan(res[mask]))

    @pytest.mark.parametrize('ties', [False, True])
    def test_against_ecdf_percentileofscore(self, ties, xp):
        rng = np.random.default_rng(853945298725498274)
        n = 50
        dtype = xp_default_dtype(xp)
        x = rng.integers(10, size=n) if ties else rng.standard_normal(size=n)
        y = rng.integers(10, size=25) if ties else rng.standard_normal(size=25)
        ref = stats.ecdf(x).cdf.evaluate(y)
        ref2 = stats.percentileofscore(x, y, 'weak')
        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.estimated_cdf(x, y, method='inverted_cdf')
        ref, ref2 = xp.asarray(ref, dtype=dtype), xp.asarray(ref2, dtype=dtype)
        xp_assert_close(res, ref)
        xp_assert_close(res, ref2 / 100)

    def test_integer_input_output_dtype(self, xp):
        x = xp.arange(10, dtype=xp.int64)
        res = stats.estimated_cdf(x, x)
        assert res.dtype == xp_default_dtype(xp)

    @pytest.mark.parametrize('nan_policy', ['propagate', 'omit', 'marray'])
    @pytest.mark.parametrize('method', _estimated_cdf_methods_list)
    def test_size_one_sample(self, nan_policy, method, xp):
        discontinuous = method in _estimated_cdf_discontinuous_methods
        x = xp.arange(10.)
        y = xp.asarray([0., -1., 1.])
        n = np.asarray(1.) if is_array_api_strict(xp) else xp.asarray(1.)
        with np.errstate(divide='ignore', invalid='ignore'):  # for method = 'linear'
            if discontinuous:
                ref = xp.asarray([1., 0., 1.])
            else:
                a, b = _estimated_cdf_continuous_methods[method]
                ref = xp.asarray([float((n - a) / (n + 1 - a - b)), 0., 1.])

        if nan_policy == 'propagate':
            x = x[:1]
            kwargs = {'nan_policy': 'propagate'}
        elif nan_policy == 'omit':
            x = xpx.at(x)[1:].set(xp.nan)
            kwargs = {'nan_policy': 'omit'}
        elif nan_policy == 'marray':
            if is_jax(xp):
                pytest.skip("JAX currently incompatible with `marray`")
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            marray = pytest.importorskip('marray')
            mxp = marray._get_namespace(xp)
            mask = (x > 0.)
            x = mxp.asarray(x, mask=mask)
            y = mxp.asarray(y)
            kwargs = {}

        with np.errstate(divide='ignore', invalid='ignore'):  # for method = 'linear'
            res = stats.estimated_cdf(x, y, method=method, **kwargs)
        res = res.data if nan_policy == 'marray' else res
        xp_assert_close(res, ref)

    # skipping marray due to mdhaber/marray#24
    @pytest.mark.parametrize('nan_policy', ['propagate', 'omit'])
    @pytest.mark.parametrize('method', _estimated_cdf_methods_list)
    def test_size_zero_sample(self, nan_policy, method, xp):
        x = xp.arange(10.)
        y = xp.asarray([0., -1., 1.])  # this should work
        ref = xp.full_like(y, xp.nan)

        if nan_policy == 'propagate':
            x = x[0:0]
            kwargs = {'nan_policy': 'propagate'}
        elif nan_policy == 'omit':
            x = xpx.at(x)[:].set(xp.nan)
            kwargs = {'nan_policy': 'omit'}
        elif nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            if is_jax(xp):
                pytest.skip("JAX currently incompatible with `marray`")
            marray = pytest.importorskip('marray')
            mxp = marray._get_namespace(xp)
            mask = (x >= 0.)
            x = mxp.asarray(x, mask=mask)
            y = mxp.asarray(y)
            kwargs = {}

        with np.errstate(divide='ignore', invalid='ignore'):  # for method = 'linear'
            res = stats.estimated_cdf(x, y, method=method, **kwargs)

        if nan_policy == 'marray':
            assert xp.all(res.mask)
        else:
            xp_assert_close(res, ref)

    @pytest.mark.parametrize('x, y, ref, kwargs',
        [
         ([], 0.5, np.nan, {}),
         ([1, 2, 3], [0.999, 3.001, np.nan], [0., 1., np.nan], {}),
         ([1, 2, 3], [], [], {}),
         ([[np.nan, 2]], 2, [np.nan, 0.5], {'nan_policy': 'omit', 'method': 'weibull'}),
         ([[], []], 0.5, np.full(2, np.nan), {'axis': -1}),
         ([[], []], 0.5, np.zeros((0,)), {'axis': 0, 'keepdims': False}),
         ([[], []], 0.5, np.zeros((1, 0)), {'axis': 0, 'keepdims': True}),
         ([], [0.5, 0.6], np.full(2, np.nan), {}),
         (np.arange(1, 28).reshape((3, 3, 3)), 14., [[[0.5]]],
          {'axis': None, 'keepdims': True}),
         ([[1, 2], [3, 4]], [1.75, 2.5, 3.25], [[0.25, 0.5, 0.75]],
          {'axis': None, 'keepdims': True}),
         ([1, 2, 3], [-np.inf, np.inf], [0.0, 1.0], {}),
         # It is our choice how much effort and computational overhead we want to put
         # into adjusting for insane edge cases like when `x` contains infinite values,
         # especially when `y` does, too.
         # One practical argument would be that y = +/- inf should produce the same
         # results as an extremely large finite number, in which case the 0th element
         # of the 'linear' result should be `0`.
         # Another argument would be for producing NaN below wherever y = +/- inf.
         # Another would be that it's not appropriate to spend significant computation
         # correcting these edge cases; we should just document what we do.
         # In any case, this is the current status.
         ([-np.inf, -1, 0, 1, np.inf], [-np.inf, -2, -1, 0, 1, 2, np.inf],
          [0.25, 0.25, 0.25, 0.5 , 0.75, 0.75, np.nan], {'method':'linear'}),
         ([-np.inf, -1, 0, 1, np.inf], [-np.inf, -2, -1, 0, 1, 2, np.inf],
          [0.2, 0.2, 0.4, 0.6, 0.8, 0.8, 1.], {'method': 'inverted_cdf'}),
        ])
    def test_edge_cases(self, x, y, ref, kwargs, xp):
        if kwargs.get('method', None) == 'weibull' and is_torch(xp):
            pytest.skip('data-apis/array-api-compat#360')
        if kwargs.get('axis', None) == -1 and is_cupy(xp):
            pytest.skip('Fails; need to investigate.')
        default_dtype = xp_default_dtype(xp)
        x, y, ref = xp.asarray(x), xp.asarray(y), xp.asarray(ref, dtype=default_dtype)
        res = stats.estimated_cdf(x, y, **kwargs)
        xp_assert_close(res, ref, rtol=1e-15)

    @pytest.mark.skip_xp_backends('jax.numpy', reason="-1e-45 is not less than 0?")
    @pytest.mark.parametrize('method', _estimated_cdf_discontinuous_methods.keys())
    def test_transition(self, method, xp):
        # test that values of discontinuous estimators are as expected around
        # transition point
        x = np.arange(8., dtype=np.float64)
        xl, xr = np.nextafter(x, -np.inf), np.nextafter(x, np.inf)
        offset = 0.5 if method == 'closest_observation' else 0.0
        ref_r = (x + 1 + offset) / 8
        ref_r[-1] = 1.0  # value is greater than or equal to the maximum observation
        ref_l = (x + offset) / 8
        ref_l[0] = 0.0  # value is less than the minimum observation
        x, xl, xr = xp.asarray(x), xp.asarray(xl), xp.asarray(xr)
        ref_l, ref_r = xp.asarray(ref_l), xp.asarray(ref_r)
        xp_assert_equal(stats.estimated_cdf(x, x, method=method), ref_r)
        xp_assert_equal(stats.estimated_cdf(x, xr, method=method), ref_r)
        xp_assert_equal(stats.estimated_cdf(x, xl, method=method), ref_l)
