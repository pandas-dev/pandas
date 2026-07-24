import pytest
import numpy as np
from numpy.testing import assert_allclose

from scipy.conftest import skip_xp_invalid_arg
import scipy._external.array_api_extra as xpx
from scipy._lib._array_api import (make_xp_test_case, xp_default_dtype, is_jax, is_cupy,
                                   eager_warns, xp_result_type, is_array_api_strict)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal
from scipy import stats
from scipy.stats._axis_nan_policy import SmallSampleWarning


lazy_xp_modules = [stats]


@make_xp_test_case(stats.chatterjeexi)
class TestChatterjeeXi:
    @pytest.mark.parametrize('case', [
        dict(y_cont=True, statistic=-0.303030303030303, pvalue=0.9351329808526656),
        dict(y_cont=False, statistic=0.07407407407407396, pvalue=0.3709859367123997)])
    @pytest.mark.parametrize('dtype', ['float32', 'float64', None])
    def test_against_R_XICOR(self, case, dtype, xp):
        # Test against R package XICOR, e.g.
        # library(XICOR)
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558,
        #       0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495)
        # y = c(0.8402081904493103, 0.5946972833914318, 0.23481606164114155,
        #       0.49754786197715384, 0.9146460831206026, 0.5848057749217579,
        #       0.7620801065573549, 0.31410063302647495, 0.7935620302236199,
        #       0.5423085761365468)
        # xicor(x, y, ties=FALSE, pvalue=TRUE)
        dtype = xp_default_dtype(xp) if dtype is None else getattr(xp, dtype)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=10)
        y = (rng.random(size=10) if case['y_cont']
             else rng.integers(0, 5, size=10))

        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.chatterjeexi(x, y, y_continuous=case['y_cont'])

        xp_assert_close(res.statistic, xp.asarray(case['statistic'], dtype=dtype))
        xp_assert_close(res.pvalue, xp.asarray(case['pvalue'], dtype=dtype))

    @pytest.mark.parametrize('y_continuous', (False, True))
    def test_permutation_asymptotic(self, y_continuous):
        # XICOR doesn't seem to perform the permutation test as advertised, so
        # compare the result of a permutation test against an asymptotic test.
        rng = np.random.default_rng(2524579827426)
        n = np.floor(rng.uniform(100, 150)).astype(int)
        shape = (2, n)
        x = rng.random(size=shape)
        y = (rng.random(size=shape) if y_continuous
             else rng.integers(0, 10, size=shape))
        method = stats.PermutationMethod(rng=rng)
        res = stats.chatterjeexi(x, y, method=method,
                                 y_continuous=y_continuous, axis=-1)
        ref = stats.chatterjeexi(x, y, y_continuous=y_continuous, axis=-1)
        np.testing.assert_allclose(res.statistic, ref.statistic, rtol=1e-15)
        np.testing.assert_allclose(res.pvalue, ref.pvalue, rtol=2e-2)

    def test_input_validation(self, xp):
        rng = np.random.default_rng(25932435798274926)
        x, y = rng.random(size=(2, 10))
        x, y = xp.asarray(x), xp.asarray(y)

        message = 'Array shapes are incompatible for broadcasting.|Incompatible shapes'
        with pytest.raises((ValueError, TypeError), match=message):
            stats.chatterjeexi(x, y[:-1])

        if not (is_jax(xp) or is_cupy(xp)):
            # jax misses out on some input validation from _axis_nan_policy decorator
            # This also fails with cupy for reasons that need to be investigated.
            message = '...axis 10 is out of bounds for array...|out of range'
            with pytest.raises((ValueError, IndexError), match=message):
                stats.chatterjeexi(x, y, axis=10)

        message = '`y_continuous` must be boolean.'
        with pytest.raises(ValueError, match=message):
            stats.chatterjeexi(x, y, y_continuous='a herring')

        message = "`method` must be 'asymptotic' or"
        with pytest.raises(ValueError, match=message):
            stats.chatterjeexi(x, y, method='ekki ekii')

    def test_special_cases(self, xp):
        message = 'One or more sample arguments is too small...'
        with pytest.warns(SmallSampleWarning, match=message):
            res = stats.chatterjeexi(xp.asarray([1]), xp.asarray([2]))

        assert xp.isnan(res.statistic)
        assert xp.isnan(res.pvalue)


@make_xp_test_case(stats.spearmanrho)
class TestSpearmanRho:
    @pytest.mark.parametrize('alternative, statistic, pvalue', [
        ('two-sided', -0.2727272727272727, 0.4458383415428),
        ('greater', -0.2727272727272727, 0.7770808292286),
        ('less', -0.2727272727272727, 0.2229191707714),
    ])
    @pytest.mark.parametrize('dtype', ['float32', 'float64', None])
    def test_against_R_cor_test(self, alternative, statistic, pvalue, dtype, xp):
        # Test against R cor.test, e.g.
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558,
        #       0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495)
        # y = c(0.8402081904493103, 0.5946972833914318, 0.23481606164114155,
        #       0.49754786197715384, 0.9146460831206026, 0.5848057749217579,
        #       0.7620801065573549, 0.31410063302647495, 0.7935620302236199,
        #       0.5423085761365468)
        # cor.test(x, y, method='spearman', alternative='t', exact=FALSE)
        dtype = xp_default_dtype(xp) if dtype is None else getattr(xp, dtype)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=10)
        y = rng.random(size=10)

        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.spearmanrho(x, y, alternative=alternative)

        xp_assert_close(res.statistic, xp.asarray(statistic, dtype=dtype))
        xp_assert_close(res.pvalue, xp.asarray(pvalue, dtype=dtype))


    @pytest.mark.parametrize('alternative, statistic, pvalue', [
        ('two-sided', -0.4857142857142857, 0.3555555555556),
        ('greater', -0.4857142857142857, 0.8513888888889),
        ('less', -0.4857142857142857, 0.1777777777778),
    ])
    def test_against_R_cor_test_exact(self, alternative, statistic, pvalue):
        xp = np  # will test for multiple backends when gh-23772 merges
        # Test against R cor.test exact=TRUE, e.g.
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558)
        # y = c(0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495, 0.8402081904493103, 0.5946972833914318)
        # cor.test(x, y, method='spearman', alternative='t', exact=TRUE)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=6)
        y = rng.random(size=6)
        x, y = xp.asarray(x), xp.asarray(y)

        method = stats.PermutationMethod()
        res = stats.spearmanrho(x, y, method=method, alternative=alternative)

        xp_assert_close(res.statistic, xp.asarray(statistic))
        xp_assert_close(res.pvalue, xp.asarray(pvalue))


    @pytest.mark.parametrize('alternative', ('two-sided', 'greater', 'less'))
    @pytest.mark.parametrize('n', [9, 99, 999])
    def test_against_scipy_spearmanr(self, alternative, n, xp):
        rng = np.random.default_rng(5982435982346983)
        x = rng.integers(n//2, size=n)
        y = rng.integers(n//2, size=n)
        dtype = xp_default_dtype(xp)

        ref = stats.spearmanr(x, y, alternative=alternative)
        ref_statistic = xp.asarray(ref.statistic, dtype=dtype)
        ref_pvalue = xp.asarray(ref.pvalue, dtype=dtype)

        x = xp.asarray(x, dtype=dtype)
        y = xp.asarray(y, dtype=dtype)
        res = stats.spearmanrho(x, y, alternative=alternative)

        xp_assert_close(res.statistic, ref_statistic)
        xp_assert_close(res.pvalue, ref_pvalue)


    @pytest.mark.parametrize('axis', [-1, 0, 1])
    def test_other_backends_nd(self, axis, xp):
        # NumPy n-d behavior is tested by `test_axis_nan_policy`;
        # check other backends against it.
        rng = np.random.default_rng(8243598346983259)
        shape = (8, 9, 10)
        x = rng.standard_normal(size=shape)
        y = rng.standard_normal(size=shape)
        res = stats.spearmanrho(xp.asarray(x), xp.asarray(y), axis=axis)
        ref = stats.spearmanrho(x, y, axis=axis)
        xp_assert_close(res.statistic, xp.asarray(ref.statistic), atol=1e-16)
        xp_assert_close(res.pvalue, xp.asarray(ref.pvalue), atol=1e-16)


    def test_input_validation(self, xp):
        rng = np.random.default_rng(25932435798274926)
        x, y = rng.random(size=(2, 10))
        x, y = xp.asarray(x), xp.asarray(y)

        msg = 'incompatible for broadcasting|Incompatible shapes|must be broadcastable'
        with pytest.raises((ValueError, TypeError), match=msg):
            stats.spearmanrho(x, y[:-1])

        if not is_jax(xp):
            # jax misses out on some input validation from _axis_nan_policy decorator
            message = '...axis 10 is out of bounds for array...|out of range'
            with pytest.raises((ValueError, IndexError), match=message):
                stats.spearmanrho(x, y, axis=10)

        message = "`alternative` must be 'less', 'greater', or 'two-sided'."
        with pytest.raises(ValueError, match=message):
            stats.spearmanrho(x, y, alternative='alternative')

        message = ("`method` must be...")
        with pytest.raises(ValueError, match=message):
            stats.spearmanrho(x, y, method='method')

    def test_special_cases(self, xp):
        def check_nan(res):
            assert xp.isnan(res.statistic)
            assert xp.isnan(res.pvalue)

        message = 'One or more sample arguments is too small...'
        with eager_warns(SmallSampleWarning, match=message, xp=xp):
            res = stats.spearmanrho(xp.asarray([1]), xp.asarray([2]))
            check_nan(res)

        x = xp.asarray([1, 1, 1, 1, 1])
        y = xp.asarray([1, 2, 3, 4, 5])
        message = 'An input array is constant; the correlation coefficient...'
        with eager_warns(stats.ConstantInputWarning, match=message, xp=xp):
            res = stats.spearmanrho(x, y)
            check_nan(res)
        with eager_warns(stats.ConstantInputWarning, match=message, xp=xp):
            res = stats.spearmanrho(y, x)
            check_nan(res)


class RobustSlopesTest:
    #  'iv, mask, warnings, consistency, gh19678, mstats'
    def test_input_validation(self, xp):
        pfun = getattr(stats, self.pfun)
        other_method = 'joint' if pfun == stats.theilslopes else 'hierarchical'
        msg = (f"method must be either '{other_method}' or 'separate'. "
               "'joint_separate' is invalid.")
        with pytest.raises(ValueError, match=msg):
            pfun(xp.asarray([0, 1, 1]), method='joint_separate')

    @skip_xp_invalid_arg
    @pytest.mark.parametrize("method", ['separate', 'other'])
    def test_mask(self, method):
        pfun = getattr(stats, self.pfun)
        if method == 'other':
            method = 'joint' if self.pfun =='theilslopes' else 'hierarchical'

        # Test for correct masking.
        mask = np.asarray([False, False, True, False])
        y = np.ma.array([0, 1, 100, 1], mask=mask)
        res = pfun(y, method=method)
        ref = pfun(y[~mask], method=method)

        assert_allclose(res.slope, ref.slope)
        assert_allclose(res.intercept, ref.intercept)

    def test_degenerate(self, xp):
        # Test with degenerate input; see gh-15943
        pfun = getattr(stats, self.pfun)
        res = pfun(xp.asarray([0, 1]), xp.asarray([0, 0]))
        assert xp.all(xp.isnan(xp.stack(res)))
        res = pfun(xp.asarray([0, 0, 0]), xp.asarray([0, 1, 0]))
        xp_assert_equal(res[0], xp.asarray(0.))
        xp_assert_equal(res[1], xp.asarray(0.))
        if pfun == stats.theilslopes:
            xp_assert_equal(res[2], xp.asarray(xp.nan))
            xp_assert_equal(res[3], xp.asarray(xp.nan))

    def test_namedtuple_consistency(self, xp):
        """
        Simple test to ensure tuple backwards-compatibility of the returned object.
        """
        pfun = getattr(stats, self.pfun)

        y = xp.asarray([1, 2, 4])
        x = xp.asarray([4, 6, 8])

        result = pfun(y, x)

        # note both returned values are distinct here
        xp_assert_equal(result[0], result.slope)
        xp_assert_equal(result[1], result.intercept)
        if pfun == stats.theilslopes:
            xp_assert_equal(result[2], result.low_slope)
            xp_assert_equal(result[3], result.high_slope)

    def test_gh19678_uint8(self, xp):
        # `theilslopes` returned unexpected results when `y` was an unsigned type.
        # Check that this is resolved.
        pfun = getattr(stats, self.pfun)
        rng = np.random.default_rng(2549824598234528)
        y = xp.asarray(rng.integers(0, 255, size=10, dtype=np.uint8))
        res = pfun(y, y)
        xp_assert_close(res.slope, xp.asarray(1.))

    @pytest.mark.parametrize("method", ['separate', 'other'])
    @pytest.mark.parametrize("case_number", [0, 1, 2, 3])
    def test_against_mstats(self, method, case_number, xp):
        pfun = getattr(stats, self.pfun)
        if method == 'other':
            method = 'joint' if self.pfun =='theilslopes' else 'hierarchical'

        rng = np.random.default_rng(349824598234528554)
        match case_number:
            case 0:  # no x
                x = None
                y = rng.random(100)
            case 1:  # no ties
                x = rng.random(100)
                y = rng.random(100)
            case 2:  # no x ties
                x = rng.random(100)
                y = rng.integers(50, size=100)
            case 3:  # ties in x and y
                x = rng.integers(25, size=100)
                y = rng.integers(50, size=100)

        ref = pfun(y, x, method=method)

        y, x = xp.asarray(y), x if x is None else xp.asarray(x)
        if (is_array_api_strict(xp) and xp.isdtype(y.dtype, 'integral')
                and xp.isdtype(x.dtype, 'real floating')):
            # array API strict doesn't allow mixed type promotion, so promote for it
            y = xp.astype(y, xp.float64)

        dtype = xp_result_type(y, x, force_floating=True, xp=xp)

        res = pfun(y, x, method=method, axis=-1)

        assert res.slope != 0  # ensure randomly generated test case isn't trivial
        xp_assert_close(res.slope, xp.asarray(ref.slope, dtype=dtype))
        xp_assert_close(res.intercept, xp.asarray(ref.intercept, dtype=dtype))
        if pfun == stats.theilslopes:
            xp_assert_close(res.low_slope, xp.asarray(ref.low_slope, dtype=dtype))
            xp_assert_close(res.high_slope, xp.asarray(ref.high_slope, dtype=dtype))


@make_xp_test_case(stats.theilslopes)
class TestTheilslopes(RobustSlopesTest):
    pfun = 'theilslopes'
    def test_theilslopes(self, xp):
        # Test for basic slope and intercept.
        y = xp.asarray([0, 1, 1])
        slope, intercept, lower, upper = stats.theilslopes(y)
        xp_assert_close(slope, xp.asarray(0.5))
        xp_assert_close(intercept, xp.asarray(0.5))

        slope, intercept, lower, upper = stats.theilslopes(y, method='joint')
        xp_assert_close(slope, xp.asarray(0.5))
        xp_assert_close(intercept, xp.asarray(0.0))

        # Test of confidence intervals from example in Sen (1968).
        x = xp.asarray([1, 2, 3, 4, 10, 12, 18])
        y = xp.asarray([9, 15, 19, 20, 45, 55, 78])
        slope, intercept, lower, upper = stats.theilslopes(y, x, 0.07)
        xp_assert_close(slope, xp.asarray(4.0))
        xp_assert_close(intercept, xp.asarray(4.0))
        xp_assert_close(upper, xp.asarray(4.38), rtol=5e-3)
        xp_assert_close(lower, xp.asarray(3.71), rtol=5e-3)

        slope, intercept, lower, upper = stats.theilslopes(y, x, 0.07,
                                                           method='joint')
        xp_assert_close(slope, xp.asarray(4.0))
        xp_assert_close(intercept, xp.asarray(6.0))
        xp_assert_close(upper, xp.asarray(4.38), rtol=5e-3)
        xp_assert_close(lower, xp.asarray(3.71), rtol=5e-3)


@make_xp_test_case(stats.siegelslopes)
class TestSiegelslopes(RobustSlopesTest):
    pfun = 'siegelslopes'
    def test_siegelslopes(self, xp):
        # method should be exact for straight line
        y = 2 * xp.arange(10.) + 0.5
        slope, intercept = stats.siegelslopes(y)
        xp_assert_close(slope, xp.asarray(2.0))
        xp_assert_close(intercept, xp.asarray(0.5))
        slope, intercept = stats.siegelslopes(y, method='separate')
        xp_assert_close(slope, xp.asarray(2.0))
        xp_assert_close(intercept, xp.asarray(0.5))

        x = 2 * xp.arange(10.)
        y = 5 * x - 3.0
        slope, intercept = stats.siegelslopes(y, x)
        xp_assert_close(slope, xp.asarray(5.0))
        xp_assert_close(intercept, xp.asarray(-3.0))
        slope, intercept = stats.siegelslopes(y, x, method='separate')
        xp_assert_close(slope, xp.asarray(5.0))
        xp_assert_close(intercept, xp.asarray(-3.0))

        # method is robust to outliers: breakdown point of 50%
        y = xpx.at(y)[:4].set(1000.)
        xp_assert_close(slope, xp.asarray(5.0))
        xp_assert_close(intercept, xp.asarray(-3.0))

        # if there are no outliers, results should be comparable to linregress
        x = np.arange(10.)
        y = -2.3 + 0.3 * x + stats.norm.rvs(size=10, random_state=231)
        slope_ols, intercept_ols, _, _, _ = stats.linregress(x, y)

        y, x = xp.asarray(y), xp.asarray(x)
        slope, intercept = stats.siegelslopes(y, x)
        xp_assert_close(slope, xp.asarray(slope_ols), rtol=0.1)
        xp_assert_close(intercept, xp.asarray(intercept_ols), rtol=0.1)

        slope, intercept = stats.siegelslopes(y, x, method='separate')
        xp_assert_close(slope, xp.asarray(slope_ols), rtol=0.1)
        xp_assert_close(intercept, xp.asarray(intercept_ols), rtol=0.1)

    def test_method(self, xp):
        # Test that distinguishes between the methods. Reference values generated with
        # SciPy 1.17 (results unchanged since at least 1.12, so assumed "correct") to
        # detect *changes*.
        y = xp.asarray([2, 1, 3, 0, 4])
        slope, intercept = stats.siegelslopes(y, method='separate')
        xp_assert_close(slope, xp.asarray(0.25))
        xp_assert_close(intercept, xp.asarray(1.75))
        slope, intercept = stats.siegelslopes(y, method='hierarchical')
        xp_assert_close(slope, xp.asarray(0.25))  # always the same between the methods
        xp_assert_close(intercept, xp.asarray(2.0))
