import copy

import matplotlib.pyplot as plt
from matplotlib.scale import (
    AsinhScale, AsinhTransform,
    FuncScale, LogitScale,
    LogTransform, InvertedLogTransform,
    SymmetricalLogTransform)
import matplotlib.scale as mscale
from matplotlib.ticker import (
    AsinhLocator, AutoLocator, LogFormatterSciNotation,
    NullFormatter, NullLocator, ScalarFormatter
)
from matplotlib.testing.decorators import check_figures_equal, image_comparison
from matplotlib.transforms import IdentityTransform

import numpy as np
from numpy.testing import assert_allclose
import io
import pytest


def test_optional_axis_signature():
    # There are three types of original signatures possible, and this only tests one
    # example class of each:
    # 1. `axis` without default: LinearScale, FuncScale, FuncScaleLog
    # 2. `axis` with default and more positional parameters: LogitScale
    # 3. `axis` with default and only keyword-only parameters: LogScale, AsinhScale,
    #    SymmetricalLogScale
    # Testing with None is sufficient as detection is purely based on the
    # signature structure; no type information is involved.
    axis = None

    # Old signature with axis positionally.
    FuncScale(axis, (lambda x: x, lambda x: x))
    FuncScale(axis, functions=(lambda x: x, lambda x: x))
    LogitScale(axis)
    LogitScale(axis, 'clip')
    LogitScale(axis, nonpositive='clip')
    LogitScale(axis, use_overline=True)
    AsinhScale(axis)
    AsinhScale(axis, linear_width=2)
    AsinhScale(axis, base=3)
    AsinhScale(axis, subs=[2, 6])
    # Old signature with axis as keyword.
    FuncScale(axis=axis, functions=(lambda x: x, lambda x: x))
    LogitScale(axis=axis)
    LogitScale(axis=axis, nonpositive='clip')
    LogitScale(axis=axis, use_overline=True)
    AsinhScale(axis=axis)
    AsinhScale(axis=axis, linear_width=2)
    AsinhScale(axis=axis, base=3)
    AsinhScale(axis=axis, subs=[2, 6])
    # New signature without axis.
    FuncScale((lambda x: x, lambda x: x))
    FuncScale(functions=(lambda x: x, lambda x: x))
    LogitScale()
    LogitScale(nonpositive='clip')
    LogitScale(use_overline=True)
    AsinhScale()
    AsinhScale(linear_width=2)
    AsinhScale(base=3)
    AsinhScale(subs=[2, 6])


@check_figures_equal()
def test_log_scales(fig_test, fig_ref):
    ax_test = fig_test.add_subplot(122, yscale='log', xscale='symlog')
    ax_test.axvline(24.1)
    ax_test.axhline(24.1)
    xlim = ax_test.get_xlim()
    ylim = ax_test.get_ylim()
    ax_ref = fig_ref.add_subplot(122, yscale='log', xscale='symlog')
    ax_ref.set(xlim=xlim, ylim=ylim)
    ax_ref.plot([24.1, 24.1], ylim, 'b')
    ax_ref.plot(xlim, [24.1, 24.1], 'b')


def test_symlog_mask_nan():
    # Use a transform round-trip to verify that the forward and inverse
    # transforms work, and that they respect nans and/or masking.
    slt = SymmetricalLogTransform(10, 2, 1)
    slti = slt.inverted()

    x = np.arange(-1.5, 5, 0.5)
    out = slti.transform_non_affine(slt.transform_non_affine(x))
    assert_allclose(out, x)
    assert type(out) is type(x)

    x[4] = np.nan
    out = slti.transform_non_affine(slt.transform_non_affine(x))
    assert_allclose(out, x)
    assert type(out) is type(x)

    x = np.ma.array(x)
    out = slti.transform_non_affine(slt.transform_non_affine(x))
    assert_allclose(out, x)
    assert type(out) is type(x)

    x[3] = np.ma.masked
    out = slti.transform_non_affine(slt.transform_non_affine(x))
    assert_allclose(out, x)
    assert type(out) is type(x)


@image_comparison(['logit_scales.png'], remove_text=True, style='_classic_test')
def test_logit_scales():
    fig, ax = plt.subplots()

    # Typical extinction curve for logit
    x = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5,
                  0.6, 0.7, 0.8, 0.9, 0.97, 0.99, 0.997, 0.999])
    y = 1.0 / x

    ax.plot(x, y)
    ax.set_xscale('logit')
    ax.grid(True)
    bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    assert np.isfinite(bbox.x0)
    assert np.isfinite(bbox.y0)


def test_log_scatter():
    """Issue #1799"""
    fig, ax = plt.subplots(1)

    x = np.arange(10)
    y = np.arange(10) - 1

    ax.scatter(x, y)

    buf = io.BytesIO()
    fig.savefig(buf, format='pdf')

    buf = io.BytesIO()
    fig.savefig(buf, format='eps')

    buf = io.BytesIO()
    fig.savefig(buf, format='svg')


def test_logscale_subs():
    fig, ax = plt.subplots()
    ax.set_yscale('log', subs=np.array([2, 3, 4]))
    # force draw
    fig.canvas.draw()


@image_comparison(['logscale_mask.png'], remove_text=True, style='_classic_test')
def test_logscale_mask():
    # Check that zero values are masked correctly on log scales.
    # See github issue 8045
    xs = np.linspace(0, 50, 1001)

    fig, ax = plt.subplots()
    ax.plot(np.exp(-xs**2))
    fig.canvas.draw()
    ax.set(yscale="log",
           yticks=10.**np.arange(-300, 0, 24))  # Backcompat tick selection.


def test_extra_kwargs_raise():
    fig, ax = plt.subplots()

    for scale in ['linear', 'log', 'symlog']:
        with pytest.raises(TypeError):
            ax.set_yscale(scale, foo='mask')


def test_logscale_invert_transform():
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    # get transformation from data to axes
    tform = (ax.transAxes + ax.transData.inverted()).inverted()

    # direct test of log transform inversion
    inverted_transform = LogTransform(base=2).inverted()
    assert isinstance(inverted_transform, InvertedLogTransform)
    assert inverted_transform.base == 2


def test_logscale_transform_repr():
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    repr(ax.transData)
    repr(LogTransform(10, nonpositive='clip'))


@image_comparison(['logscale_nonpos_values.png'],
                  remove_text=True, tol=0.02, style='mpl20')
def test_logscale_nonpos_values():
    np.random.seed(19680801)
    xs = np.random.normal(size=int(1e3))
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.hist(xs, range=(-5, 5), bins=10)
    ax1.set_yscale('log')
    ax2.hist(xs, range=(-5, 5), bins=10)
    ax2.set_yscale('log', nonpositive='mask')

    xdata = np.arange(0, 10, 0.01)
    ydata = np.exp(-xdata)
    edata = 0.2*(10-xdata)*np.cos(5*xdata)*np.exp(-xdata)

    ax3.fill_between(xdata, ydata - edata, ydata + edata)
    ax3.set_yscale('log')

    x = np.logspace(-1, 1)
    y = x ** 3
    yerr = x**2
    ax4.errorbar(x, y, yerr=yerr)

    ax4.set_yscale('log')
    ax4.set_xscale('log')
    ax4.set_yticks([1e-2, 1, 1e+2])  # Backcompat tick selection.


def test_invalid_log_lims():
    # Check that invalid log scale limits are ignored
    fig, ax = plt.subplots()
    ax.scatter(range(0, 4), range(0, 4))

    ax.set_xscale('log')
    original_xlim = ax.get_xlim()
    with pytest.warns(UserWarning):
        ax.set_xlim(left=0)
    assert ax.get_xlim() == original_xlim
    with pytest.warns(UserWarning):
        ax.set_xlim(right=-1)
    assert ax.get_xlim() == original_xlim

    ax.set_yscale('log')
    original_ylim = ax.get_ylim()
    with pytest.warns(UserWarning):
        ax.set_ylim(bottom=0)
    assert ax.get_ylim() == original_ylim
    with pytest.warns(UserWarning):
        ax.set_ylim(top=-1)
    assert ax.get_ylim() == original_ylim


@image_comparison(['function_scales.png'], remove_text=True, style='mpl20')
def test_function_scale():
    def inverse(x):
        return x**2

    def forward(x):
        return x**(1/2)

    fig, ax = plt.subplots()

    x = np.arange(1, 1000)

    ax.plot(x, x)
    ax.set_xscale('function', functions=(forward, inverse))
    ax.set_xlim(1, 1000)


def test_pass_scale():
    # test passing a scale object works...
    fig, ax = plt.subplots()
    scale = mscale.LogScale(axis=None)
    ax.set_xscale(scale)
    scale = mscale.LogScale(axis=None)
    ax.set_yscale(scale)
    assert ax.xaxis.get_scale() == 'log'
    assert ax.yaxis.get_scale() == 'log'


def test_scale_deepcopy():
    sc = mscale.LogScale(axis='x', base=10)
    sc2 = copy.deepcopy(sc)
    assert str(sc.get_transform()) == str(sc2.get_transform())
    assert sc._transform is not sc2._transform


class TestAsinhScale:
    def test_transforms(self):
        a0 = 17.0
        a = np.linspace(-50, 50, 100)

        forward = AsinhTransform(a0)
        inverse = forward.inverted()
        invinv = inverse.inverted()

        a_forward = forward.transform_non_affine(a)
        a_inverted = inverse.transform_non_affine(a_forward)
        assert_allclose(a_inverted, a)

        a_invinv = invinv.transform_non_affine(a)
        assert_allclose(a_invinv, a0 * np.arcsinh(a / a0))

    def test_init(self):
        fig, ax = plt.subplots()

        s = AsinhScale(axis=None, linear_width=23.0)
        assert s.linear_width == 23
        assert s._base == 10
        assert s._subs == (2, 5)

        tx = s.get_transform()
        assert isinstance(tx, AsinhTransform)
        assert tx.linear_width == s.linear_width

    def test_base_init(self):
        fig, ax = plt.subplots()

        s3 = AsinhScale(axis=None, base=3)
        assert s3._base == 3
        assert s3._subs == (2,)

        s7 = AsinhScale(axis=None, base=7, subs=(2, 4))
        assert s7._base == 7
        assert s7._subs == (2, 4)

    def test_fmtloc(self):
        class DummyAxis:
            def __init__(self):
                self.fields = {}
            def set(self, **kwargs):
                self.fields.update(**kwargs)
            def set_major_formatter(self, f):
                self.fields['major_formatter'] = f

        ax0 = DummyAxis()
        s0 = AsinhScale(axis=ax0, base=0)
        s0.set_default_locators_and_formatters(ax0)
        assert isinstance(ax0.fields['major_locator'], AsinhLocator)
        assert isinstance(ax0.fields['major_formatter'], str)

        ax5 = DummyAxis()
        s7 = AsinhScale(axis=ax5, base=5)
        s7.set_default_locators_and_formatters(ax5)
        assert isinstance(ax5.fields['major_locator'], AsinhLocator)
        assert isinstance(ax5.fields['major_formatter'],
                          LogFormatterSciNotation)

    def test_bad_scale(self):
        fig, ax = plt.subplots()

        with pytest.raises(ValueError):
            AsinhScale(axis=None, linear_width=0)
        with pytest.raises(ValueError):
            AsinhScale(axis=None, linear_width=-1)
        s0 = AsinhScale(axis=None, )
        s1 = AsinhScale(axis=None, linear_width=3.0)


def test_custom_scale_without_axis():
    """
    Test that one can register and use custom scales that don't take an *axis* param.
    """
    class CustomTransform(IdentityTransform):
        pass

    class CustomScale(mscale.ScaleBase):
        name = "custom"

        # Important: __init__ has no *axis* parameter
        def __init__(self):
            self._transform = CustomTransform()

        def get_transform(self):
            return self._transform

        def set_default_locators_and_formatters(self, axis):
            axis.set_major_locator(AutoLocator())
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_locator(NullLocator())
            axis.set_minor_formatter(NullFormatter())

    try:
        mscale.register_scale(CustomScale)
        fig, ax = plt.subplots()
        ax.set_xscale('custom')
        assert isinstance(ax.xaxis.get_transform(), CustomTransform)
    finally:
        # cleanup - there's no public unregister_scale()
        del mscale._scale_mapping["custom"]
        del mscale._scale_has_axis_parameter["custom"]


def test_custom_scale_with_axis():
    """
    Test that one can still register and use custom scales with an *axis*
    parameter, but that registering issues a pending-deprecation warning.
    """
    class CustomTransform(IdentityTransform):
        pass

    class CustomScale(mscale.ScaleBase):
        name = "custom"

        # Important: __init__ still has the *axis* parameter
        def __init__(self, axis):
            self._transform = CustomTransform()

        def get_transform(self):
            return self._transform

        def set_default_locators_and_formatters(self, axis):
            axis.set_major_locator(AutoLocator())
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_locator(NullLocator())
            axis.set_minor_formatter(NullFormatter())

    try:
        with pytest.warns(
                PendingDeprecationWarning,
                match=r"'axis' parameter .* is pending-deprecated"):
            mscale.register_scale(CustomScale)
        fig, ax = plt.subplots()
        ax.set_xscale('custom')
        assert isinstance(ax.xaxis.get_transform(), CustomTransform)
    finally:
        # cleanup - there's no public unregister_scale()
        del mscale._scale_mapping["custom"]
        del mscale._scale_has_axis_parameter["custom"]


def test_val_in_range():

    test_cases = [
        # LinearScale: Always True (even for Inf/NaN)
        ('linear', 10.0, True),
        ('linear', -10.0, True),
        ('linear', 0.0, True),
        ('linear', np.inf, False),
        ('linear', np.nan, False),

        # LogScale: Only positive values (> 0)
        ('log', 1.0, True),
        ('log', 1e-300, True),
        ('log', 0.0, False),
        ('log', -1.0, False),
        ('log', np.inf, False),
        ('log', np.nan, False),

        # LogitScale: Strictly between 0 and 1
        ('logit', 0.5, True),
        ('logit', 0.0, False),
        ('logit', 1.0, False),
        ('logit', -0.1, False),
        ('logit', 1.1, False),
        ('logit', np.inf, False),
        ('logit', np.nan, False),

        # SymmetricalLogScale: Valid for all real numbers
        # Uses ScaleBase fallback. NaN returns False since NaN != NaN
        ('symlog', 10.0, True),
        ('symlog', -10.0, True),
        ('symlog', 0.0, True),
        ('symlog', np.inf, False),
        ('symlog', np.nan, False),
    ]

    for name, val, expected in test_cases:
        scale_cls = mscale._scale_mapping[name]
        s = scale_cls(axis=None)

        result = s.val_in_range(val)
        assert result is expected, (
            f"Failed {name}.val_in_range({val})."
            f"Expected {expected}, got {result}"
        )


def test_val_in_range_base_fallback():
    # Directly test the ScaleBase fallback for custom scales.
    # ScaleBase.limit_range_for_scale returns values unchanged by default
    s = mscale.ScaleBase(axis=None)

    # Normal values should be True
    assert s.val_in_range(1.0) is True
    assert s.val_in_range(-5.5) is True

    # NaN and Inf returns False since they cannot be drawn in a plot
    assert s.val_in_range(np.nan) is False
    assert s.val_in_range(np.inf) is False
    assert s.val_in_range(-np.inf) is False


def test_val_in_range_array():
    # Vectorized: scalar in -> scalar bool, array in -> bool array.
    arr = np.array([0.5, -1.0, 0.0, np.nan, np.inf, 0.25])
    cases = {
        'linear': [True, True, True, False, False, True],
        'log':    [True, False, False, False, False, True],
        'symlog': [True, True, True, False, False, True],
        'asinh':  [True, True, True, False, False, True],
        'logit':  [True, False, False, False, False, True],
    }
    for name, expected in cases.items():
        s = mscale._scale_mapping[name](axis=None)
        np.testing.assert_array_equal(s.val_in_range(arr), expected)

    # 2D shape is preserved.
    out = mscale._scale_mapping['log'](axis=None).val_in_range(
        np.array([[1.0, -1.0], [0.5, np.nan]]))
    np.testing.assert_array_equal(out, [[True, False], [True, False]])
