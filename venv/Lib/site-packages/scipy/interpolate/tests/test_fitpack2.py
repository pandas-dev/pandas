# Created by Pearu Peterson, June 2003
import itertools
import sys
import warnings

import numpy as np
import pytest
from pytest import raises as assert_raises
from scipy._lib._array_api import (
    xp_assert_equal, xp_assert_close, assert_almost_equal, assert_array_almost_equal
)

from numpy import array, diff, linspace, meshgrid, ones, pi, shape
from scipy.interpolate._fitpack_py import bisplrep, bisplev, splrep, spalde
from scipy.interpolate._fitpack2 import (UnivariateSpline,
        LSQUnivariateSpline, InterpolatedUnivariateSpline,
        LSQBivariateSpline, SmoothBivariateSpline, RectBivariateSpline,
        LSQSphereBivariateSpline, SmoothSphereBivariateSpline,
        RectSphereBivariateSpline)

from scipy._lib._testutils import _run_concurrent_barrier

from scipy.interpolate import make_splrep, NdBSpline
from scipy.interpolate._regrid import (_regrid,
        _ndbspline_call_like_bivariate)

def convert_to_ndbspline(lut):
    tx, ty = lut.get_knots()
    kx, ky = lut.degrees
    nx, ny = len(tx), len(ty)
    c = lut.get_coeffs().reshape((nx - kx - 1, ny - ky - 1))
    return NdBSpline((tx, ty), c, (kx, ky))

class TestUnivariateSpline:
    def test_linear_constant(self):
        x = [1,2,3]
        y = [3,3,3]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(), [1, 3])
        assert_array_almost_equal(lut.get_coeffs(), [3, 3])
        assert abs(lut.get_residual()) < 1e-10
        assert_array_almost_equal(lut([1, 1.5, 2]), [3, 3, 3])

    @pytest.mark.parametrize("bc_type", [None, 'periodic'])
    def test_linear_constant_periodic(self, bc_type):
        x = [1,2,3]
        y = [3,3,3]
        lut = UnivariateSpline(x,y,k=1)

        spl = make_splrep(x, y, k=1, s=len(x), bc_type=bc_type)
        xp_assert_close(spl.t[1:-1], lut.get_knots(), atol=1e-15)
        xp_assert_close(spl.c, lut.get_coeffs(), atol=1e-15)

    def test_preserve_shape(self):
        x = [1, 2, 3]
        y = [0, 2, 4]
        lut = UnivariateSpline(x, y, k=1)
        arg = 2
        assert shape(arg) == shape(lut(arg))
        assert shape(arg) == shape(lut(arg, nu=1))
        arg = [1.5, 2, 2.5]
        assert shape(arg) == shape(lut(arg))
        assert shape(arg) == shape(lut(arg, nu=1))

    def test_linear_1d(self):
        x = [1,2,3]
        y = [0,2,4]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[0,4])
        assert abs(lut.get_residual()) < 1e-15
        assert_array_almost_equal(lut([1,1.5,2]),[0,1,2])

    def test_subclassing(self):
        # See #731

        class ZeroSpline(UnivariateSpline):
            def __call__(self, x):
                return 0*array(x)

        sp = ZeroSpline([1,2,3,4,5], [3,2,3,2,3], k=2)
        xp_assert_equal(sp([1.5, 2.5]), [0., 0.])

    def test_empty_input(self):
        # Test whether empty input returns an empty output. Ticket 1014
        x = [1,3,5,7,9]
        y = [0,4,9,12,21]
        spl = UnivariateSpline(x, y, k=3)
        xp_assert_equal(spl([]), array([]))

    def test_roots(self):
        x = [1, 3, 5, 7, 9]
        y = [0, 4, 9, 12, 21]
        spl = UnivariateSpline(x, y, k=3)
        assert_almost_equal(spl.roots()[0], 1.050290639101332)

    def test_roots_length(self): # for gh18335
        x = np.linspace(0, 50 * np.pi, 1000)
        y = np.cos(x)
        spl = UnivariateSpline(x, y, s=0)
        assert len(spl.roots()) == 50

    def test_derivatives(self):
        x = [1, 3, 5, 7, 9]
        y = [0, 4, 9, 12, 21]
        spl = UnivariateSpline(x, y, k=3)
        assert_almost_equal(spl.derivatives(3.5),
                            [5.5152902, 1.7146577, -0.1830357, 0.3125])

    def test_derivatives_2(self):
        x = np.arange(8)
        y = x**3 + 2.*x**2

        tck = splrep(x, y, s=0)
        ders = spalde(3, tck)
        xp_assert_close(ders, [45.,   # 3**3 + 2*(3)**2
                               39.,   # 3*(3)**2 + 4*(3)
                               22.,   # 6*(3) + 4
                               6.],   # 6*3**0
                        atol=1e-15)
        spl = UnivariateSpline(x, y, s=0, k=3)
        xp_assert_close(spl.derivatives(3),
                        ders,
                        atol=1e-15)

    def test_resize_regression(self):
        """Regression test for #1375."""
        x = [-1., -0.65016502, -0.58856235, -0.26903553, -0.17370892,
             -0.10011001, 0., 0.10011001, 0.17370892, 0.26903553, 0.58856235,
             0.65016502, 1.]
        y = [1.,0.62928599, 0.5797223, 0.39965815, 0.36322694, 0.3508061,
             0.35214793, 0.3508061, 0.36322694, 0.39965815, 0.5797223,
             0.62928599, 1.]
        w = [1.00000000e+12, 6.88875973e+02, 4.89314737e+02, 4.26864807e+02,
             6.07746770e+02, 4.51341444e+02, 3.17480210e+02, 4.51341444e+02,
             6.07746770e+02, 4.26864807e+02, 4.89314737e+02, 6.88875973e+02,
             1.00000000e+12]
        spl = UnivariateSpline(x=x, y=y, w=w, s=None)
        desired = array([0.35100374, 0.51715855, 0.87789547, 0.98719344])
        xp_assert_close(spl([0.1, 0.5, 0.9, 0.99]), desired, atol=5e-4)

    def test_out_of_range_regression(self):
        # Test different extrapolation modes. See ticket 3557
        x = np.arange(5, dtype=float)
        y = x**3

        xp = linspace(-8, 13, 100)
        xp_zeros = xp.copy()
        xp_zeros[np.logical_or(xp_zeros < 0., xp_zeros > 4.)] = 0
        xp_clip = xp.copy()
        xp_clip[xp_clip < x[0]] = x[0]
        xp_clip[xp_clip > x[-1]] = x[-1]

        for cls in [UnivariateSpline, InterpolatedUnivariateSpline]:
            spl = cls(x=x, y=y)
            for ext in [0, 'extrapolate']:
                xp_assert_close(spl(xp, ext=ext), xp**3, atol=1e-16)
                xp_assert_close(cls(x, y, ext=ext)(xp), xp**3, atol=1e-16)
            for ext in [1, 'zeros']:
                xp_assert_close(spl(xp, ext=ext), xp_zeros**3, atol=1e-16)
                xp_assert_close(cls(x, y, ext=ext)(xp), xp_zeros**3, atol=1e-16)
            for ext in [2, 'raise']:
                assert_raises(ValueError, spl, xp, **dict(ext=ext))
            for ext in [3, 'const']:
                xp_assert_close(spl(xp, ext=ext), xp_clip**3, atol=2e-16)
                xp_assert_close(cls(x, y, ext=ext)(xp), xp_clip**3, atol=2e-16)

        # also test LSQUnivariateSpline [which needs explicit knots]
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        spl = LSQUnivariateSpline(x, y, t)
        xp_assert_close(spl(xp, ext=0), xp**3, atol=1e-16)
        xp_assert_close(spl(xp, ext=1), xp_zeros**3, atol=1e-16)
        assert_raises(ValueError, spl, xp, **dict(ext=2))
        xp_assert_close(spl(xp, ext=3), xp_clip**3, atol=1e-16)

        # also make sure that unknown values for `ext` are caught early
        for ext in [-1, 'unknown']:
            spl = UnivariateSpline(x, y)
            assert_raises(ValueError, spl, xp, **dict(ext=ext))
            assert_raises(ValueError, UnivariateSpline,
                    **dict(x=x, y=y, ext=ext))

    def test_lsq_fpchec(self):
        xs = np.arange(100) * 1.
        ys = np.arange(100) * 1.
        knots = np.linspace(0, 99, 10)
        bbox = (-1, 101)
        assert_raises(ValueError, LSQUnivariateSpline, xs, ys, knots,
                      bbox=bbox)

    def test_derivative_and_antiderivative(self):
        # Thin wrappers to splder/splantider, so light smoke test only.
        x = np.linspace(0, 1, 70)**3
        y = np.cos(x)

        spl = UnivariateSpline(x, y, s=0)
        spl2 = spl.antiderivative(2).derivative(2)
        xp_assert_close(spl(0.3), spl2(0.3))

        spl2 = spl.antiderivative(1)
        xp_assert_close(spl2(0.6) - spl2(0.2),
                        spl.integral(0.2, 0.6))

    def test_derivative_extrapolation(self):
        # Regression test for gh-10195: for a const-extrapolation spline
        # its derivative evaluates to zero for extrapolation
        x_values = [1, 2, 4, 6, 8.5]
        y_values = [0.5, 0.8, 1.3, 2.5, 5]
        f = UnivariateSpline(x_values, y_values, ext='const', k=3)

        x = [-1, 0, -0.5, 9, 9.5, 10]
        xp_assert_close(f.derivative()(x), np.zeros_like(x), atol=1e-15)

    def test_integral_out_of_bounds(self):
        # Regression test for gh-7906: .integral(a, b) is wrong if both
        # a and b are out-of-bounds
        x = np.linspace(0., 1., 7)
        for ext in range(4):
            f = UnivariateSpline(x, x, s=0, ext=ext)
            for (a, b) in [(1, 1), (1, 5), (2, 5),
                           (0, 0), (-2, 0), (-2, -1)]:
                assert abs(f.integral(a, b)) < 1e-15

    def test_nan(self):
        # bail out early if the input data contains nans
        x = np.arange(10, dtype=float)
        y = x**3
        w = np.ones_like(x)
        # also test LSQUnivariateSpline [which needs explicit knots]
        spl = UnivariateSpline(x, y, check_finite=True)
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        y_end = y[-1]
        for z in [np.nan, np.inf, -np.inf]:
            y[-1] = z
            assert_raises(ValueError, UnivariateSpline,
                    **dict(x=x, y=y, check_finite=True))
            assert_raises(ValueError, InterpolatedUnivariateSpline,
                    **dict(x=x, y=y, check_finite=True))
            assert_raises(ValueError, LSQUnivariateSpline,
                    **dict(x=x, y=y, t=t, check_finite=True))
            y[-1] = y_end  # check valid y but invalid w
            w[-1] = z
            assert_raises(ValueError, UnivariateSpline,
                    **dict(x=x, y=y, w=w, check_finite=True))
            assert_raises(ValueError, InterpolatedUnivariateSpline,
                    **dict(x=x, y=y, w=w, check_finite=True))
            assert_raises(ValueError, LSQUnivariateSpline,
                    **dict(x=x, y=y, t=t, w=w, check_finite=True))

    def test_strictly_increasing_x(self):
        # Test the x is required to be strictly increasing for
        # UnivariateSpline if s=0 and for InterpolatedUnivariateSpline,
        # but merely increasing for UnivariateSpline if s>0
        # and for LSQUnivariateSpline; see gh-8535
        xx = np.arange(10, dtype=float)
        yy = xx**3
        x = np.arange(10, dtype=float)
        x[1] = x[0]
        y = x**3
        w = np.ones_like(x)
        # also test LSQUnivariateSpline [which needs explicit knots]
        spl = UnivariateSpline(xx, yy, check_finite=True)
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        UnivariateSpline(x=x, y=y, w=w, s=1, check_finite=True)
        LSQUnivariateSpline(x=x, y=y, t=t, w=w, check_finite=True)
        assert_raises(ValueError, UnivariateSpline,
                **dict(x=x, y=y, s=0, check_finite=True))
        assert_raises(ValueError, InterpolatedUnivariateSpline,
                **dict(x=x, y=y, check_finite=True))

    def test_increasing_x(self):
        # Test that x is required to be increasing, see gh-8535
        xx = np.arange(10, dtype=float)
        yy = xx**3
        x = np.arange(10, dtype=float)
        x[1] = x[0] - 1.0
        y = x**3
        w = np.ones_like(x)
        # also test LSQUnivariateSpline [which needs explicit knots]
        spl = UnivariateSpline(xx, yy, check_finite=True)
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        assert_raises(ValueError, UnivariateSpline,
                **dict(x=x, y=y, check_finite=True))
        assert_raises(ValueError, InterpolatedUnivariateSpline,
                **dict(x=x, y=y, check_finite=True))
        assert_raises(ValueError, LSQUnivariateSpline,
                **dict(x=x, y=y, t=t, w=w, check_finite=True))

    def test_invalid_input_for_univariate_spline(self):

        with assert_raises(ValueError) as info:
            x_values = [1, 2, 4, 6, 8.5]
            y_values = [0.5, 0.8, 1.3, 2.5]
            UnivariateSpline(x_values, y_values)
        assert "x and y should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            x_values = [1, 2, 4, 6, 8.5]
            y_values = [0.5, 0.8, 1.3, 2.5, 2.8]
            w_values = [-1.0, 1.0, 1.0, 1.0]
            UnivariateSpline(x_values, y_values, w=w_values)
        assert "x, y, and w should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            bbox = (-1)
            UnivariateSpline(x_values, y_values, bbox=bbox)
        assert "bbox shape should be (2,)" in str(info.value)

        with assert_raises(ValueError) as info:
            UnivariateSpline(x_values, y_values, k=6)
        assert "k should be 1 <= k <= 5" in str(info.value)

        with assert_raises(ValueError) as info:
            UnivariateSpline(x_values, y_values, s=-1.0)
        assert "s should be s >= 0.0" in str(info.value)

    def test_invalid_input_for_interpolated_univariate_spline(self):

        with assert_raises(ValueError) as info:
            x_values = [1, 2, 4, 6, 8.5]
            y_values = [0.5, 0.8, 1.3, 2.5]
            InterpolatedUnivariateSpline(x_values, y_values)
        assert "x and y should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            x_values = [1, 2, 4, 6, 8.5]
            y_values = [0.5, 0.8, 1.3, 2.5, 2.8]
            w_values = [-1.0, 1.0, 1.0, 1.0]
            InterpolatedUnivariateSpline(x_values, y_values, w=w_values)
        assert "x, y, and w should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            bbox = (-1)
            InterpolatedUnivariateSpline(x_values, y_values, bbox=bbox)
        assert "bbox shape should be (2,)" in str(info.value)

        with assert_raises(ValueError) as info:
            InterpolatedUnivariateSpline(x_values, y_values, k=6)
        assert "k should be 1 <= k <= 5" in str(info.value)

    def test_invalid_input_for_lsq_univariate_spline(self):

        x_values = [1, 2, 4, 6, 8.5]
        y_values = [0.5, 0.8, 1.3, 2.5, 2.8]
        spl = UnivariateSpline(x_values, y_values, check_finite=True)
        t_values = spl.get_knots()[3:4]  # interior knots w/ default k=3

        with assert_raises(ValueError) as info:
            x_values = [1, 2, 4, 6, 8.5]
            y_values = [0.5, 0.8, 1.3, 2.5]
            LSQUnivariateSpline(x_values, y_values, t_values)
        assert "x and y should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            x_values = [1, 2, 4, 6, 8.5]
            y_values = [0.5, 0.8, 1.3, 2.5, 2.8]
            w_values = [1.0, 1.0, 1.0, 1.0]
            LSQUnivariateSpline(x_values, y_values, t_values, w=w_values)
        assert "x, y, and w should have a same length" in str(info.value)

        message = "Interior knots t must satisfy Schoenberg-Whitney conditions"
        with assert_raises(ValueError, match=message) as info:
            bbox = (100, -100)
            LSQUnivariateSpline(x_values, y_values, t_values, bbox=bbox)

        with assert_raises(ValueError) as info:
            bbox = (-1)
            LSQUnivariateSpline(x_values, y_values, t_values, bbox=bbox)
        assert "bbox shape should be (2,)" in str(info.value)

        with assert_raises(ValueError) as info:
            LSQUnivariateSpline(x_values, y_values, t_values, k=6)
        assert "k should be 1 <= k <= 5" in str(info.value)

    def test_array_like_input(self):
        x_values = np.array([1, 2, 4, 6, 8.5])
        y_values = np.array([0.5, 0.8, 1.3, 2.5, 2.8])
        w_values = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        bbox = np.array([-100, 100])
        # np.array input
        spl1 = UnivariateSpline(x=x_values, y=y_values, w=w_values,
                                bbox=bbox)
        # list input
        spl2 = UnivariateSpline(x=x_values.tolist(), y=y_values.tolist(),
                                w=w_values.tolist(), bbox=bbox.tolist())

        xp_assert_close(spl1([0.1, 0.5, 0.9, 0.99]),
                        spl2([0.1, 0.5, 0.9, 0.99]))

    def test_fpknot_oob_crash(self):
        # https://github.com/scipy/scipy/issues/3691
        x = range(109)
        y = [0., 0., 0., 0., 0., 10.9, 0., 11., 0.,
             0., 0., 10.9, 0., 0., 0., 0., 0., 0.,
             10.9, 0., 0., 0., 11., 0., 0., 0., 10.9,
             0., 0., 0., 10.5, 0., 0., 0., 10.7, 0.,
             0., 0., 11., 0., 0., 0., 0., 0., 0.,
             10.9, 0., 0., 10.7, 0., 0., 0., 10.6, 0.,
             0., 0., 10.5, 0., 0., 10.7, 0., 0., 10.5,
             0., 0., 11.5, 0., 0., 0., 10.7, 0., 0.,
             10.7, 0., 0., 10.9, 0., 0., 10.8, 0., 0.,
             0., 10.7, 0., 0., 10.6, 0., 0., 0., 10.4,
             0., 0., 10.6, 0., 0., 10.5, 0., 0., 0.,
             10.7, 0., 0., 0., 10.4, 0., 0., 0., 10.8, 0.]
        msg = r"does not satisfy the condition abs\(fp-s\)/s < tol"
        with pytest.warns(UserWarning, match=msg):
            UnivariateSpline(x, y, k=1)

    def test_concurrency(self):
        # Check that no segfaults appear with concurrent access to
        # UnivariateSpline
        xx = np.arange(100, dtype=float)
        yy = xx**3
        x = np.arange(100, dtype=float)
        x[1] = x[0]
        spl = UnivariateSpline(xx, yy, check_finite=True)

        def worker_fn(_, interp, x):
            interp(x)

        _run_concurrent_barrier(10, worker_fn, spl, x)

    def test_curfit_2d_array_handling(self):
        """Test that 2D arrays are flattened correctly in Fortran order."""

        # Create 1D data
        x_1d = np.array([0., 1., 2., 3., 4., 5.])
        y_1d = np.array([0., 1., 4., 9., 16., 25.])  # y = x^2

        # Column vector (shape (n, 1)); should work and give same result as 1D
        x_col = x_1d.reshape(-1, 1)
        y_col = y_1d.reshape(-1, 1)

        spl_1d = InterpolatedUnivariateSpline(x_1d, y_1d)
        spl_col = InterpolatedUnivariateSpline(x_col, y_col)

        test_points = np.array([0.5, 1.5, 2.5, 3.5])
        np.testing.assert_allclose(spl_1d(test_points), spl_col(test_points))

        # Row vector (shape (1, n)) - should work and give same result as 1D
        x_row = x_1d.reshape(1, -1)
        y_row = y_1d.reshape(1, -1)

        spl_row = InterpolatedUnivariateSpline(x_row, y_row)
        np.testing.assert_allclose(spl_1d(test_points), spl_row(test_points))

        # Verify F-order flattening for truly 2D case. Not a likely input but
        # for general testing purposes.
        #
        # Create a 2x3 array where F-order and C-order give different results
        # F-order: [1,4,2,5,3,6], C-order: [1,2,3,4,5,6]
        x_2d = np.array([[1., 2., 3.],
                        [4., 5., 6.]])
        y_2d = np.array([[1., 4., 9.],
                        [16., 25., 36.]])

        # With F-order flattening: x = [1,4,2,5,3,6], y = [1,16,4,25,9,36]
        # This would NOT be sorted, so the spline should fail or give different results
        # than C-order which would give sorted x = [1,2,3,4,5,6]

        # The F-order flattened x is not monotonic, so FITPACK returns ier=10
        # and emits a UserWarning about erronous input.
        with pytest.warns(UserWarning, match="x\\[0\\]<x\\[1\\]<"):
            InterpolatedUnivariateSpline(x_2d, y_2d)


class TestLSQBivariateSpline:
    # NOTE: The systems in this test class are rank-deficient
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with pytest.warns(UserWarning, match="\nThe coefficients of the spline") as r:
            lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)
            assert len(r) == 1

        assert_almost_equal(lut(2, 2), np.asarray(3.))

    def test_bilinearity(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,7,8,3,4,7,1,3,4]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with pytest.warns(UserWarning, match="\nThe coefficients of the spline"):
            # This seems to fail (ier=1, see ticket 1642).
            lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)

        tx, ty = lut.get_knots()
        for xa, xb in zip(tx[:-1], tx[1:]):
            for ya, yb in zip(ty[:-1], ty[1:]):
                for t in [0.1, 0.5, 0.9]:
                    for s in [0.3, 0.4, 0.7]:
                        xp = xa*(1-t) + xb*t
                        yp = ya*(1-s) + yb*s
                        zp = (+ lut(xa, ya)*(1-t)*(1-s)
                              + lut(xb, ya)*t*(1-s)
                              + lut(xa, yb)*(1-t)*s
                              + lut(xb, yb)*t*s)
                        assert_almost_equal(lut(xp,yp), zp)

    def test_integral(self):
        x = [1,1,1,2,2,2,8,8,8]
        y = [1,2,3,1,2,3,1,2,3]
        z = array([0,7,8,3,4,7,1,3,4])

        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with pytest.warns(UserWarning, match="\nThe coefficients of the spline") as r:
            lut = LSQBivariateSpline(x, y, z, tx, ty, kx=1, ky=1)
            assert len(r) == 1
        tx, ty = lut.get_knots()
        tz = lut(tx, ty)
        trpz = .25*(diff(tx)[:,None]*diff(ty)[None,:]
                    * (tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()

        assert_almost_equal(np.asarray(lut.integral(tx[0], tx[-1], ty[0], ty[-1])),
                            np.asarray(trpz))

    def test_empty_input(self):
        # Test whether empty inputs returns an empty output. Ticket 1014
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with pytest.warns(UserWarning, match="\nThe coefficients of the spline") as r:
            lut = LSQBivariateSpline(x, y, z, tx, ty, kx=1, ky=1)
            assert len(r) == 1

        xp_assert_equal(lut([], []), np.zeros((0,0)))
        xp_assert_equal(lut([], [], grid=False), np.zeros((0,)))

    def test_invalid_input(self):
        s = 0.1
        tx = [1 + s, 3 - s]
        ty = [1 + s, 3 - s]

        with assert_raises(ValueError) as info:
            x = np.linspace(1.0, 10.0)
            y = np.linspace(1.0, 10.0)
            z = np.linspace(1.0, 10.0, num=10)
            LSQBivariateSpline(x, y, z, tx, ty)
        assert "x, y, and z should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            x = np.linspace(1.0, 10.0)
            y = np.linspace(1.0, 10.0)
            z = np.linspace(1.0, 10.0)
            w = np.linspace(1.0, 10.0, num=20)
            LSQBivariateSpline(x, y, z, tx, ty, w=w)
        assert "x, y, z, and w should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            w = np.linspace(-1.0, 10.0)
            LSQBivariateSpline(x, y, z, tx, ty, w=w)
        assert "w should be positive" in str(info.value)

        with assert_raises(ValueError) as info:
            bbox = (-100, 100, -100)
            LSQBivariateSpline(x, y, z, tx, ty, bbox=bbox)
        assert "bbox shape should be (4,)" in str(info.value)

        with assert_raises(ValueError) as info:
            LSQBivariateSpline(x, y, z, tx, ty, kx=10, ky=10)
        assert "The length of x, y and z should be at least (kx+1) * (ky+1)" in \
               str(info.value)

        with assert_raises(ValueError) as exc_info:
            LSQBivariateSpline(x, y, z, tx, ty, eps=0.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            LSQBivariateSpline(x, y, z, tx, ty, eps=1.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

    def test_array_like_input(self):
        s = 0.1
        tx = np.array([1 + s, 3 - s])
        ty = np.array([1 + s, 3 - s])
        x = np.linspace(1.0, 10.0)
        y = np.linspace(1.0, 10.0)
        z = np.linspace(1.0, 10.0)
        w = np.linspace(1.0, 10.0)
        bbox = np.array([1.0, 10.0, 1.0, 10.0])

        with pytest.warns(UserWarning, match="\nThe coefficients of the spline") as r:
            # np.array input
            spl1 = LSQBivariateSpline(x, y, z, tx, ty, w=w, bbox=bbox)
            # list input
            spl2 = LSQBivariateSpline(x.tolist(), y.tolist(), z.tolist(),
                                      tx.tolist(), ty.tolist(), w=w.tolist(),
                                      bbox=bbox)
            xp_assert_close(spl1(2.0, 2.0), spl2(2.0, 2.0))
            assert len(r) == 2

    def test_unequal_length_of_knots(self):
        """Test for the case when the input knot-location arrays in x and y are
        of different lengths.
        """
        x, y = np.mgrid[0:100, 0:100]
        x = x.ravel()
        y = y.ravel()
        z = 3.0 * np.ones_like(x)
        tx = np.linspace(0.1, 98.0, 29)
        ty = np.linspace(0.1, 98.0, 33)
        with pytest.warns(UserWarning, match="\nThe coefficients of the spline") as r:
            lut = LSQBivariateSpline(x,y,z,tx,ty)
            assert len(r) == 1

        assert_almost_equal(lut(x, y, grid=False), z)


class TestSmoothBivariateSpline:
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        for t in lut.get_knots():
            assert_array_almost_equal(t, [1, 1, 3, 3])

        assert_array_almost_equal(lut.get_coeffs(), [3, 3, 3, 3])
        assert abs(lut.get_residual()) < 1e-15
        assert_array_almost_equal(lut([1, 1.5, 2], [1, 1.5]), [[3, 3], [3, 3], [3, 3]])

    def test_linear_1d(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,0,0,2,2,2,4,4,4]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        for t in lut.get_knots():
            xp_assert_close(t, np.asarray([1.0, 1, 3, 3]))
        assert_array_almost_equal(lut.get_coeffs(), [0, 0, 4, 4])
        assert abs(lut.get_residual()) < 1e-15
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[0,0],[1,1],[2,2]])

    def test_integral(self):
        x = [1,1,1,2,2,2,4,4,4]
        y = [1,2,3,1,2,3,1,2,3]
        z = array([0,7,8,3,4,7,1,3,4])

        with warnings.catch_warnings():
            # This seems to fail (ier=1, see ticket 1642).
            warnings.filterwarnings(
                "ignore", "\nThe required storage space", UserWarning)
            lut = SmoothBivariateSpline(x, y, z, kx=1, ky=1, s=0)

        tx = [1,2,4]
        ty = [1,2,3]

        tz = lut(tx, ty)
        trpz = .25*(diff(tx)[:,None]*diff(ty)[None,:]
                    * (tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()
        assert_almost_equal(np.asarray(lut.integral(tx[0], tx[-1], ty[0], ty[-1])),
                            np.asarray(trpz))

        lut2 = SmoothBivariateSpline(x, y, z, kx=2, ky=2, s=0)
        assert_almost_equal(np.asarray(lut2.integral(tx[0], tx[-1], ty[0], ty[-1])),
                            np.asarray(trpz),
                            decimal=0)  # the quadratures give 23.75 and 23.85

        tz = lut(tx[:-1], ty[:-1])
        trpz = .25*(diff(tx[:-1])[:,None]*diff(ty[:-1])[None,:]
                    * (tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()
        assert_almost_equal(np.asarray(lut.integral(tx[0], tx[-2], ty[0], ty[-2])),
                            np.asarray(trpz))

    def test_rerun_lwrk2_too_small(self):
        # in this setting, lwrk2 is too small in the default run. Here we
        # check for equality with the bisplrep/bisplev output because there,
        # an automatic re-run of the spline representation is done if ier>10.
        x = np.linspace(-2, 2, 80)
        y = np.linspace(-2, 2, 80)
        z = x + y
        xi = np.linspace(-1, 1, 100)
        yi = np.linspace(-2, 2, 100)
        tck = bisplrep(x, y, z)
        res1 = bisplev(xi, yi, tck)
        interp_ = SmoothBivariateSpline(x, y, z)
        res2 = interp_(xi, yi)
        assert_almost_equal(res1, res2)

    def test_invalid_input(self):

        with assert_raises(ValueError) as info:
            x = np.linspace(1.0, 10.0)
            y = np.linspace(1.0, 10.0)
            z = np.linspace(1.0, 10.0, num=10)
            SmoothBivariateSpline(x, y, z)
        assert "x, y, and z should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            x = np.linspace(1.0, 10.0)
            y = np.linspace(1.0, 10.0)
            z = np.linspace(1.0, 10.0)
            w = np.linspace(1.0, 10.0, num=20)
            SmoothBivariateSpline(x, y, z, w=w)
        assert "x, y, z, and w should have a same length" in str(info.value)

        with assert_raises(ValueError) as info:
            w = np.linspace(-1.0, 10.0)
            SmoothBivariateSpline(x, y, z, w=w)
        assert "w should be positive" in str(info.value)

        with assert_raises(ValueError) as info:
            bbox = (-100, 100, -100)
            SmoothBivariateSpline(x, y, z, bbox=bbox)
        assert "bbox shape should be (4,)" in str(info.value)

        with assert_raises(ValueError) as info:
            SmoothBivariateSpline(x, y, z, kx=10, ky=10)
        assert "The length of x, y and z should be at least (kx+1) * (ky+1)" in\
               str(info.value)

        with assert_raises(ValueError) as info:
            SmoothBivariateSpline(x, y, z, s=-1.0)
        assert "s should be s >= 0.0" in str(info.value)

        with assert_raises(ValueError) as exc_info:
            SmoothBivariateSpline(x, y, z, eps=0.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            SmoothBivariateSpline(x, y, z, eps=1.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

    def test_array_like_input(self):
        x = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3])
        y = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3])
        z = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3])
        w = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])
        bbox = np.array([1.0, 3.0, 1.0, 3.0])
        # np.array input
        spl1 = SmoothBivariateSpline(x, y, z, w=w, bbox=bbox, kx=1, ky=1)
        # list input
        spl2 = SmoothBivariateSpline(x.tolist(), y.tolist(), z.tolist(),
                                     bbox=bbox.tolist(), w=w.tolist(),
                                     kx=1, ky=1)
        xp_assert_close(spl1(0.1, 0.5), spl2(0.1, 0.5))


class TestLSQSphereBivariateSpline:
    def setup_method(self):
        # define the input data and coordinates
        ntheta, nphi = 70, 90
        theta = linspace(0.5/(ntheta - 1), 1 - 0.5/(ntheta - 1), ntheta) * pi
        phi = linspace(0.5/(nphi - 1), 1 - 0.5/(nphi - 1), nphi) * 2. * pi
        data = ones((theta.shape[0], phi.shape[0]))
        # define knots and extract data values at the knots
        knotst = theta[::5]
        knotsp = phi[::5]
        knotdata = data[::5, ::5]
        # calculate spline coefficients
        lats, lons = meshgrid(theta, phi)
        lut_lsq = LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                           data.T.ravel(), knotst, knotsp)
        self.lut_lsq = lut_lsq
        self.data = knotdata
        self.new_lons, self.new_lats = knotsp, knotst

    def test_linear_constant(self):
        assert abs(self.lut_lsq.get_residual()) < 1e-15
        assert_array_almost_equal(self.lut_lsq(self.new_lats, self.new_lons),
                                  self.data)

    def test_empty_input(self):
        assert_array_almost_equal(self.lut_lsq([], []), np.zeros((0,0)))
        assert_array_almost_equal(self.lut_lsq([], [], grid=False), np.zeros((0,)))

    def test_invalid_input(self):
        ntheta, nphi = 70, 90
        theta = linspace(0.5 / (ntheta - 1), 1 - 0.5 / (ntheta - 1),
                         ntheta) * pi
        phi = linspace(0.5 / (nphi - 1), 1 - 0.5 / (nphi - 1), nphi) * 2. * pi
        data = ones((theta.shape[0], phi.shape[0]))
        # define knots and extract data values at the knots
        knotst = theta[::5]
        knotsp = phi[::5]

        with assert_raises(ValueError) as exc_info:
            invalid_theta = linspace(-0.1, 1.0, num=ntheta) * pi
            invalid_lats, lons = meshgrid(invalid_theta, phi)
            LSQSphereBivariateSpline(invalid_lats.ravel(), lons.ravel(),
                                     data.T.ravel(), knotst, knotsp)
        assert "theta should be between [0, pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_theta = linspace(0.1, 1.1, num=ntheta) * pi
            invalid_lats, lons = meshgrid(invalid_theta, phi)
            LSQSphereBivariateSpline(invalid_lats.ravel(), lons.ravel(),
                                     data.T.ravel(), knotst, knotsp)
        assert "theta should be between [0, pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_phi = linspace(-0.1, 1.0, num=ntheta) * 2.0 * pi
            lats, invalid_lons = meshgrid(theta, invalid_phi)
            LSQSphereBivariateSpline(lats.ravel(), invalid_lons.ravel(),
                                     data.T.ravel(), knotst, knotsp)
        assert "phi should be between [0, 2pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_phi = linspace(0.0, 1.1, num=ntheta) * 2.0 * pi
            lats, invalid_lons = meshgrid(theta, invalid_phi)
            LSQSphereBivariateSpline(lats.ravel(), invalid_lons.ravel(),
                                     data.T.ravel(), knotst, knotsp)
        assert "phi should be between [0, 2pi]" in str(exc_info.value)

        lats, lons = meshgrid(theta, phi)

        with assert_raises(ValueError) as exc_info:
            invalid_knotst = np.copy(knotst)
            invalid_knotst[0] = -0.1
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                     data.T.ravel(), invalid_knotst, knotsp)
        assert "tt should be between (0, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_knotst = np.copy(knotst)
            invalid_knotst[0] = pi
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                     data.T.ravel(), invalid_knotst, knotsp)
        assert "tt should be between (0, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_knotsp = np.copy(knotsp)
            invalid_knotsp[0] = -0.1
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                     data.T.ravel(), knotst, invalid_knotsp)
        assert "tp should be between (0, 2pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_knotsp = np.copy(knotsp)
            invalid_knotsp[0] = 2 * pi
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                     data.T.ravel(), knotst, invalid_knotsp)
        assert "tp should be between (0, 2pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_w = array([-1.0, 1.0, 1.5, 0.5, 1.0, 1.5, 0.5, 1.0, 1.0])
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(), data.T.ravel(),
                                     knotst, knotsp, w=invalid_w)
        assert "w should be positive" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(), data.T.ravel(),
                                     knotst, knotsp, eps=0.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            LSQSphereBivariateSpline(lats.ravel(), lons.ravel(), data.T.ravel(),
                                     knotst, knotsp, eps=1.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

    def test_array_like_input(self):
        ntheta, nphi = 70, 90
        theta = linspace(0.5 / (ntheta - 1), 1 - 0.5 / (ntheta - 1),
                         ntheta) * pi
        phi = linspace(0.5 / (nphi - 1), 1 - 0.5 / (nphi - 1),
                       nphi) * 2. * pi
        lats, lons = meshgrid(theta, phi)
        data = ones((theta.shape[0], phi.shape[0]))
        # define knots and extract data values at the knots
        knotst = theta[::5]
        knotsp = phi[::5]
        w = ones(lats.ravel().shape[0])

        # np.array input
        spl1 = LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                        data.T.ravel(), knotst, knotsp, w=w)
        # list input
        spl2 = LSQSphereBivariateSpline(lats.ravel().tolist(),
                                        lons.ravel().tolist(),
                                        data.T.ravel().tolist(),
                                        knotst.tolist(),
                                        knotsp.tolist(), w=w.tolist())
        assert_array_almost_equal(spl1(1.0, 1.0), spl2(1.0, 1.0))


class TestSmoothSphereBivariateSpline:
    def setup_method(self):
        theta = array([.25*pi, .25*pi, .25*pi, .5*pi, .5*pi, .5*pi, .75*pi,
                       .75*pi, .75*pi])
        phi = array([.5 * pi, pi, 1.5 * pi, .5 * pi, pi, 1.5 * pi, .5 * pi, pi,
                     1.5 * pi])
        r = array([3, 3, 3, 3, 3, 3, 3, 3, 3])
        self.lut = SmoothSphereBivariateSpline(theta, phi, r, s=1E10)

    def test_linear_constant(self):
        assert abs(self.lut.get_residual()) < 1e-15
        assert_array_almost_equal(self.lut([1, 1.5, 2],[1, 1.5]),
                                  [[3, 3], [3, 3], [3, 3]])

    def test_empty_input(self):
        assert_array_almost_equal(self.lut([], []), np.zeros((0,0)))
        assert_array_almost_equal(self.lut([], [], grid=False), np.zeros((0,)))

    def test_invalid_input(self):
        theta = array([.25 * pi, .25 * pi, .25 * pi, .5 * pi, .5 * pi, .5 * pi,
                       .75 * pi, .75 * pi, .75 * pi])
        phi = array([.5 * pi, pi, 1.5 * pi, .5 * pi, pi, 1.5 * pi, .5 * pi, pi,
                     1.5 * pi])
        r = array([3, 3, 3, 3, 3, 3, 3, 3, 3])

        with assert_raises(ValueError) as exc_info:
            invalid_theta = array([-0.1 * pi, .25 * pi, .25 * pi, .5 * pi,
                                   .5 * pi, .5 * pi, .75 * pi, .75 * pi,
                                   .75 * pi])
            SmoothSphereBivariateSpline(invalid_theta, phi, r, s=1E10)
        assert "theta should be between [0, pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_theta = array([.25 * pi, .25 * pi, .25 * pi, .5 * pi,
                                   .5 * pi, .5 * pi, .75 * pi, .75 * pi,
                                   1.1 * pi])
            SmoothSphereBivariateSpline(invalid_theta, phi, r, s=1E10)
        assert "theta should be between [0, pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_phi = array([-.1 * pi, pi, 1.5 * pi, .5 * pi, pi, 1.5 * pi,
                                 .5 * pi, pi, 1.5 * pi])
            SmoothSphereBivariateSpline(theta, invalid_phi, r, s=1E10)
        assert "phi should be between [0, 2pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_phi = array([1.0 * pi, pi, 1.5 * pi, .5 * pi, pi, 1.5 * pi,
                                 .5 * pi, pi, 2.1 * pi])
            SmoothSphereBivariateSpline(theta, invalid_phi, r, s=1E10)
        assert "phi should be between [0, 2pi]" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            invalid_w = array([-1.0, 1.0, 1.5, 0.5, 1.0, 1.5, 0.5, 1.0, 1.0])
            SmoothSphereBivariateSpline(theta, phi, r, w=invalid_w, s=1E10)
        assert "w should be positive" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            SmoothSphereBivariateSpline(theta, phi, r, s=-1.0)
        assert "s should be positive" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            SmoothSphereBivariateSpline(theta, phi, r, eps=-1.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            SmoothSphereBivariateSpline(theta, phi, r, eps=1.0)
        assert "eps should be between (0, 1)" in str(exc_info.value)

    def test_array_like_input(self):
        theta = np.array([.25 * pi, .25 * pi, .25 * pi, .5 * pi, .5 * pi,
                          .5 * pi, .75 * pi, .75 * pi, .75 * pi])
        phi = np.array([.5 * pi, pi, 1.5 * pi, .5 * pi, pi, 1.5 * pi, .5 * pi,
                        pi, 1.5 * pi])
        r = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3])
        w = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

        # np.array input
        spl1 = SmoothSphereBivariateSpline(theta, phi, r, w=w, s=1E10)

        # list input
        spl2 = SmoothSphereBivariateSpline(theta.tolist(), phi.tolist(),
                                           r.tolist(), w=w.tolist(), s=1E10)
        assert_array_almost_equal(spl1(1.0, 1.0), spl2(1.0, 1.0))


class TestRectBivariateSpline:
    @pytest.mark.parametrize("k", [3, 4])
    def test_defaults(self, k):
        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1],
                   [1, 2, 1, 2, 1],
                   [1, 2, 3, 2, 1],
                   [1, 2, 2, 2, 1],
                   [1, 2, 1, 2, 1]])

        lut = RectBivariateSpline(x, y, z, kx=k, ky=k)
        assert_array_almost_equal(lut(x, y), z)

        lut_custom = _regrid(x, y, z, kx=k, ky=k)
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom, x, y), z.astype(np.float64))

    @pytest.mark.parametrize("k", [3, 4])
    def test_interpolated_offgrid_points(self, k):
        # Setup
        x = np.linspace(0, 4, 5)
        y = np.linspace(0, 4, 5)
        z = np.sin(x[:, None]) * np.cos(y[None, :])

        xi = np.linspace(0, 4, 50)  # much finer grid
        yi = np.linspace(0, 4, 50)

        ref = RectBivariateSpline(x, y, z, kx=k, ky=k)
        ref_vals = ref(xi, yi)

        custom = _regrid(x, y, z, kx=k, ky=k)
        custom_vals = _ndbspline_call_like_bivariate(custom, xi, yi)
        # Compare interpolated values
        xp_assert_close(custom_vals, ref_vals,
                        rtol=0, atol=1.5*10**(-6))

    @pytest.mark.parametrize("k", [3, 4])
    def test_midpoints(self, k):
        # Midpoint interpolation
        x = np.arange(1, k + 2)
        y = np.arange(1, k + 2)
        z = x[:, None] + y[None, :]  # simple known function

        xi = (x[:-1] + x[1:]) / 2
        yi = (y[:-1] + y[1:]) / 2

        ref = RectBivariateSpline(x, y, z, kx=k, ky=k)
        custom = _regrid(x, y, z, kx=k, ky=k)

        ref_vals = ref(xi, yi)
        custom_vals = _ndbspline_call_like_bivariate(custom, xi, yi)

        xp_assert_close(custom_vals, ref_vals,
                        rtol=0, atol=1.5*10**(-6))

    @pytest.mark.parametrize("s_val", [0.1, 1.0, 10.0, 100.0])
    @pytest.mark.parametrize("k", [3, 4])
    def test_custom_matches_reference(self, s_val, k):
        x = np.linspace(0, 4, 20)
        y = np.linspace(0, 4, 20)
        X, Y = np.meshgrid(x, y, indexing="ij")
        z = X + Y

        x_mid = (x[:-1] + x[1:]) / 2
        y_mid = (y[:-1] + y[1:]) / 2
        xi_dense = np.linspace(0, 4, 50)
        yi_dense = np.linspace(0, 4, 50)

        query_sets = [
            ("original grid", x, y),
            ("midpoints", x_mid, y_mid),
            ("dense grid", xi_dense, yi_dense),
        ]

        ref = RectBivariateSpline(x, y, z, s=s_val, kx=k, ky=k)
        custom = _regrid(x, y, z, s=s_val, kx=k, ky=k)

        for _, xi, yi in query_sets:
            ref_vals = ref(xi, yi)
            custom_vals = _ndbspline_call_like_bivariate(custom, xi, yi)

            xp_assert_close(custom_vals, ref_vals,
                            rtol=2e-2, atol=2e-2)

    @pytest.mark.parametrize("k", [3, 4])
    def test_evaluate(self, k):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z,kx=k,ky=k)
        lut_custom = _regrid(x, y, z, kx=k, ky=k)

        xi = [1, 2.3, 5.3, 0.5, 3.3, 1.2, 3]
        yi = [1, 3.3, 1.2, 4.0, 5.0, 1.0, 3]
        zi = lut.ev(xi, yi)
        zi2 = array([lut(xp, yp)[0,0] for xp, yp in zip(xi, yi)])

        assert_almost_equal(zi, zi2)

        zi_custom = _ndbspline_call_like_bivariate(lut_custom, xi, yi, grid=False)
        zi2_custom = array([_ndbspline_call_like_bivariate(lut_custom, xp, yp)[0, 0]
                            for xp, yp in zip(xi, yi)])
        xp_assert_close(zi_custom, zi2_custom,
                        rtol=0, atol=1.5*10**(-6))

    def test_derivatives_grid(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        dx = array([[0,0,-20,0,0],[0,0,13,0,0],[0,0,4,0,0],
            [0,0,-11,0,0],[0,0,4,0,0]])/6.
        dy = array([[4,-1,0,1,-4],[4,-1,0,1,-4],[0,1.5,0,-1.5,0],
            [2,.25,0,-.25,-2],[4,-1,0,1,-4]])
        dxdy = array([[40,-25,0,25,-40],[-26,16.25,0,-16.25,26],
            [-8,5,0,-5,8],[22,-13.75,0,13.75,-22],[-8,5,0,-5,8]])/6.

        lut = RectBivariateSpline(x,y,z)
        lut_custom = _regrid(x,y,z)

        assert_array_almost_equal(lut(x,y,dx=1),dx)
        assert_array_almost_equal(lut(x,y,dy=1),dy)
        assert_array_almost_equal(lut(x,y,dx=1,dy=1),dxdy)
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom, x,y,dx=1),dx,
            rtol=0, atol=1.5*10**(-6))
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom, x,y,dy=1),dy,
            rtol=0, atol=1.5*10**(-6))
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom, x,y,dx=1,dy=1),dxdy,
            rtol=0, atol=1.5*10**(-6))

    def test_derivatives(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        dx = array([0,0,2./3,0,0])
        dy = array([4,-1,0,-.25,-4])
        dxdy = array([160,65,0,55,32])/24.
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y,dx=1,grid=False),dx)
        assert_array_almost_equal(lut(x,y,dy=1,grid=False),dy)
        assert_array_almost_equal(lut(x,y,dx=1,dy=1,grid=False),dxdy)

        lut_custom = _regrid(x,y,z)
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom,x,y,dx=1,grid=False),
            dx,rtol=0, atol=1.5*10**(-6)
        )
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom,x,y,dy=1,grid=False),
            dy,rtol=0, atol=1.5*10**(-6)
        )
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom,x,y,dx=1,dy=1,grid=False),
            dxdy,rtol=0, atol=1.5*10**(-6)
        )

    def make_pair_grid(self, x, y):
        """
        Create an array of (xi, yi) pairs for all xi in x and yi in y,
        and reshape it to the desired shape.

        Parameters
        ----------
        x : array_like
            1D array of x-values.
        y : array_like
            1D array of y-values.
        dest_shape : tuple
            Desired output shape.

        Returns
        -------
        np.ndarray
            Reshaped array of (x, y) pairs.
        """
        return np.array([[xi, yi] for xi in x for yi in y])

    def test_partial_derivative_method_grid(self):
        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1],
                   [1, 2, 1, 2, 1],
                   [1, 2, 3, 2, 1],
                   [1, 2, 2, 2, 1],
                   [1, 2, 1, 2, 1]])
        dx = array([[0, 0, -20, 0, 0],
                    [0, 0, 13, 0, 0],
                    [0, 0, 4, 0, 0],
                    [0, 0, -11, 0, 0],
                    [0, 0, 4, 0, 0]]) / 6.
        dy = array([[4, -1, 0, 1, -4],
                    [4, -1, 0, 1, -4],
                    [0, 1.5, 0, -1.5, 0],
                    [2, .25, 0, -.25, -2],
                    [4, -1, 0, 1, -4]])
        dxdy = array([[40, -25, 0, 25, -40],
                      [-26, 16.25, 0, -16.25, 26],
                      [-8, 5, 0, -5, 8],
                      [22, -13.75, 0, 13.75, -22],
                      [-8, 5, 0, -5, 8]]) / 6.
        lut = RectBivariateSpline(x, y, z)
        lut_ndbspline = convert_to_ndbspline(lut)
        for orders, expected in [([1, 0], dx), ([0, 1], dy), ([1, 1], dxdy)]:
            actual_rect = lut.partial_derivative(*orders)(x, y)
            actual_ndb = lut_ndbspline.derivative(orders)(
                self.make_pair_grid(x, y)
            ).reshape(expected.shape)

            assert_array_almost_equal(actual_rect, expected)
            assert_array_almost_equal(actual_ndb, expected)

        lut_custom = _regrid(x, y, z)
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom.derivative([1, 0]), x, y),
            dx,rtol=0, atol=1.5*10**(-6)
        )
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom.derivative([0, 1]), x, y),
            dy,rtol=0, atol=1.5*10**(-6)
        )
        xp_assert_close(
            _ndbspline_call_like_bivariate(lut_custom.derivative([1, 1]), x, y),
            dxdy,rtol=0, atol=1.5*10**(-6)
        )

    def test_partial_derivative_method(self):
        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1],
                   [1, 2, 1, 2, 1],
                   [1, 2, 3, 2, 1],
                   [1, 2, 2, 2, 1],
                   [1, 2, 1, 2, 1]])
        expected = {
            (1, 0): array([0, 0, 2./3, 0, 0]), # dx
            (0, 1): array([4, -1, 0, -.25, -4]), # dy
            (1, 1): array([160, 65, 0, 55, 32]) / 24. # dxdy
        }

        lut = RectBivariateSpline(x, y, z)
        lut_ndbspline = convert_to_ndbspline(lut)
        lut_custom = _regrid(x, y, z)

        points = self.make_pair_grid(x, y)  # shape: (25, 2)

        # Evaluate only the diagonal points: (x[i], y[i])
        diag_idx = np.arange(len(x))
        diag_points = points[diag_idx * len(y) + diag_idx]

        for orders, expected_vals in expected.items():
            dx, dy = orders
            # RectBivariateSpline result
            actual_rbs = lut.partial_derivative(dx, dy)(x, y, grid=False)
            assert_array_almost_equal(actual_rbs, expected_vals)

            # NdBSpline result
            actual_ndb = lut_ndbspline.derivative([dx, dy])(diag_points)
            assert_array_almost_equal(actual_ndb, expected_vals)

            # regrid result
            actual_rp = _ndbspline_call_like_bivariate(lut_custom.derivative([dx, dy]),
                x, y, grid=False)
            xp_assert_close(actual_rp, expected_vals,
                            rtol=0, atol=1.5*10**(-6))

    @pytest.mark.parametrize("k", [3, 4])
    def test_partial_derivative_order_too_large(self, k):
        x = array([0, 1, 2, 3, 4], dtype=float)
        y = x.copy()
        z = ones((x.size, y.size))
        lut = RectBivariateSpline(x, y, z)
        lut_ndbspline = convert_to_ndbspline(lut)
        with assert_raises(ValueError):
            lut.partial_derivative(k + 1, 1)

        assert (lut_ndbspline.derivative([k + 1, 1]).c == 0.0).all()

        lut = RectBivariateSpline(x, y, z, kx=k, ky=k)
        with assert_raises(ValueError):
            lut.partial_derivative(k + 1, 1)

        lut_custom = _regrid(x, y, z, kx=k, ky=k)
        assert (lut_custom.derivative([k + 1, 1]).c == 0.0).all()

    @pytest.mark.parametrize("k", [3, 4])
    def test_broadcast(self, k):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z,kx=k,ky=k)
        xp_assert_close(lut(x, y), lut(x[:,None], y[None,:], grid=False))

    def test_invalid_input(self):

        x = array([6, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                    [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])

        with assert_raises(ValueError) as info:
            RectBivariateSpline(x, y, z)
        assert "x must be strictly increasing" in str(info.value)

        with assert_raises(ValueError) as info_custom:
            _regrid(x, y, z)
        assert "x must be strictly increasing" in str(info_custom.value)

        x = array([1, 2, 3, 4, 5])
        y = array([2, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                    [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])
        with assert_raises(ValueError) as info:
            RectBivariateSpline(x, y, z)
        assert "y must be strictly increasing" in str(info.value)

        with assert_raises(ValueError) as info_custom:
            _regrid(x, y, z)
        assert "y must be strictly increasing" in str(info_custom.value)

        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                    [1, 2, 2, 2, 1]])
        with assert_raises(ValueError) as info:
            RectBivariateSpline(x, y, z)
        assert "x dimension of z must have same number of elements as x"\
               in str(info.value)

        with assert_raises(ValueError) as info_custom:
            _regrid(x, y, z)
        assert "x dimension of z must have same number of elements as x"\
               in str(info_custom.value)

        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2], [1, 2, 1, 2], [1, 2, 3, 2],
                    [1, 2, 2, 2], [1, 2, 1, 2]])

        with assert_raises(ValueError) as info:
            RectBivariateSpline(x, y, z)
        assert "y dimension of z must have same number of elements as y"\
               in str(info.value)

        with assert_raises(ValueError) as info_custom:
            _regrid(x, y, z)
        assert "y dimension of z must have same number of elements as y"\
               in str(info_custom.value)

        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                    [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])
        bbox = (-100, 100, -100)

        with assert_raises(ValueError) as info:
            RectBivariateSpline(x, y, z, bbox=bbox)
        assert "bbox shape should be (4,)" in str(info.value)

        with assert_raises(ValueError) as info_custom:
            _regrid(x, y, z, bbox=bbox)
        assert "bbox shape should be (4,)" in str(info_custom.value)

        with assert_raises(ValueError) as info:
            RectBivariateSpline(x, y, z, s=-1.0)
        assert "s should be s >= 0.0" in str(info.value)

        with assert_raises(ValueError) as info_custom:
            _regrid(x, y, z, s=-1.0)
        assert "s should be s >= 0.0" in str(info_custom.value)

    @pytest.mark.parametrize("k", [3, 4])
    def test_array_like_input(self, k):
        x = array([1, 2, 3, 4, 5])
        y = array([1, 2, 3, 4, 5])
        z = array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                   [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])
        bbox = array([1, 5, 1, 5])

        spl1 = RectBivariateSpline(x, y, z, bbox=bbox, kx=k, ky=k)
        spl1_custom = _regrid(x, y, z, bbox=bbox, kx=k, ky=k)
        spl2 = RectBivariateSpline(x.tolist(), y.tolist(), z.tolist(),
                                   bbox=bbox.tolist(), kx=k, ky=k)
        spl2_custom = _regrid(x.tolist(), y.tolist(), z.tolist(),
                                   bbox=bbox.tolist(), kx=k, ky=k)
        assert_array_almost_equal(spl1(1.0, 1.0), spl2(1.0, 1.0))
        xp_assert_close(spl1(1.0, 1.0),
                        _ndbspline_call_like_bivariate(spl1_custom, 1.0, 1.0),
                        rtol=0, atol=1.5*10**(-6))
        xp_assert_close(spl2(1.0, 1.0),
                        _ndbspline_call_like_bivariate(spl2_custom, 1.0, 1.0),
                        rtol=0, atol=1.5*10**(-6))
        xp_assert_close(_ndbspline_call_like_bivariate(spl1_custom, 1.0, 1.0),
                        _ndbspline_call_like_bivariate(spl2_custom, 1.0, 1.0),
                        rtol=0, atol=1.5*10**(-6))

    def test_not_increasing_input(self):
        # gh-8565
        NSamp = 20
        Theta = np.random.uniform(0, np.pi, NSamp)
        Phi = np.random.uniform(0, 2 * np.pi, NSamp)
        Data = np.ones(NSamp)

        Interpolator = SmoothSphereBivariateSpline(Theta, Phi, Data, s=3.5)

        NLon = 6
        NLat = 3
        GridPosLats = np.arange(NLat) / NLat * np.pi
        GridPosLons = np.arange(NLon) / NLon * 2 * np.pi

        # No error
        Interpolator(GridPosLats, GridPosLons)

        nonGridPosLats = GridPosLats.copy()
        nonGridPosLats[2] = 0.001
        with assert_raises(ValueError) as exc_info:
            Interpolator(nonGridPosLats, GridPosLons)
        assert "x must be strictly increasing" in str(exc_info.value)

        nonGridPosLons = GridPosLons.copy()
        nonGridPosLons[2] = 0.001
        with assert_raises(ValueError) as exc_info:
            Interpolator(GridPosLats, nonGridPosLons)
        assert "y must be strictly increasing" in str(exc_info.value)

    def _sample_large_2d_data(self, nx, ny):
        rng = np.random.default_rng(1)
        x = np.arange(nx)
        y = np.arange(ny)
        z = rng.integers(0, 100, (nx, ny))

        return x, y, z.astype(np.float64)

    # Thin wrapper used as the ``RectBivariateSpline`` entry in the
    # ``spl_apis`` parametrize dimension of ``test_spline_large_2d`` and
    # ``test_spline_large_2d_maxit``.  A standalone callable is required
    # because ``RectBivariateSpline.__call__`` cannot be passed directly as a
    # parametrize value alongside ``_ndbspline_call_like_bivariate``.
    def _RectBivariateSplineEval(spl, x, y):
        return spl(x, y)

    @pytest.mark.slow()
    @pytest.mark.parametrize('shape', [(350, 850), (2000, 170)])
    @pytest.mark.parametrize('s_tols', [(0, 1e-12, 1e-7),
                                        (1, 7e-3, 1e-4),
                                        (3, 2e-2, 1e-4)])
    @pytest.mark.parametrize("spl_apis", ((RectBivariateSpline,
                                           _RectBivariateSplineEval),
                                          (_regrid,
                                           _ndbspline_call_like_bivariate)))
    def test_spline_large_2d(self, shape, s_tols, spl_apis):
        # Reference - https://github.com/scipy/scipy/issues/17787
        #
        # ``spl_apis`` is intentionally a parametrize dimension
        # so that pytest-timeout's ``--timeout`` limit is applied
        # individually to each (API, shape, s_tols) combination.  If both APIs
        # were called sequentially inside a single test function they would
        # share one timeout budget, making it much easier to exceed the limit
        # on slow platforms (e.g. Windows CI) and triggering the
        # ``--timeout-method=thread`` worker crash described in gh-25142.
        nx, ny = shape
        s, atol, rtol = s_tols
        x, y, z = self._sample_large_2d_data(nx, ny)

        spl_construct, spl_eval = spl_apis

        spl = spl_construct(x, y, z, s=s)
        z_spl = spl_eval(spl, x, y)
        assert(not np.isnan(z_spl).any())
        xp_assert_close(z_spl, z, atol=atol, rtol=rtol)

    @pytest.mark.xslow()
    @pytest.mark.skipif(sys.maxsize <= 2**32, reason="Segfaults on 32-bit system "
                                                     "due to large input data")
    @pytest.mark.parametrize("k", [3, 4])
    @pytest.mark.parametrize("spl_apis", ((RectBivariateSpline,
                                           _RectBivariateSplineEval),
                                          (_regrid,
                                           _ndbspline_call_like_bivariate)))
    def test_spline_large_2d_maxit(self, k, spl_apis):
        # Reference - for https://github.com/scipy/scipy/issues/17787
        #
        # ``spl_apis`` is intentionally a parametrize dimension
        # for the same reason as in ``test_spline_large_2d``: so
        # that pytest-timeout's ``--timeout`` limit is applied individually to
        # each (API, k) combination instead of being shared across both APIs
        # in a single call, which would risk a ``--timeout-method=thread``
        # worker crash on Windows CI (see gh-25142).
        nx, ny = 1000, 1700
        s, atol, rtol = 2, 2e-2, 1e-12
        x, y, z = self._sample_large_2d_data(nx, ny)

        spl_construct, spl_eval = spl_apis

        spl = spl_construct(x, y, z, s=s,
                            maxit=30, kx=k, ky=k)
        z_spl = spl_eval(spl, x, y)
        assert(not np.isnan(z_spl).any())
        xp_assert_close(z_spl, z, atol=atol, rtol=rtol)

    @pytest.mark.parametrize("spl_apis", ((RectBivariateSpline,
                                           _RectBivariateSplineEval),
                                          (_regrid,
                                           _ndbspline_call_like_bivariate)))
    def test_spline_synthetic_data(self, spl_apis):
        """
        Test regrid with synthetic smooth data (mixed frequencies + noise).
        """
        # ``spl_apis`` is parametrized for consistency with
        # ``test_spline_large_2d`` and ``test_spline_large_2d_maxit``,
        # ensuring both APIs are exercised under the same test structure.
        #
        # Create strictly-increasing axes and smooth test surface
        nx, ny = 64, 64
        kx, ky = 3, 3
        s = 0.3

        rng = np.random.default_rng(0)
        x = np.linspace(0.0, 4.0, nx, dtype=float)
        y = np.linspace(0.0, 4.0, ny, dtype=float)

        # Smooth, nontrivial surface with mixed frequencies + mild noise
        X, Y = np.meshgrid(x, y, indexing="ij")
        z = (
            np.sin(X) * np.cos(Y)
            + 0.2 * np.sin(2 * X + 0.5) * np.cos(1.5 * Y - 0.3)
            + 0.05 * rng.normal(size=(nx, ny))
        ).astype(float)

        spl_construct, spl_eval = spl_apis

        spl = spl_construct(x, y, z, kx=kx, ky=ky, s=s)
        z_spl = spl_eval(spl, x, y)
        assert not np.isnan(z_spl).any()
        xp_assert_close(z_spl, z, atol=0.1, rtol=0.1)


class TestRectSphereBivariateSpline:
    def test_defaults(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])
        lut = RectSphereBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y),z)

    def test_evaluate(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])
        lut = RectSphereBivariateSpline(x,y,z)
        yi = [0.2, 1, 2.3, 2.35, 3.0, 3.99, 5.25]
        xi = [1.5, 0.4, 1.1, 0.45, 0.2345, 1., 0.0001]
        zi = lut.ev(xi, yi)
        zi2 = array([lut(xp, yp)[0,0] for xp, yp in zip(xi, yi)])
        assert_almost_equal(zi, zi2)

    def test_invalid_input(self):
        data = np.dot(np.atleast_2d(90. - np.linspace(-80., 80., 18)).T,
                      np.atleast_2d(180. - np.abs(np.linspace(0., 350., 9)))).T

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(-1, 170, 9) * np.pi / 180.
            lons = np.linspace(0, 350, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "u should be between (0, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 181, 9) * np.pi / 180.
            lons = np.linspace(0, 350, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "u should be between (0, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 170, 9) * np.pi / 180.
            lons = np.linspace(-181, 10, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "v[0] should be between [-pi, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 170, 9) * np.pi / 180.
            lons = np.linspace(-10, 360, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "v[-1] should be v[0] + 2pi or less" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 170, 9) * np.pi / 180.
            lons = np.linspace(10, 350, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data, s=-1)
        assert "s should be positive" in str(exc_info.value)

    def test_derivatives_grid(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])

        lut = RectSphereBivariateSpline(x,y,z)

        y = linspace(0.02, 2*pi-0.02, 7)
        x = linspace(0.02, pi-0.02, 7)

        xp_assert_close(lut(x, y, dtheta=1), _numdiff_2d(lut, x, y, dx=1),
                        rtol=1e-4, atol=1e-4)
        xp_assert_close(lut(x, y, dphi=1), _numdiff_2d(lut, x, y, dy=1),
                        rtol=1e-4, atol=1e-4)
        xp_assert_close(lut(x, y, dtheta=1, dphi=1),
                        _numdiff_2d(lut, x, y, dx=1, dy=1, eps=1e-6),
                        rtol=1e-3, atol=1e-3)

        xp_assert_equal(lut(x, y, dtheta=1),
                           lut.partial_derivative(1, 0)(x, y))
        xp_assert_equal(lut(x, y, dphi=1),
                           lut.partial_derivative(0, 1)(x, y))
        xp_assert_equal(lut(x, y, dtheta=1, dphi=1),
                           lut.partial_derivative(1, 1)(x, y))

        xp_assert_equal(lut(x, y, dtheta=1, grid=False),
                           lut.partial_derivative(1, 0)(x, y, grid=False))
        xp_assert_equal(lut(x, y, dphi=1, grid=False),
                           lut.partial_derivative(0, 1)(x, y, grid=False))
        xp_assert_equal(lut(x, y, dtheta=1, dphi=1, grid=False),
                           lut.partial_derivative(1, 1)(x, y, grid=False))

    def test_derivatives(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])

        lut = RectSphereBivariateSpline(x,y,z)

        y = linspace(0.02, 2*pi-0.02, 7)
        x = linspace(0.02, pi-0.02, 7)

        assert lut(x, y, dtheta=1, grid=False).shape == x.shape
        xp_assert_close(lut(x, y, dtheta=1, grid=False),
                        _numdiff_2d(lambda x,y: lut(x,y,grid=False), x, y, dx=1),
                        rtol=1e-4, atol=1e-4)
        xp_assert_close(lut(x, y, dphi=1, grid=False),
                        _numdiff_2d(lambda x,y: lut(x,y,grid=False), x, y, dy=1),
                        rtol=1e-4, atol=1e-4)
        xp_assert_close(lut(x, y, dtheta=1, dphi=1, grid=False),
                        _numdiff_2d(lambda x,y: lut(x,y,grid=False),
                                    x, y, dx=1, dy=1, eps=1e-6),
                        rtol=1e-3, atol=1e-3)

    def test_invalid_input_2(self):
        data = np.dot(np.atleast_2d(90. - np.linspace(-80., 80., 18)).T,
                      np.atleast_2d(180. - np.abs(np.linspace(0., 350., 9)))).T

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(0, 170, 9) * np.pi / 180.
            lons = np.linspace(0, 350, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "u should be between (0, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 180, 9) * np.pi / 180.
            lons = np.linspace(0, 350, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "u should be between (0, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 170, 9) * np.pi / 180.
            lons = np.linspace(-181, 10, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "v[0] should be between [-pi, pi)" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 170, 9) * np.pi / 180.
            lons = np.linspace(-10, 360, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data)
        assert "v[-1] should be v[0] + 2pi or less" in str(exc_info.value)

        with assert_raises(ValueError) as exc_info:
            lats = np.linspace(10, 170, 9) * np.pi / 180.
            lons = np.linspace(10, 350, 18) * np.pi / 180.
            RectSphereBivariateSpline(lats, lons, data, s=-1)
        assert "s should be positive" in str(exc_info.value)

    def test_array_like_input(self):
        y = linspace(0.01, 2 * pi - 0.01, 7)
        x = linspace(0.01, pi - 0.01, 7)
        z = array([[1, 2, 1, 2, 1, 2, 1], [1, 2, 1, 2, 1, 2, 1],
                   [1, 2, 3, 2, 1, 2, 1],
                   [1, 2, 2, 2, 1, 2, 1], [1, 2, 1, 2, 1, 2, 1],
                   [1, 2, 2, 2, 1, 2, 1],
                   [1, 2, 1, 2, 1, 2, 1]])
        # np.array input
        spl1 = RectSphereBivariateSpline(x, y, z)
        # list input
        spl2 = RectSphereBivariateSpline(x.tolist(), y.tolist(), z.tolist())
        assert_array_almost_equal(spl1(x, y), spl2(x, y))

    def test_negative_evaluation(self):
        lats = np.array([25, 30, 35, 40, 45])
        lons = np.array([-90, -85, -80, -75, 70])
        mesh = np.meshgrid(lats, lons)
        data = mesh[0] + mesh[1]  # lon + lat value
        lat_r = np.radians(lats)
        lon_r = np.radians(lons)
        interpolator = RectSphereBivariateSpline(lat_r, lon_r, data)
        query_lat = np.radians(np.array([35, 37.5]))
        query_lon = np.radians(np.array([-80, -77.5]))
        data_interp = interpolator(query_lat, query_lon)
        ans = np.array([[-45.0, -42.480862],
                        [-49.0625, -46.54315]])
        assert_array_almost_equal(data_interp, ans)

    def test_pole_continuity_gh_14591(self):
        # regression test for https://github.com/scipy/scipy/issues/14591
        # with pole_continuty=(True, True), the internal work array size
        # was too small, leading to a FITPACK data validation error.

        # The reproducer in gh-14591 was using a NetCDF4 file with
        # 361x507 arrays, so here we trivialize array sizes to a minimum
        # which still demonstrates the issue.
        u = np.arange(1, 10) * np.pi / 10
        v = np.arange(1, 10) * np.pi / 10
        r = np.zeros((9, 9))
        for p in [(True, True), (True, False), (False, False)]:
            RectSphereBivariateSpline(u, v, r, s=0, pole_continuity=p)


def _numdiff_2d(func, x, y, dx=0, dy=0, eps=1e-8):
    if dx == 0 and dy == 0:
        return func(x, y)
    elif dx == 1 and dy == 0:
        return (func(x + eps, y) - func(x - eps, y)) / (2*eps)
    elif dx == 0 and dy == 1:
        return (func(x, y + eps) - func(x, y - eps)) / (2*eps)
    elif dx == 1 and dy == 1:
        return (func(x + eps, y + eps) - func(x - eps, y + eps)
                - func(x + eps, y - eps) + func(x - eps, y - eps)) / (2*eps)**2
    else:
        raise ValueError("invalid derivative order")


class Test_DerivedBivariateSpline:
    """Test the creation, usage, and attribute access of the (private)
    _DerivedBivariateSpline class.
    """
    def setup_method(self):
        x = np.concatenate(list(zip(range(10), range(10))))
        y = np.concatenate(list(zip(range(10), range(1, 11))))
        z = np.concatenate((np.linspace(3, 1, 10), np.linspace(1, 3, 10)))
        with pytest.warns(UserWarning, match="\nThe coefficients of the spline"):
            self.lut_lsq = LSQBivariateSpline(x, y, z,
                                              linspace(0.5, 19.5, 4),
                                              linspace(1.5, 20.5, 4),
                                              eps=1e-2)
        self.lut_smooth = SmoothBivariateSpline(x, y, z)
        xx = linspace(0, 1, 20)
        yy = xx + 1.0
        zz = array([np.roll(z, i) for i in range(z.size)])
        self.lut_rect = RectBivariateSpline(xx, yy, zz, kx=4, ky=4)
        self.lut_rect_custom = _regrid(xx, yy, zz, kx=4, ky=4)
        self.orders = list(itertools.product(range(3), range(3)))

    def test_creation_from_LSQ(self):
        for nux, nuy in self.orders:
            lut_der = self.lut_lsq.partial_derivative(nux, nuy)
            a = lut_der(3.5, 3.5, grid=False)
            b = self.lut_lsq(3.5, 3.5, dx=nux, dy=nuy, grid=False)
            assert a == b

    def test_creation_from_Smooth(self):
        for nux, nuy in self.orders:
            lut_der = self.lut_smooth.partial_derivative(nux, nuy)
            a = lut_der(5.5, 5.5, grid=False)
            b = self.lut_smooth(5.5, 5.5, dx=nux, dy=nuy, grid=False)
            assert a == b

    def test_creation_from_Rect(self):
        for nux, nuy in self.orders:
            lut_der = self.lut_rect.partial_derivative(nux, nuy)
            a = lut_der(0.5, 1.5, grid=False)
            b = self.lut_rect(0.5, 1.5, dx=nux, dy=nuy, grid=False)
            assert a == b

        for nux, nuy in self.orders:
            lut_der = self.lut_rect_custom.derivative([nux, nuy])
            a = _ndbspline_call_like_bivariate(lut_der, 0.5, 1.5)
            b = _ndbspline_call_like_bivariate(
                self.lut_rect_custom, 0.5, 1.5, dx=nux, dy=nuy)
            xp_assert_close(a, b, rtol=0, atol=1.5*10**(-6))

        for nux, nuy in self.orders:
            lut_der = self.lut_rect.partial_derivative(nux, nuy)
            lut_ndspline = convert_to_ndbspline(self.lut_rect)
            lut_der_ndbspline = lut_ndspline.derivative((nux, nuy))
            a = lut_der(0.5, 1.5, grid=False)
            a_ndbspline = lut_der_ndbspline([(0.5, 1.5)])
            b = self.lut_rect(0.5, 1.5, dx=nux, dy=nuy, grid=False)
            b_ndbspline = lut_ndspline([(0.5, 1.5)], nu=(nux, nuy))
            assert a == b
            assert_almost_equal(a_ndbspline, b_ndbspline)

    def test_invalid_attribute_fp(self):
        der = self.lut_rect.partial_derivative(1, 1)
        with assert_raises(AttributeError):
            der.fp

        der = self.lut_rect_custom.derivative([1, 1])
        with assert_raises(AttributeError):
            der.fp

        lut_ndbspline = convert_to_ndbspline(self.lut_rect)
        der = lut_ndbspline.derivative(nu=[1, 1])
        with assert_raises(AttributeError):
            der.fp

    def test_invalid_attribute_get_residual(self):
        der = self.lut_smooth.partial_derivative(1, 1)
        with assert_raises(AttributeError):
            der.get_residual()
