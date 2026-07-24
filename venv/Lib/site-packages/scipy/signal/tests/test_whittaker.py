import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.signal import whittaker_henderson  # type: ignore[attr-defined]
from scipy.signal._whittaker import (
    _logdet_difference_matrix, _polynomial_fit, _reml, _solveh_banded, _solve_WH_banded,
)
from scipy.stats import special_ortho_group


def _solve_WH_order2_fast(y, lamb):
    """Efficiently solve WH of order 2 according to Weinert.

    Needs order = 2 and n >= 3 (data points).

    Weinert (2007)
    "Efficient computation for Whittaker-Henderson smoothing".
    Computational Statistics and Data Analysis 52:959-74.
    https://doi.org/10.1016/j.csda.2006.11.038
    """
    n = y.shape[0]
    # Convert penalty to convention of Weinert (2007), i.e. A = lambda I + D'D.
    # Note that Weinert denotes the difference matrix M instead of D.
    lamb = 1 / lamb
    # A = LDL' decompositon, i.e. L is unit lower triangular and banded
    # (bandwith=2) and D diagonal.
    # First subdiagonal of L is (-e_1, .., -e_{n-1}) and 2nd subdiagonal
    # (f_1, .., f_{n-2}). Diagonal of D is (d_1, .., d_n).
    # The equation A @ x = lamb * y becomes
    # LD @ b = lamb * y  (I)
    # L.T @ x = b        (II)
    # We shift Weinert's 1-based indices to 0-based indices, Eq. 2.2-2.6
    # Solve problem (I)
    b = np.empty(n)
    e = np.empty(n)
    f = np.empty(n)
    # i=0
    d = 1 + lamb
    f[0] = 1 / d
    mu = 2
    e[0] = mu * f[0]
    b[0] = f[0] * lamb * y[0]
    mu_old = mu
    # i=1
    if n == 3:
        d = 4 + lamb - mu_old * e[0]
        mu = 2 - e[0]
    else:
        d = 5 + lamb - mu_old * e[0]
        mu = 4 - e[0]
    f[1] = 1 / d
    e[1] = mu * f[1]
    b[1] = f[1] * (lamb * y[1] + mu_old * b[0])
    mu_old = mu
    for i in range(2, n-2):
        d = 6 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        mu = 4 - e[i-1]
        e[i] = mu * f[i]
        b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])
        mu_old = mu
    # i=n-2
    if n >= 4:
        i = n - 2
        d = 5 + lamb - mu_old * e[i-1] - f[i-2]
        f[i] = 1 / d
        mu = 2 - e[i-1]
        e[i] = mu * f[i]
        b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])
        mu_old = mu
    # i=n-1
    i = n - 1
    d = 1 + lamb - mu_old * e[i-1] - f[i-2]
    f[i] = 1 / d
    b[i] = f[i] * (lamb * y[i] + mu_old * b[i-1] - b[i-2])

    # Solve problem (II)
    x = np.empty(n)
    x[n-1] = b[n-1]
    x[n-2] = b[n-2] + e[n-2] * x[n-1]
    for i in range(n-3, -1, -1):
        x[i] = b[i] + e[i] * x[i+1] - f[i] * x[i+2]
    return x


def test_solveh_banded():
    n = 10
    # construct positive definite matrix A
    Q = special_ortho_group(n, seed=1234).rvs()
    A = Q @ np.diag(1 + np.arange(n)) @ Q.T
    b = np.arange(n)
    ab = np.zeros((n, n))
    ab[0, :] = np.diagonal(A, 0)
    for i in range(1, n):
        ab[i, :-i] = np.diagonal(A, i)
    x, logdet, _ = _solveh_banded(ab, b, calc_logdet=True)
    assert_allclose(x, np.linalg.solve(A, b), rtol=1e-11)
    assert_allclose(logdet, np.log(np.linalg.det(A)), rtol=1e-11)

    # tridiagonal case
    x, logdet, _ = _solveh_banded(ab[:2], b, calc_logdet=True)
    A_tri = (
        np.diag(np.diag(A)) + np.diag(np.diag(A, -1), -1) + np.diag(np.diag(A, 1), 1)
    )
    assert_allclose(logdet, np.log(np.linalg.det(A_tri)), rtol=1e-11)


@pytest.mark.parametrize(
        ["order", "n"],
        [
            (1, 2), (1, 10),
            (2, 3), (2, 4), (2, 5), (2, 20),
            (5, 6), (5, 7), (5, 10), (5, 20),
            (6, 7), (6, 10), (6, 20),
        ])
def test_logdet_difference_matrix(order, n):
    p = order
    logdet = _logdet_difference_matrix(order=p, n=n)
    D = np.diff(np.eye(n), n=p, axis=0)  # shape (n-p, n)
    assert_allclose(logdet, np.log(np.linalg.det(D @ D.T)), rtol=1e-11)
    eigvals = np.linalg.eigvals(D.T @ D)
    assert_allclose(logdet, np.sum(np.log(eigvals[eigvals > 1e-8])), rtol=1e-11)


@pytest.mark.parametrize("calc_logdet", [False, True])
@pytest.mark.parametrize("lamb", [0.1, 1e20])
def test_polynomial_fit(calc_logdet, lamb):
    """Test that _polynomial_fit works as expected."""
    n = 20
    y = np.sin(2 * np.pi * np.linspace(0, 1, n))
    poly, logdet = _polynomial_fit(y, lamb, calc_logdet=calc_logdet)

    if not calc_logdet:
        assert logdet == 0
    
    poly2, logdet2 = _polynomial_fit(
        y, lamb, weights=1.23 * np.ones(n), calc_logdet=calc_logdet
    )
    assert logdet2 == pytest.approx(logdet, rel=1e-12)
    assert_allclose(poly2, poly)


@pytest.mark.parametrize(
        ["signal", "lamb", "order", "weights", "err", "msg"],
        [
            ([[1, 2, 3] * 3], 1, 2, None, ValueError,
             "Input array signal must be of shape \\(n,\\)"),
            (np.zeros(2), 1, 2, None, ValueError,
             "Input array signal must be at least of shape"),
            (np.arange(10), -0.9, 2, None, ValueError,
             "Parameter lamb must be string"),
            ([1, 2, 3], "XXX", 1, None, ValueError,
             "Parameter lamb must be string 'reml'"),
            ([1, 2, 3], 1, 1.5, None, TypeError,
             "order must be an integer."),
            ([1, 2, 3], 1, 0, None, ValueError,
             "order must be an integer not less than 1"),
            ([1, 2, 3], 1, 2, [0, 1], ValueError,
             "Input array weights must have the same shape as the signal array."),
            ([1, 2, 3], 1, 2, [[0]], ValueError,
             "Input array weights must have the same shape as the signal array."),
            ([1, 2, np.nan], 1, 2, None, ValueError,
             "Input array signal must be finite."),
            ([1, 2, 3], 1, 2, [1, 2, np.nan], ValueError,
             "Input array weights must be finite."),
            ([1, 2, np.nan], 1, 2, [1, 2, 3], ValueError,
             "Input array weights must be zero for all non-finite"),
        ])
def test_whittaker_raises(signal, lamb, order, weights, err, msg):
    """Test that whittaker raises errors."""
    with pytest.raises(err, match=msg):
        whittaker_henderson(signal, lamb=lamb, order=order, weights=weights)


def test_whittaker_small_data():
    """Test that whittaker works on a few data points."""
    # Should work on order + 1 data points. The first 2*order+1 are special.
    la = 1
    y = np.arange(2)
    wh = whittaker_henderson(y, order=1, lamb=la)
    # Analytical solution for order=1 and n=2
    res = (np.array([[1 + la, la], [la, 1 + la]]) / (1 + 2 * la)) @ y
    assert_allclose(wh.x, res, atol=1e-15)
    whittaker_henderson(np.arange(3), order=1, lamb=la)

    y = np.arange(3)
    wh = whittaker_henderson(y, order=2, lamb=la)
    # Analytical solution for order=2 and n=3:
    # import numpy as np
    # from sympy import eye, Symbol
    # n, order = 3, 2
    # la = Symbol("la")
    # D = np.diff(eye(n), n=order, axis=0)
    # A = eye(n) + la * D.T @ D
    # A.inv()
    res = (
        np.array([
            [1 + 5 * la,     2 * la,        -la],
            [    2 * la, 1 + 2 * la,     2 * la],
            [       -la,     2 * la, 1 + 5 * la],
        ]) / (1 + 6 * la)) @ y
    assert_allclose(wh.x, res, atol=1e-15)
    assert wh.lamb == la
    
    y = np.arange(4)
    wh = whittaker_henderson(y, order=2, lamb=la)
    la2 = la ** 2
    res =(np.array([
        [14*la2 + 11*la + 1, 8*la2 + 2*la    , 2*la2 -   la,     -4*la2            ],
        [ 8*la2 +  2*la    , 6*la2 + 7*la + 1, 4*la2 + 4*la,      2*la2 -  la      ],
        [ 2*la2 -    la    , 4*la2 + 4*la    , 6*la2 + 7*la + 1,  8*la2 + 2*la     ],
        [-4*la2            , 2*la2 -   la    , 8*la2 + 2*la    , 14*la2 + 11*la + 1],
    ]) /(20*la2 + 12*la + 1)) @ y
    assert_allclose(wh.x, res, atol=1e-15)
    whittaker_henderson(np.arange(5), order=2, lamb=la)

    whittaker_henderson(np.arange(4), order=3, lamb=la)
    whittaker_henderson(np.arange(5), order=3, lamb=la)
    whittaker_henderson(np.arange(6), order=3, lamb=la)
    whittaker_henderson(np.arange(7), order=3, lamb=la)


@pytest.mark.parametrize("n", [3, 4, 5, 100])
def test_whittaker_direct_vs_fast_order2(n):
    """Test equivalent results"""
    rng = np.random.default_rng(42)
    signal = np.sin(2 * np.pi * np.linspace(0, 1, n)) + rng.standard_normal(n)
    x1, _ = _solve_WH_banded(signal, lamb=1.23, order=2)
    x2 = _solve_WH_order2_fast(signal, lamb=1.23)
    assert_allclose(x1, x2, rtol=1e-11)

    wh = whittaker_henderson(signal, lamb=1.23, order=2)
    assert_allclose(wh.x, x1, rtol=1e-11)


@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_whittaker_limit_zero_penalty(order):
    """Test whittaker for penalty lamb close to zero."""
    rng = np.random.default_rng(42)
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    noise = rng.standard_normal(n)
    signal = y + noise

    # In the limit of a zero penalty, the smoothing results in reproducing the signal.
    wh = whittaker_henderson(signal, lamb=1e-7, order=order)
    assert_allclose(wh.x, signal, rtol=1e-4, atol=1e-5)
    assert wh.lamb == 1e-7


@pytest.mark.parametrize("weights", [None, "ones"])
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_whittaker_limit_infinity_penalty(weights, order):
    """Test whittaker for penalty lamb close to infinity."""
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    w = None if weights is None else np.ones_like(y)

    # In the limit of an infinite penalty, the smoothing results in an interpolation
    # polynom of degree = penalty order - 1 (order=2 => linear)
    wh = whittaker_henderson(y, lamb=1e11, order=order, weights=w)
    x_poly = np.arange(len(y))
    poly = np.polynomial.Polynomial.fit(x=x_poly, y=y, deg=order - 1)
    assert_allclose(wh.x, poly(x_poly), rtol=10**(-6 + order), atol=1e-5)
    if order == 2:
        # Linear interpolation:
        # As the sine is positive fom 0 to pi and negative from pi to 2*pi, we expect
        # a negative slope.
        assert np.diff(wh.x)[0] < 0

    if order > 2:
        # Very large lambda should trigger a user warning and produce the exact
        # polynomial fit (separate code path where linear solver breaks down).
        with pytest.warns(
            UserWarning,
            match="The linear solver in Whittaker-Henderson smoothing detected",
        ):
            wh = whittaker_henderson(y, lamb=1e15, order=order, weights=w)
    else:
        wh = whittaker_henderson(y, lamb=1e20, order=order, weights=w)
    assert_allclose(wh.x, poly(x_poly), rtol=1e-10)


def test_whittaker_unpenalized():
    """Test whittaker for lamb=0."""
    n = 10
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    wh = whittaker_henderson(y, lamb=0)
    assert_allclose(wh.x, y, rtol=1e-11)
    assert not np.may_share_memory(wh.x, y)


def test_whittaker_weights():
    """Test that whittaker with weights of 1 is same as without weights."""
    rng = np.random.default_rng(42)
    n = 100
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    noise = rng.standard_normal(n)
    signal = y + noise
    wh1 = whittaker_henderson(signal, lamb=1)
    w = np.ones_like(signal)
    wh2 = whittaker_henderson(signal, lamb=1, weights=w)
    assert_allclose(wh1.x, wh2.x, rtol=1e-11)

    # Multiplying penalty and weights by the same number does not change the result.
    wh3 = whittaker_henderson(signal, lamb=3, weights=3*w)
    assert_allclose(wh3.x, wh1.x, rtol=1e-11)


@pytest.mark.parametrize("missing_value", [np.nan, np.inf, 1e-20])
def test_whittaker_zero_weight_interpolation(missing_value):
    """Test that whittaker interpolates where weights are zero."""
    n = 100
    signal = np.sin(2*np.pi * np.linspace(0, 1, n))
    signal[50:] += 2
    w = np.ones_like(signal)
    signal[40:60] = missing_value  # value does not matter
    w[40:60] = 0
    # Note: interpolation is a polynomial of degree = 2 * order - 1.

    # order = 1 => linear interpolation
    wh = whittaker_henderson(signal, lamb=1, order=1, weights=w)
    interp = np.interp(np.arange(40, 60), [39, 60], [wh.x[39], wh.x[60]])
    assert_allclose(wh.x[40:60], interp, rtol=1e-11)

    # order = 2 => cubic interpolation
    wh = whittaker_henderson(signal, lamb=1, order=2, weights=w)
    poly = np.polynomial.Polynomial.fit(
        x=[38, 39, 60, 61], y=[wh.x[38], wh.x[39], wh.x[60], wh.x[61]], deg=3
    )
    assert_allclose(wh.x[40:60], poly(np.arange(40, 60)), rtol=1e-11)


def test_whittaker_zero_weight_extrapolation():
    """Test that whittaker extrapolates where weights are zero."""
    n = 100
    signal = np.sin(2*np.pi * np.linspace(0, 1, n))
    signal[50:] += 2
    w = np.ones_like(signal)
    signal[80:] = np.nan  # value does not matter, np.nan to check arg validation
    w[80:] = 0
    # Note: extrapolation is a polynomial of degree = order - 1.

    # order = 1 => constant extrapolation
    wh = whittaker_henderson(signal, lamb=1, order=1, weights=w)
    assert_allclose(wh.x[80:], wh.x[79], rtol=1e-11)

    # order 2 => linear extrapolation
    wh = whittaker_henderson(signal, lamb=1, order=2, weights=w)
    poly = np.polynomial.Polynomial.fit(x=[78, 79], y=wh.x[78:80], deg=1)
    assert_allclose(wh.x[80:], poly(np.arange(80, 100)), rtol=1e-11)

    # order = 3 => quadratic extrapolation
    wh = whittaker_henderson(signal, lamb=1, order=3, weights=w)
    poly = np.polynomial.Polynomial.fit(x=[77, 78, 79], y=wh.x[77:80], deg=2)
    assert_allclose(wh.x[80:], poly(np.arange(80, 100)), rtol=1e-9)


@pytest.mark.parametrize("order", [1, 2, 3])
def test_reml_criterion(order):
    la = 1.2345
    n = 10
    y = np.sin(2*np.pi * np.linspace(0, 1, n))
    D = np.diff(np.eye(n), n=order, axis=0)
    M = D.T @ D

    def test_reml(la, y):
        # See docstring and comment of _reml.
        A = np.eye(n) + la * M
        x = np.linalg.solve(A, y)
        resid = y - x
        r2 = resid @ resid + la * x @ M @ x
        return -0.5 * ((n-order) * (1 + np.log(r2 / (n - order)))
                       - np.log(np.linalg.det(la * D @ D.T))
                       + np.log(np.linalg.det(A)))

    r1 = _reml(lamb=la, y=y, order=order)
    r2 = test_reml(la=la, y=y)
    assert_allclose(r1, r2, rtol=1e-11)


@pytest.mark.parametrize("weights", [None, "ones"])
def test_whittaker_reml(weights):
    """Test that whittaker with REML criterion works."""
    n = 10
    order = 2
    signal = np.sin(2*np.pi * np.linspace(0, 1, n))
    y = signal.copy()
    # add some (deterministic) noise
    y[::2] += 0.2
    y[1::2] -= 0.2
    w = None if weights is None else np.ones_like(y)

    wh = whittaker_henderson(y, lamb="reml", order=order, weights=w)
    assert wh.lamb > 0

    # Validate against R mgcv
    # library(mgcv)
    # n = 10
    # order = 2
    # signal = sin(2*pi * 0:(n-1) / (n-1))
    # y = signal
    # # add some (deterministic) noise
    # y[seq(1, n, 2)] = y[seq(1, n, 2)] + 0.2
    # y[seq(2, n, 2)] = y[seq(2, n, 2)] - 0.2
    # x <- 1:n
    # b <- gam(y ~ s(x, bs="ps", k=n, m=c(0, 2)), family = gaussian(),
    #          knots = list(x=0:11), method="REML",
    #          control = list(epsilon = 1e-12, mgcv.tol = 1e-12, efs.tol = 1e-12))
    # options(digits=12)
    # b$fitted.values
    expected = [
        0.209959334732, 0.589996568382, 0.938965773497, 0.766582291904, 0.349462653088,
        -0.349462653088, -0.766582291904, -0.938965773497, -0.589996568382,
        -0.209959334732,
    ]
    assert_allclose(wh.x, expected, rtol=1e-7)
    # b$sp
    expected_sp = 5.12904626508
    expected_sp /= 16  # not 100% sure how to justify this 16.
    assert wh.lamb == pytest.approx(expected_sp, rel=1e-7)

    # Double weights should result in double lamb.
    wh2 = whittaker_henderson(y, lamb="reml", order=order, weights=2 * np.ones_like(y))
    assert wh2.lamb == pytest.approx(2 * wh.lamb, rel=1e-7)
