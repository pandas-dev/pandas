# Licensed under a 3-clause BSD style license - see LICENSE.rst
import random

import pytest

from asv import step_detect
from asv.step_detect import (detect_regressions, detect_steps, golden_search, median,
                             rolling_median_dev, solve_potts, solve_potts_approx,
                             solve_potts_autogamma)

try:
    import numpy as np
    np.random.seed
    HAVE_NUMPY = True
except (ImportError, NameError):
    HAVE_NUMPY = False


@pytest.mark.skipif(not HAVE_NUMPY, reason="test needs numpy")
def test_solve_potts(use_rangemedian):
    np.random.seed(1234)

    # Easy case, exact solver
    y = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    right, values, dists = solve_potts(y, w=[1] * len(y), gamma=0.1)
    assert right == [3, 6, 9]
    assert np.allclose(values, [1, 2, 3], atol=0)

    right, values, dists = solve_potts(y, w=[1] * len(y), gamma=8.0)
    assert right == [9]
    assert np.allclose(values, [2], atol=0)

    # l1 norm
    right, values, dists = solve_potts(y, w=[1] * len(y), gamma=0.1)
    assert right == [3, 6, 9]
    assert np.allclose(values, [1, 2, 3], atol=0)

    # Bigger case, exact solver
    t = np.arange(100)
    y0 = (+ 0.4 * (t >= 5) +
          0.9 * (t >= 10) -
          0.2 * (t >= 20) +
          0.2 * (t >= 50) +
          1.1 * (t >= 70))
    y = y0 + 0.05 * np.random.rand(y0.size)
    right, values, dists = solve_potts(y.tolist(), w=[1] * len(y), gamma=0.1)
    assert right == [5, 10, 20, 50, 70, 100]

    # Bigger case, approximative solver
    right, values, dists = solve_potts_approx(y.tolist(), w=[1] * len(y), gamma=0.1)
    assert right == [5, 10, 20, 50, 70, 100]

    # Larger case
    t = np.arange(3000)
    y0 = (+ 0.4 * (t >= 10) +
          0.9 * (t >= 30) -
          0.2 * (t >= 200) +
          0.2 * (t >= 600) +
          1.1 * (t >= 2500) -
          0.5 * (t >= 2990))

    # Small amount of noise shouldn't disturb step finding
    y = y0 + 0.05 * np.random.randn(y0.size)
    right, values, dists, gamma = solve_potts_autogamma(y.tolist(), w=[1] * len(y))
    assert right == [10, 30, 200, 600, 2500, 2990, 3000]

    # Large noise should prevent finding any steps
    y = y0 + 5.0 * np.random.randn(y0.size)
    right, values, dists, gamma = solve_potts_autogamma(y.tolist(), w=[1] * len(y))
    assert right == [3000]

    # The routine shouldn't choke on datasets with 10k points.
    # Appending noisy data to weakly noisy data should retain the
    # steps in the former
    y = y0 + 0.025 * np.random.rand(y.size)
    ypad = 0.025 * np.random.randn(10000 - 3000)
    right, values, dists, gamma = solve_potts_autogamma(y.tolist() + ypad.tolist(),
                                                        w=[1] * (len(y) + len(ypad)))
    assert right == [10, 30, 200, 600, 2500, 2990, 3000, 10000]


@pytest.mark.skipif(not HAVE_NUMPY, reason="test needs numpy")
def test_autocorrelated():
    # Check that a low-amplitude cosine signal is not interpreted as
    # multiple steps
    j = np.arange(1000)
    y = 0.2 * np.cos(j / 100.0) + 1.0 * (j >= 500)
    right, values, dists, gamma = solve_potts_autogamma(y.tolist(), w=[1] * len(y))
    assert right == [500, 1000]


def test_zero_variance():
    # Should not choke on this data
    y = [1.0] * 1000
    right, values, dists, gamma = solve_potts_autogamma(y, w=[1] * len(y))
    assert right == [1000]


@pytest.mark.skipif(not HAVE_NUMPY, reason="test needs numpy")
def test_weighted():
    np.random.seed(1234)

    t = np.arange(100)
    y0 = (+ 0.4 * (t >= 5) +
          0.9 * (t >= 10) -
          0.2 * (t >= 20) +
          0.2 * (t >= 50) +
          1.1 * (t >= 70))
    y = y0 + 0.05 * np.random.rand(y0.size)

    y = y.tolist()
    w = [1] * len(y)

    y[15] = 2
    right, values, dists = solve_potts(y, w=w, gamma=0.1)
    assert right == [5, 10, 15, 16, 20, 50, 70, 100]

    steps = detect_steps(y, w=w)
    steps_pos = [s[0] for s in steps]
    assert steps_pos == [0, 5, 10, 15, 16, 20, 50, 70]

    w[15] = 0.1
    right, values, dists = solve_potts(y, w=w, gamma=0.1)
    assert right == [5, 10, 20, 50, 70, 100]

    steps = detect_steps(y, w=w)
    steps_pos = [s[0] for s in steps]
    assert steps_pos == [0, 5, 10, 20, 50, 70]

    # Missing weights and data should be handled properly
    y[35] = None
    w[23] = None

    steps = detect_steps(y, w=w)
    steps_pos = [s[0] for s in steps]
    assert steps_pos == [0, 5, 10, 20, 50, 70]


@pytest.mark.skipif(not HAVE_NUMPY, reason="test needs numpy")
def test_detect_regressions(use_rangemedian):
    np.random.seed(1234)

    for seed in [1234, 5678, 8901, 2345]:
        np.random.seed(seed)
        t = np.arange(4000)
        y = 0.7 * np.random.rand(t.size)

        y -= 0.3 * (t >= 1000)
        y += 2.0 * (t >= 3234 + (seed % 123))
        y += 2.0 * ((t >= 3264 + (seed % 53)) & (t < 3700))
        y -= 2.7 * ((t >= 3350) & (t < 3500 + (seed % 71)))

        y = y.tolist()
        y[123] = None
        y[1234] = np.nan
        steps = detect_steps(y)

        steps_lr = [(l, r) for l, r, _, _, _ in steps]
        k = steps[0][1]
        assert 990 <= k <= 1010
        assert steps_lr == [(0, k),
                            (k, 3234 + (seed % 123)),
                            (3234 + (seed % 123), 3264 + (seed % 53)),
                            (3264 + (seed % 53), 3350),
                            (3350, 3500 + (seed % 71)),
                            (3500 + (seed % 71), 3700),
                            (3700, 4000)]
        steps_v = [x[2] for x in steps]
        assert np.allclose(steps_v, [0.35, 0.05, 2.05, 4.05, 1.15, 4.05, 2.05], rtol=0.3)

        # The expected mean error is 0.7 <|U(0,1) - 1/2|> = 0.7/4
        steps_err = [x[4] for x in steps]
        assert np.allclose(steps_err, [0.7 / 4] * 7, rtol=0.3)

        # Check detect_regressions
        new_value, best_value, regression_pos = detect_regressions(steps)
        assert regression_pos == [(3233 + (seed % 123), (3233 + (seed % 123) + 1),
                                   steps_v[1], min(steps_v[2:])),
                                  (3499 + (seed % 71), 3499 + (seed % 71) + 1,
                                   steps_v[4], min(steps_v[5:]))]
        assert np.allclose(best_value, 0.7 / 2 - 0.3, rtol=0.3, atol=0)
        assert np.allclose(new_value, 0.7 / 2 - 0.3 + 2, rtol=0.3, atol=0)


def test_golden_search():
    def f(x):
        return 1 + x**3 + x**4
    x = golden_search(f, -1, -0.25, xatol=1e-5, ftol=0)
    assert abs(x - (-3 / 4)) < 1e-4
    x = golden_search(f, -0.25, 0.25, xatol=1e-5, ftol=0)
    assert abs(x - (-0.25)) < 1e-4


def rolling_median_dev_naive(items):
    for j in range(1, len(items)):
        m = median(items[:j])
        d = sum(abs(x - m) for x in items[:j])
        yield m, d


def test_rolling_median():
    random.seed(1)

    datasets = [
        [1, 1, 10, 3, 5, 1, -16, -3, 4, 9],
        [random.gauss(0, 1) for j in range(500)]
    ]

    for x in datasets:
        got = list(rolling_median_dev(x))
        expected = rolling_median_dev_naive(x)
        for j, b in enumerate(expected):
            a = got[j]
            assert abs(a[0] - b[0]) < 1e-10, (a, b)
            assert abs(a[1] - b[1]) < 1e-10, (a, b)


def test_l1dist(use_rangemedian):
    random.seed(1)

    datasets = [
        ([1, 1, 10, 3, 5, 1, -16, -3, 4, 9], [1] * 10),
        ([random.gauss(0, 1) for j in range(50)], [1] * 50),
        ([1, 2, 3, 4], [1, 2, 1, 2])
    ]

    def median_iter(y, w):
        if all(ww == w[0] for ww in w):
            for m, d in rolling_median_dev(y):
                yield m, d
        else:
            for j in range(1, len(y) + 1):
                m = step_detect.weighted_median(y[:j], w[:j])
                d = sum(ww * abs(yy - m) for yy, ww in zip(y[:j], w[:j]))
                yield m, d

    for y, w in datasets:
        dist = step_detect.get_mu_dist(y, w)

        for i in range(len(y)):
            for p, (m2, d2) in enumerate(median_iter(y[i:], w[i:])):
                j = i + p

                m = dist.mu(i, j)
                d = dist.dist(i, j)

                assert m == m2, (i, j)
                assert abs(d - d2) < 1e-10, (i, j)


def test_regression_threshold():
    steps = [(0, 1, 1.0, 1.0, 0.0),
             (1, 2, 1.1, 1.1, 0.0),
             (2, 3, 2.0, 2.0, 0.0)]

    latest, best, pos = detect_regressions(steps, threshold=0.05, min_size=1)
    assert latest == 2
    assert best == 1
    assert pos == [(0, 1, 1.0, 1.1), (1, 2, 1.1, 2.0)]

    latest, best, pos = detect_regressions(steps, threshold=0.2, min_size=1)
    assert latest == 2
    assert best == 1
    assert pos == [(1, 2, 1.1, 2.0)]

    latest, best, pos = detect_regressions(steps, threshold=0.8, min_size=1)
    assert latest == 2
    assert best == 1
    assert pos == [(1, 2, 1.1, 2.0)]

    latest, best, pos = detect_regressions(steps, threshold=1.1, min_size=1)
    assert latest is None
    assert best is None
    assert pos is None

    steps = [(0, 1, 1.0, 1.0, 0.0),
             (1, 2, 1.3, 1.3, 0.0),
             (2, 3, 1.1, 1.1, 0.0)]

    latest, best, pos = detect_regressions(steps, threshold=0.2, min_size=1)
    assert latest is None
    assert best is None
    assert pos is None

    # Gradual change should result to a regression detected somewhere,
    # even if the individual steps are smaller than the threshold
    steps = [(0, 1, 1.0, 1.0, 0.0),
             (1, 2, 1.04, 1.04, 0.0),
             (2, 3, 1.08, 1.08, 0.0)]

    latest, best, pos = detect_regressions(steps, threshold=0.05)
    assert pos == [(0, 1, 1.0, 1.04)]


def test_zero_weight():
    t = list(range(50))
    y = [1.0 * (tt > 25) for tt in t]
    w = [0 if tt >= 20 and tt < 30 else 1 for tt in t]

    # Due to zero weight, part of the data is not assigned to belong to a step
    steps = detect_steps(y, w=w)
    steps_pos = [s[:2] for s in steps]
    assert steps_pos == [(0, 20), (30, 50)]

    # A small weight allows to find the actual step
    w = [1e-6 if ww == 0 else 1 for ww in w]
    steps = detect_steps(y, w=w)
    steps_pos = [s[:2] for s in steps]
    assert steps_pos == [(0, 26), (26, 50)]


def test_regression_min_size():
    steps = [(0, 2, 1.0, 1.0, 0.0),
             (2, 3, 0.0, 0.0, 0.0),
             (3, 5, 1.0, 1.0, 0.0)]

    latest, best, pos = detect_regressions(steps)
    assert latest is None and best is None and pos is None

    latest, best, pos = detect_regressions(steps, min_size=1)
    assert latest == 1.0
    assert best == 0.0
    assert pos == [(2, 3, 0.0, 1.0)]

    steps = [(0, 2, 1.0, 1.0, 0.0),
             (2, 3, 0.0, 0.0, 0.0),
             (3, 5, 2.0, 1.0, 0.0)]

    latest, best, pos = detect_regressions(steps)
    assert latest == 2.0
    assert best == 0.0
    assert pos == [(2, 3, 0.0, 2.0)]


def test_solve_potts_approx_bug():
    y = [2.9, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1,
         3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1,
         3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.2, 3.1, 3.1, 3.1, 3.1,
         3.1, 3.1, 3.1, 3.1, 2.9, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1, 3.1]
    w = [86.1, 1.0, 1.0, 1.0, 0.8, 1.0, 0.8, 1.0, 1.0, 0.9, 0.9, 1.0,
         0.8, 1.0, 1.0, 1.0, 1.0, 0.9, 0.6, 0.9, 0.5, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 1.0, 1.0, 0.7, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 0.8, 0.8, 50.0, 0.8, 1.0, 1.0, 0.6, 0.8, 1.0, 1.0]
    gamma = 0.3

    r0, v0, d0 = solve_potts(y, w, gamma=gamma)
    r, v, d = solve_potts_approx(y, w, gamma=gamma)
    assert r == r0
