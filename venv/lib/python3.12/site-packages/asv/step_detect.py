# Licensed under a 3-clause BSD style license - see LICENSE.rst

import collections
import heapq
import math
from statistics import median

try:
    from . import _rangemedian
except ImportError:
    _rangemedian = None


#
# Detecting regressions
#


def detect_steps(y, w=None):
    """
    Detect steps in a (noisy) signal.

    Parameters
    ----------
    y : list of float, none or nan
        Single benchmark result series, with possible missing data
    w : list of float, none or nan
        Data point relative weights. Missing weights are set equal
        to the median weight.

    Returns
    -------
    steps : list of (left_pos, right_pos, value, min_value, err_est)
        List containing a decomposition of the input data to a piecewise
        constant function. Each element contains the left (inclusive) and
        right (exclusive) bounds of a segment, the average value on
        the segment, the minimum value in the segment, and the l1 error
        estimate, :math:`|Y - avg|`. Missing data points are not necessarily
        contained in any segment; right_pos-1 is the last non-missing data
        point.

    """

    index_map = {}
    y_filtered = []
    for j, x in enumerate(y):
        if x is None or x != x:
            # None or NaN: missing data
            continue
        if w is not None and w[j] is not None and (w[j] <= 0):
            # non-positive weight: consider as missing data
            continue
        index_map[len(y_filtered)] = j
        y_filtered.append(x)

    # Weights
    if w is None:
        w_filtered = [1] * len(y_filtered)
    else:
        # Fill-in and normalize weights
        w_valid = [ww for ww in w if ww is not None and ww == ww]
        if w_valid:
            w_median = median(w_valid)
            if w_median == 0:
                w_median = 1.0
        else:
            w_median = 1.0

        w_filtered = [1.0] * len(y_filtered)
        for j in range(len(w_filtered)):
            jj = index_map[j]
            if w[jj] is not None and w[jj] == w[jj]:
                w_filtered[j] = w[jj] / w_median

    # Find piecewise segments
    right, values, dists, gamma = solve_potts_autogamma(y_filtered, w=w_filtered)

    # Extract the steps, mapping indices back etc.
    steps = []
    l = 0
    for r, v, d in zip(right, values, dists):
        steps.append(
            (index_map[l], index_map[r - 1] + 1, v, min(y_filtered[l:r]), abs(d / (r - l)))
        )
        l = r
    return steps


def detect_regressions(steps, threshold=0, min_size=2):
    """Detect regressions in a (noisy) signal.

    A regression means an upward step in the signal.  The value
    'before' a regression is the value immediately preceding the
    upward step.  The value 'after' a regression is the minimum of
    values after the upward step.

    Parameters
    ----------
    steps : list of (left, right, value, min, error)
        List of steps computed by detect_steps, or equivalent
    threshold : float
        Relative threshold for reporting regressions. Filter out jumps
        whose relative size is smaller than threshold, if they are not
        necessary to explain the difference between the best and the latest
        values.
    min_size : int
        Minimum number of commits in a regression to consider it.

    Returns
    -------
    latest_value
        Latest value
    best_value
        Best value
    regression_pos : list of (before, after, value_before, best_value_after)
        List of positions between which the value increased. The first item
        corresponds to the last position at which the best value was obtained.
        The last item indicates the best value found after the regression
        (which is not always the value immediately following the regression).

    """
    if not steps:
        # No data: no regressions
        return None, None, None

    regression_pos = []

    last_v = steps[-1][2]
    best_v = last_v
    thresholded_best_v = last_v
    thresholded_best_err = steps[-1][4]
    prev_l = None
    short_prev = None

    # Find upward steps that resulted to worsened value afterward
    for l, r, cur_v, cur_min, cur_err in reversed(steps):
        threshold_step = max(cur_err, thresholded_best_err, threshold * cur_v)

        if thresholded_best_v > cur_v + threshold_step:
            if r - l < min_size:
                # Accept short intervals conditionally
                short_prev = (thresholded_best_v, thresholded_best_err)

            regression_pos.append((r - 1, prev_l, cur_v, best_v))

            thresholded_best_v = cur_v
            thresholded_best_err = cur_err
        elif short_prev is not None:
            # Ignore the previous short interval, if the level
            # is now back to where it was
            if short_prev[0] <= cur_v + threshold_step:
                regression_pos.pop()
                thresholded_best_v, thresholded_best_err = short_prev
            short_prev = None

        prev_l = l

        if cur_v < best_v:
            best_v = cur_v

    regression_pos.reverse()

    # Return results
    if regression_pos:
        return (last_v, best_v, regression_pos)
    else:
        return (None, None, None)


#
# Fitting piecewise constant functions to noisy data
#


def solve_potts(y, w, gamma, min_size=1, max_size=None, min_pos=None, max_pos=None, mu_dist=None):
    """Fit penalized stepwise constant function (Potts model) to data.

    Given a time series y = {y_1, ..., y_n}, fit series x = {x_1, ..., x_n}
    by minimizing the cost functional::

        F[x] = gamma * J(x) + sum(|y - x|**p)

    where J(x) is the number of jumps (x_{j+1} != x_j) in x.

    The algorithm used is described in
    :cite:t:`ptts-friedrichComplexityPenalizedMEstimation2008`, it uses dynamic
    programming to find an exact solution to the problem (within the constraints
    specified).

    The computational cost is ~ O(n**2 log n).

    Parameters
    ----------
    y : list of floats
        Input data series
    gamma : float
        Penalty parameter.
    min_size : int, optional
        Minimum interval size to consider
    max_size : int, optional
        Maximum interval size to consider
    mu_dist : Dist, optional
        Precomputed interval means/medians and cost function values
    min_pos : int, optional
        Start point (inclusive) for the interval grid
    max_pos : int, optional
        End point (exclusive) for the interval grid

    Returns
    -------
    right : list
        List of (exclusive) right bounds of the intervals
    values : list
        List of values of the intervals
    dist : list
        List of ``sum(|y - x|**p)`` for each interval.
    mu_dist : Dist
        Precomputed interval means/medians and cost function values

    References
    ----------
    .. bibliography::
       :filter: docname in docnames
       :labelprefix: PTTS_
       :keyprefix: ptts-

    """

    if len(y) == 0:
        return [], [], []

    if min_pos is None:
        min_pos = 0

    if len(y) != len(w):
        raise ValueError("y and w must have same size")

    if max_pos is None:
        max_pos = len(y)

    if mu_dist is None:
        mu_dist = get_mu_dist(y, w)

    if max_size is None:
        max_size = len(y)

    mu, dist = mu_dist.mu, mu_dist.dist

    if min_size >= max_pos - min_pos:
        return [len(y)], [mu(0, len(y) - 1)], [dist(0, len(y) - 1)]

    # Perform the Bellman recursion for the optimal partition.
    # Routine "Find best partition" in [1]
    #
    # Computes:
    #
    # p : list, length n
    #     Set of intervals, represented as follows:
    #     For interval (inclusive) right edge r in {0, ..., n-1},
    #     the best (exclusive) left edge is at l=p[r].
    #     Where intervals overlap, the rightmost one has priority.

    if hasattr(mu_dist, 'find_best_partition'):
        p = mu_dist.find_best_partition(gamma, min_size, max_size, min_pos, max_pos)
    else:
        i0 = min_pos
        i1 = max_pos

        B = [-gamma] * (i1 - i0 + 1)
        p = [0] * (i1 - i0)
        for r in range(i0, i1):
            B[r + 1 - i0] = math.inf
            a = max(r + 1 - max_size, i0)
            b = max(r + 1 - min_size + 1, i0)
            for l in range(a, b):
                b = B[l - i0] + gamma + dist(l, r)
                if b <= B[r + 1 - i0]:
                    B[r + 1 - i0] = b
                    p[r - i0] = l - 1

            mu_dist.cleanup_cache()

    # Routine "Segmentation from partition" in [1]
    # Convert interval representation computed above
    # to a list of intervals and values.
    r = len(p) - 1 + min_pos
    l = p[r - min_pos]
    right = []
    values = []
    dists = []
    while r >= min_pos:
        right.append(r + 1)
        values.append(mu((l + 1), r))
        dists.append(dist((l + 1), r))
        r = l
        l = p[r - min_pos]
    right.reverse()
    values.reverse()
    dists.reverse()

    return right, values, dists


def solve_potts_autogamma(y, w, beta=None, **kw):
    """Solve Potts problem with automatically determined gamma.

    The optimal value is determined by minimizing the information measure::

        f(gamma) = beta J(x(gamma)) + log sum(abs(x(gamma) - y)**p)

    where x(gamma) is the solution to the Potts problem for a fixed
    gamma. The minimization is only performed rather roughly.

    Parameters
    ----------
    beta : float or 'bic'
         Penalty parameter. Default is 4*ln(n)/n, similar to Bayesian
         information criterion for gaussian model with unknown variance
         assuming 4 DOF per breakpoint.

    """
    n = len(y)

    if n == 0:
        return [], [], [], None

    mu_dist = get_mu_dist(y, w)
    dist = mu_dist.dist

    if beta is None:
        beta = 4 * math.log(n) / n

    gamma_0 = dist(0, n - 1)

    if gamma_0 == 0:
        # Zero variance
        gamma_0 = 1.0

    best_r = [None]
    best_v = [None]
    best_d = [None]
    best_obj = [math.inf]
    best_gamma = [None]

    def f(x):
        gamma = gamma_0 * math.exp(x)
        r, v, d = solve_potts_approx(y, w, gamma=gamma, mu_dist=mu_dist, **kw)

        # MLE fit noise correlation
        def sigma_star(rights, values, rho):
            """
            |E_0| + sum_{j>0} |E_j - rho E_{j-1}|
            """
            l = 1
            E_prev = y[0] - values[0]
            s = abs(E_prev) * w[0]
            for r, v in zip(rights, values):
                for yv, wv in zip(y[l:r], w[l:r]):
                    E = yv - v
                    s += abs(E - rho * E_prev) * wv
                    E_prev = E
                l = r
            return s

        rho_best = golden_search(
            lambda rho: sigma_star(r, v, rho), -1, 1, xatol=0.05, expand_bounds=True
        )

        # Measurement noise floor
        if len(v) > 2:
            absdiff = [abs(v[j + 1] - v[j]) for j in range(len(v) - 1)]
            sigma_0 = 0.1 * min(absdiff)
        else:
            absv = [abs(z) for z in v]
            sigma_0 = 0.001 * min(absv)
        sigma_0 = max(1e-300, sigma_0)

        # Objective function
        s = sigma_star(r, v, rho_best)
        obj = beta * len(r) + math.log(sigma_0 + s)

        # Done
        if obj < best_obj[0]:
            best_r[0] = r
            best_v[0] = v
            best_d[0] = d
            best_gamma[0] = gamma
            best_obj[0] = obj
        return obj

    # Try to find best gamma (golden section search on log-scale); we
    # don't need an accurate value for it however
    a = math.log(0.1 / n)
    b = 0.0
    golden_search(f, a, b, xatol=abs(a) * 0.1, ftol=0, expand_bounds=True)
    return best_r[0], best_v[0], best_d[0], best_gamma[0]


def solve_potts_approx(y, w, gamma=None, min_size=1, **kw):
    """
    Fit penalized stepwise constant function (Potts model) to data
    approximately, in linear time.

    Do this by running the exact solver using a small maximum interval
    size, and then combining consecutive intervals together if it
    decreases the cost function.
    """
    n = len(y)

    if n == 0:
        return [], [], []

    mu_dist = kw.get('mu_dist')
    if mu_dist is None:
        mu_dist = get_mu_dist(y, w)
        kw['mu_dist'] = mu_dist

    if gamma is None:
        dist = mu_dist.dist
        gamma = 3 * dist(0, n - 1) * math.log(n) / n

    if min_size < 10:
        max_size = 20
    else:
        max_size = min_size + 50

    right, values, dists = solve_potts(y, w, gamma, min_size=min_size, max_size=max_size, **kw)
    return merge_pieces(gamma, right, values, dists, mu_dist, max_size=max_size)


def merge_pieces(gamma, right, values, dists, mu_dist, max_size):
    """
    Combine consecutive intervals in Potts model solution, if doing
    that reduces the cost function.
    """
    mu, dist = mu_dist.mu, mu_dist.dist

    right = list(right)

    # Combine consecutive intervals, if it results to decrease of cost
    # function
    while True:
        min_change = 0
        min_change_j = len(right)

        l = 0
        for j in range(1, len(right)):
            if min_change_j < j - 2:
                break

            # Check whether merging consecutive intervals results to
            # decrease in the cost function
            change = dist(l, right[j] - 1) - (
                dist(l, right[j - 1] - 1) + dist(right[j - 1], right[j] - 1) + gamma
            )
            if change <= min_change:
                min_change = change
                min_change_j = j - 1
            l = right[j - 1]

        if min_change_j < len(right):
            del right[min_change_j]
        else:
            break

    # Check whether perturbing boundary positions leads to improvement
    # in the cost function. The restricted Potts minimization can
    # return sub-optimal boundaries due to the interval maximum size
    # restriction.
    l = 0
    for j in range(1, len(right)):
        prev_score = dist(l, right[j - 1] - 1) + dist(right[j - 1], right[j] - 1)
        new_off = 0
        for off in range(-max_size, max_size + 1):
            if right[j - 1] + off - 1 < l or right[j - 1] + off > right[j] - 1 or off == 0:
                continue
            new_score = dist(l, right[j - 1] + off - 1) + dist(right[j - 1] + off, right[j] - 1)
            if new_score < prev_score:
                new_off = off
                prev_score = new_score

        if new_off != 0:
            right[j - 1] += new_off

        l = right[j - 1]

    # Rebuild values and dists lists
    l = 0
    values = []
    dists = []
    for j in range(len(right)):
        dists.append(dist(l, right[j] - 1))
        values.append(mu(l, right[j] - 1))
        l = right[j]

    return right, values, dists


class L1Dist:
    r"""
    Fast computations for the L1 distance measures.

    This computes:

    .. code-block::

        mu(l, r) = median(y[l:r+1], weights=w[l:r+1])
        dist(l, r) = sum(w*abs(x - mu(l, r)) for x, w in zip(y[l:r+1], weights[l:r+1]))

    We do not use here an approach that has asymptotically optimal performance;
    at least :math:`O(n^2 * \log(n))` would be achievable, whereas we have here
    :math:`O(n^3)`.  The asymptotic performance does not matter for
    :py:func:`asv.step_detect.solve_potts_approx`, which only looks at small
    windows of the data. It is more important to try to optimize the constant
    prefactors, which for Python means minimal code.

    """

    def __init__(self, y, w):
        self.y = y
        self.w = w

        class mu_dict(collections.defaultdict):
            def __missing__(self, a):
                l, r = a
                v = weighted_median(y[l : r + 1], w[l : r + 1])
                self[a] = v
                return v

        mu = mu_dict()

        class dist_dict(collections.defaultdict):
            def __missing__(self, a):
                l, r = a
                m = mu[l, r]
                v = sum(wx * abs(x - m) for x, wx in zip(y[l : r + 1], w[l : r + 1]))
                self[a] = v
                return v

        self.mu_memo = mu
        self.dist_memo = dist_dict()

    def mu(self, *a):
        return self.mu_memo[a]

    def dist(self, *a):
        return self.dist_memo[a]

    def cleanup_cache(self):
        # Reset cache if it is too big
        if len(self.mu_memo) < 500000:
            return

        self.mu_memo.clear()
        self.dist_memo.clear()


def get_mu_dist(y, w):
    if _rangemedian is not None:
        return _rangemedian.RangeMedian(y, w)
    else:
        return L1Dist(y, w)


def rolling_median_dev(items):
    """
    Compute median(items[:j]), deviation[j]) for j in range(1, len(items))
    in O(n log n) time.

    deviation[j] == sum(abs(x - median(items[:j])) for x in items[:j])
    """
    min_heap = []
    max_heap = []
    min_heap_sum = 0  # equal to -sum(min_heap)
    max_heap_sum = 0  # equal to sum(max_heap)
    s = iter(items)
    try:
        while True:
            # Odd
            v = next(s)
            min_heap_sum += v
            v = -heapq.heappushpop(min_heap, -v)
            min_heap_sum -= v
            heapq.heappush(max_heap, v)
            max_heap_sum += v
            # Ensure d >= 0 despite rounding error
            d = max(0, max_heap_sum - min_heap_sum - max_heap[0])
            yield (max_heap[0], d)

            # Even
            v = next(s)
            max_heap_sum += v
            v = heapq.heappushpop(max_heap, v)
            max_heap_sum -= v
            heapq.heappush(min_heap, -v)
            min_heap_sum += v
            d = max(0, max_heap_sum - min_heap_sum)
            yield ((max_heap[0] - min_heap[0]) / 2, d)
    except StopIteration:
        return


def weighted_median(y, w):
    """
    Compute weighted median of `y` with weights `w`.
    """
    items = sorted(zip(y, w))
    midpoint = sum(w) / 2

    yvals = []
    wsum = 0

    for yy, ww in items:
        wsum += ww
        if wsum > midpoint:
            yvals.append(yy)
            break
        elif wsum == midpoint:
            yvals.append(yy)
    else:
        yvals = y

    return sum(yvals) / len(yvals)


def golden_search(f, a, b, xatol=1e-6, ftol=1e-8, expand_bounds=False):
    """
    Find minimum of a function on interval [a, b]
    using golden section search.

    If expand_bounds=True, expand the interval so that the function is
    first evaluated at x=a and x=b.
    """

    ratio = 2 / (1 + math.sqrt(5))

    if not expand_bounds:
        x0 = a
        x3 = b
    else:
        x0 = (ratio * a - (1 - ratio) * b) / (2 * ratio - 1)
        x3 = (ratio * b - (1 - ratio) * a) / (2 * ratio - 1)

    x1 = ratio * x0 + (1 - ratio) * x3
    x2 = (1 - ratio) * x0 + ratio * x3

    f1 = f(x1)
    f2 = f(x2)

    f0 = max(abs(f1), abs(f2))

    while True:
        if abs(x0 - x3) < xatol or abs(f1 - f2) < ftol * f0:
            break

        if f2 < f1:
            x0 = x1
            x1 = x2
            x2 = ratio * x1 + (1 - ratio) * x3
            f1 = f2
            f2 = f(x2)
        else:
            x3 = x2
            x2 = x1
            x1 = ratio * x2 + (1 - ratio) * x0
            f2 = f1
            f1 = f(x1)

    if f2 < f1:
        return x2
    else:
        return x1


def _plot_potts(x, sol):
    import matplotlib.pyplot as plt
    import numpy as np

    t = np.arange(len(x))

    plt.clf()
    plt.plot(t, x, 'k.')

    l = 0
    for r, v in zip(sol[0], sol[1]):
        plt.plot([l, r - 1], [v, v], 'b-o', hold=1)
        l = r
