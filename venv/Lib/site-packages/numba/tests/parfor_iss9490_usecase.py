"""
This is a testcase for https://github.com/numba/numba/issues/9490.
The bug is very sensitive to the control-flow and variable uses.
It is impossible to shrink the reproducer in any meaningful way.

The test is also sensitive to PYTHONHASHSEED.
PYTHONHASHSEED=1 will trigger the bug.

Example of traceback:

  File "/numba/parfors/parfor.py", line 2070, in _arrayexpr_to_parfor
    index_vars, loopnests = _mk_parfor_loops(pass_states.typemap, size_vars,
                                             scope, loc)
  File "/numba/parfors/parfor.py", line 1981, in _mk_parfor_loops
    for size_var in size_vars:
TypeError: Failed in nopython mode pipeline (step: convert to parfors)
'NoneType' object is not iterable
"""

import numba
import numpy as np


@numba.jit(nopython=True, parallel=True)
def stable_fit(X, y, threshold=3):
    min_obs = int(X.shape[1] * 1.5)
    beta = np.zeros((X.shape[1], y.shape[1]), dtype=np.float64)
    residuals = np.full_like(y, np.nan)
    stable = np.empty((y.shape[1]))
    for idx in numba.prange(y.shape[1]):
        y_sub = y[:, idx]
        isna = np.isnan(y_sub)
        X_sub = X[~isna]
        y_sub = y_sub[~isna]
        is_stable = False

        # Run until minimum observations
        # or until stability is reached
        for jdx in range(len(y_sub), min_obs - 1, -2):
            # Timeseries gets reduced by two elements
            # each iteration
            y_ = y_sub[-jdx:]
            X_ = X_sub[-jdx:]
            beta_sub = np.linalg.solve(np.dot(X_.T, X_), np.dot(X_.T, y_))
            resid_sub = np.dot(X_, beta_sub) - y_
            # Check for stability
            rmse = np.sqrt(np.mean(resid_sub ** 2))
            first = np.fabs(resid_sub[0]) / rmse < threshold
            last = np.fabs(resid_sub[-1]) / rmse < threshold
            slope = np.fabs(beta_sub[1]) / rmse < threshold
            # Break if stability is reached
            is_stable = slope & first & last
            if is_stable:
                break

        beta[:, idx] = beta_sub
        residuals[-jdx:, idx] = resid_sub
        stable[idx] = is_stable
    return beta, residuals, stable.astype(np.bool_)


def check():
    np.random.seed(0)
    X_shape = (10, 4)
    y_shape = (10, 5)
    X = np.random.random(X_shape)
    y = np.random.random(y_shape)

    got_beta, got_residuals, got_stable = stable_fit(X, y)
    exp_beta, exp_residuals, exp_stable = stable_fit.py_func(X, y)

    np.testing.assert_allclose(got_beta, exp_beta)
    np.testing.assert_allclose(got_residuals, exp_residuals)
    np.testing.assert_allclose(got_stable, exp_stable)


if __name__ == "__main__":
    check()
