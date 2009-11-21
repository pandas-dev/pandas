"""
Unit test suite for OLS and PanelOLS classes
"""

# pylint: disable-msg=W0212

from __future__ import division

from datetime import datetime
import unittest

from numpy.testing import dec
import numpy as np
import scikits.statsmodels as sm
import scikits.statsmodels.datasets as datasets
from scikits.statsmodels import tools

from pandas.core.panel import LongPanel
from pandas.core.api import DataMatrix, Index, Series
from pandas.stats.api import *
from pandas.stats.plm import NonPooledPanelOLS
from pandas.stats.tests.common import assert_almost_equal, BaseTest

class TestOLS(BaseTest):

    FIELDS = ['beta', 'df', 'df_model', 'df_resid', 'f_stat', 'p_value',
              'r2', 'r2_adj', 'resid', 'rmse', 'std_err', 't_stat',
              'var_beta', 'y_fitted']

    # TODO: Add tests for OLS y predict
    # TODO: Right now we just check for consistency between full-sample and
    # rolling/expanding results of the panel OLS.  We should also cross-check
    # with trusted implementations of panel OLS (e.g. R).
    # TODO: Add tests for non pooled OLS.

    def testOLSWithDatasets(self):
        self.checkDataSet(datasets.ccard.Load(), skip_moving=True)
        self.checkDataSet(datasets.cpunish.Load(), skip_moving=True)
        self.checkDataSet(datasets.longley.Load(), skip_moving=True)
        self.checkDataSet(datasets.stackloss.Load(), skip_moving=True)

        self.checkDataSet(datasets.ccard.Load(), 39, 49) # one col in X all 0s
        self.checkDataSet(datasets.copper.Load())
        self.checkDataSet(datasets.scotland.Load())

    def checkDataSet(self, dataset, start=None, end=None, skip_moving=False):
        exog = dataset.exog[start : end]
        endog = dataset.endog[start : end]
        x = DataMatrix(exog, index=np.arange(exog.shape[0]),
                       columns=np.arange(exog.shape[1]))
        y = Series(endog, index=np.arange(len(endog)))

        self.checkOLS(exog, endog, x, y)

        if not skip_moving:
            self.checkMovingOLS(ROLLING, x, y)
            self.checkMovingOLS(ROLLING, x, y, nw_lags=0)
            self.checkMovingOLS(EXPANDING, x, y, nw_lags=0)
            self.checkMovingOLS(ROLLING, x, y, nw_lags=1)
            self.checkMovingOLS(EXPANDING, x, y, nw_lags=1)
            self.checkMovingOLS(EXPANDING, x, y, nw_lags=1, nw_overlap=True)

    def checkOLS(self, exog, endog, x, y):
        reference = sm.OLS(endog, sm.add_constant(exog)).fit()
        result = ols(y=y, x=x)

        assert_almost_equal(reference.params, result._beta_raw)
        assert_almost_equal(reference.df_model, result._df_model_raw)
        assert_almost_equal(reference.df_resid, result._df_resid_raw)
        assert_almost_equal(reference.fvalue, result._f_stat_raw[0])
        assert_almost_equal(reference.pvalues, result._p_value_raw)
        assert_almost_equal(reference.rsquared, result._r2_raw)
        assert_almost_equal(reference.rsquared_adj, result._r2_adj_raw)
        assert_almost_equal(reference.resid, result._resid_raw)
        assert_almost_equal(reference.bse, result._std_err_raw)
        assert_almost_equal(reference.t(), result._t_stat_raw)
        assert_almost_equal(reference.cov_params(), result._var_beta_raw)
        assert_almost_equal(reference.fittedvalues, result._y_fitted_raw)

        _check_non_raw_results(result)

    def checkMovingOLS(self, window_type, x, y, **kwds):
        window = tools.rank(x.values) + 2

        moving = ols(y=y, x=x, window_type=window_type, window=window,
                     **kwds)

        if isinstance(moving.y, Series):
            index = moving.y.index
        elif isinstance(moving.y, LongPanel):
            index = moving.y.major_axis

        time = len(index)

        reference_last_only = ['resid', 'y_fitted']

        for i in xrange(time - window + 1):
            if window_type == ROLLING:
                start = index[i]
            else:
                start = index[0]

            end = index[i + window - 1]

            x2 = {}
            for k, v in x.iteritems():
                x2[k] = v.truncate(start, end)
            y2 = y.truncate(start, end)

            static = ols(y=y2, x=x2, **kwds)

            self.compare(static, moving, reference_last_only, i)

            # y-predict (just non-null check)
            self.assertTrue(np.isfinite(moving._y_predict_raw).all())

        _check_non_raw_results(moving)

    def compare(self, reference, result, reference_last_only=None,
                result_index=None):
        for field in self.FIELDS:
            attr = '_%s_raw' % field

            ref = getattr(reference, attr)

            if (reference_last_only is not None
                and field in reference_last_only):
                ref = ref[-1]

            res = getattr(result, attr)

            if result_index is not None:
                res = res[result_index]

            assert_almost_equal(ref, res)

class TestPanelOLS(BaseTest):


    FIELDS = ['beta', 'df', 'df_model', 'df_resid', 'f_stat',
              'p_value', 'r2', 'r2_adj', 'rmse', 'std_err',
              't_stat', 'var_beta']

    _other_fields = ['resid', 'y_fitted']

    def testFiltering(self):
        result = ols(y=self.panel_y2, x=self.panel_x2)

        x = result._x
        index = [x.major_axis[i] for i in x.index.major_labels]
        index = Index(sorted(set(index)))
        exp_index = Index([datetime(2000, 1, 1), datetime(2000, 1, 3)])
        self.assertTrue(exp_index.equals(index))

        index = [x.minor_axis[i] for i in x.index.minor_labels]
        index = Index(sorted(set(index)))
        exp_index = Index(['A', 'B'])
        self.assertTrue(exp_index.equals(index))

        x = result._x_filtered
        index = [x.major_axis[i] for i in x.index.major_labels]
        index = Index(sorted(set(index)))
        exp_index = Index([datetime(2000, 1, 1),
                           datetime(2000, 1, 3),
                           datetime(2000, 1, 4)])
        self.assertTrue(exp_index.equals(index))

        assert_almost_equal(result._y.values.flat, [1, 4, 5])

        exp_x = [[6, 14, 1],
                 [9, 17, 1],
                 [30, 48, 1]]
        assert_almost_equal(exp_x, result._x.values)

        exp_x_filtered = [[6, 14, 1], [9, 17, 1], [30, 48, 1], [11, 20, 1],
                          [12, 21, 1]]
        assert_almost_equal(exp_x_filtered, result._x_filtered.values)

        self.assertTrue(result._x_filtered.major_axis.equals(
            result.y_fitted.index))

    def testWithWeights(self):
        data = np.arange(10).reshape((5, 2))
        index = [datetime(2000, 1, 1),
                 datetime(2000, 1, 2),
                 datetime(2000, 1, 3),
                 datetime(2000, 1, 4),
                 datetime(2000, 1, 5)]
        cols = ['A', 'B']
        weights = DataMatrix(data, index=index, columns=cols)

        result = ols(y=self.panel_y2, x=self.panel_x2, weights=weights)

        assert_almost_equal(result._y.values.flat, [0, 16, 25])

        exp_x = [[0, 0, 0],
                 [36, 68, 4],
                 [150, 240, 5]]
        assert_almost_equal(result._x.values, exp_x)

        exp_x_filtered = [[0, 0, 0],
                          [36, 68, 4],
                          [150, 240, 5],
                          [66, 120, 6],
                          [84, 147, 7]]

        assert_almost_equal(result._x_filtered.values, exp_x_filtered)

        # _check_non_raw_results(result)

    def testWithTimeEffects(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, time_effects=True)

        assert_almost_equal(result._y_trans.values.flat, [0, -0.5, 0.5])

        exp_x = [[0, 0], [-10.5, -15.5], [10.5, 15.5]]
        assert_almost_equal(result._x_trans.values, exp_x)

        # _check_non_raw_results(result)

    def testWithEntityEffects(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, entity_effects=True)

        assert_almost_equal(result._y.values.flat, [1, 4, 5])
        exp_x = [[6, 14, 0, 1], [9, 17, 0, 1], [30, 48, 1, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1', 'x2', 'fe_B', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.items))

        # _check_non_raw_results(result)

    def testWithEntityEffectsAndDroppedDummies(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, entity_effects=True,
                     dropped_dummies={'entity' : 'B'})

        assert_almost_equal(result._y.values.flat, [1, 4, 5])
        exp_x = [[6, 14, 1, 1], [9, 17, 1, 1], [30, 48, 0, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1', 'x2', 'fe_A', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.items))

        # _check_non_raw_results(result)

    def testWithXEffects(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, x_effects=['x1'])

        assert_almost_equal(result._y.values.flat, [1, 4, 5])
        exp_x = [[0, 0, 14, 1], [0, 1, 17, 1], [1, 0, 48, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1_30', 'x1_9', 'x2', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.items))

        # _check_non_raw_results(result)

    def testWithXEffectsAndDroppedDummies(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, x_effects=['x1'],
                     dropped_dummies={'x1' : 30})

        assert_almost_equal(result._y.values.flat, [1, 4, 5])
        exp_x = [[1, 0, 14, 1], [0, 1, 17, 1], [0, 0, 48, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1_6', 'x1_9', 'x2', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.items))

        # _check_non_raw_results(result)

    def testWithXEffectsAndConversion(self):
        result = ols(y=self.panel_y3, x=self.panel_x3, x_effects=['x1', 'x2'])

        assert_almost_equal(result._y.values.flat, [1, 2, 3, 4])
        exp_x = [[0, 0, 0, 1, 1], [1, 0, 0, 0, 1], [0, 1, 1, 0, 1],
                 [0, 0, 0, 1, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1_B', 'x1_C', 'x2_2.65', 'x2_3.14', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.items))

        # _check_non_raw_results(result)

    def testWithXEffectsAndConversionAndDroppedDummies(self):
        result = ols(y=self.panel_y3, x=self.panel_x3, x_effects=['x1', 'x2'],
                     dropped_dummies={'x2' : '3.14'})

        assert_almost_equal(result._y.values.flat, [1, 2, 3, 4])
        exp_x = [[0, 0, 0, 0, 1], [1, 0, 1, 0, 1], [0, 1, 0, 1, 1],
                 [0, 0, 0, 0, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1_B', 'x1_C', 'x2_1.59', 'x2_2.65', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.items))

        # _check_non_raw_results(result)

    def testForSeries(self):
        self.checkForSeries(self.series_panel_x, self.series_panel_y,
                            self.series_x, self.series_y)

        self.checkForSeries(self.series_panel_x, self.series_panel_y,
                            self.series_x, self.series_y, nw_lags=0)

        self.checkForSeries(self.series_panel_x, self.series_panel_y,
                            self.series_x, self.series_y, nw_lags=1,
                            nw_overlap=True)

    def testRollingWithWeights(self):
        weights = self.panel_y.copy()

        weights.values = np.random.standard_normal(weights.values.shape)
        self.checkRollingOLS(self.panel_x,
                            self.panel_y, weights=weights)

    def testRolling(self):
        self.checkRollingOLS(self.panel_x, self.panel_y)

    def testRollingWithFixedEffects(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            entity_effects=True)

    def testRollingWithTimeEffects(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            time_effects=True)

    def testRollingWithNeweyWest(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            nw_lags=1)

    def testRollingWithEntityCluster(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            cluster=ENTITY)

    def testRollingWithTimeEffectsAndEntityCluster(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            time_effects=True, cluster=ENTITY)

    def testRollingWithTimeCluster(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            cluster=TIME)

    def testRollingWithNeweyWestAndEntityCluster(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            nw_lags=1, cluster=ENTITY)

    def testRollingWithNeweyWestAndTimeEffectsAndEntityCluster(self):
        self.checkRollingOLS(self.panel_x, self.panel_y,
                            nw_lags=1, cluster=ENTITY, time_effects=True)

    def testExpanding(self):
        self.checkRollingOLS(self.panel_x, self.panel_y, window_type=EXPANDING)

    def testNonPooled(self):
        self.checkNonPooled(y=self.panel_y, x=self.panel_x)
        self.checkNonPooled(y=self.panel_y, x=self.panel_x,
                                    window_type=ROLLING, window=25)

    def checkNonPooled(self, x, y, **kwds):
        # For now, just check that it doesn't crash
        result = ols(y=y, x=x, pool=False, **kwds)
        print result
        for attr in NonPooledPanelOLS.ATTRIBUTES:
            print getattr(result, attr)

    def checkRollingOLS(self, x, y, window_type=ROLLING, **kwds):
        window = 25  # must be larger than rank of x

        moving = ols(y=y, x=x, window_type=window_type,
                     window=window, **kwds)

        if isinstance(moving.y, Series):
            index = moving.y.index
        elif isinstance(moving.y, LongPanel):
            index = moving.y.major_axis

        time_periods = moving._window_time_obs

        for n, i in enumerate(moving._valid_indices):
            if window_type == ROLLING:
                prior_date = index[i - time_periods[i] + 1]
            else:
                prior_date = index[0]

            date = index[i]

            x_iter = {}
            for k, v in x.iteritems():
                x_iter[k] = v.truncate(before=prior_date, after=date)
            y_iter = y.truncate(before=prior_date, after=date)

            static = ols(y=y_iter, x=x_iter, **kwds)

            self.compare(static, moving, event_index=i,
                         result_index=n)

        _check_non_raw_results(moving)

    def checkForSeries(self, x, y, series_x, series_y, **kwds):
        # Consistency check with simple OLS.
        result = ols(y=y, x=x, **kwds)
        reference = ols(y=series_y, x=series_x, **kwds)

        self.compare(reference, result)

    def compare(self, static, moving, event_index=None,
                result_index=None):

        # Check resid if we have a time index specified
        if event_index is not None:
            staticSlice = _period_slice(static, -1)
            movingSlice = _period_slice(moving, event_index)

            ref = static._resid_raw[staticSlice]
            res = moving._resid_raw[movingSlice]

            assert_almost_equal(ref, res)

            ref = static._y_fitted_raw[staticSlice]
            res = moving._y_fitted_raw[movingSlice]

            assert_almost_equal(ref, res)

        # Check y_fitted

        for field in self.FIELDS:
            attr = '_%s_raw' % field

            ref = getattr(static, attr)
            res = getattr(moving, attr)

            if result_index is not None:
                res = res[result_index]

            assert_almost_equal(ref, res)

def _check_non_raw_results(model):
    print model
    print model.resid
    print model.summary_as_matrix
    print model.y_fitted
    print model.y_predict

def _period_slice(panelModel, i):
    index = panelModel._x_trans.index
    period = index.major_axis[i]

    L, R = index.getMajorBounds(period, period)

    return slice(L, R)

if __name__ == '__main__':
    unittest.main()
