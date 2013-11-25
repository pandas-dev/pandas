"""
Unit test suite for OLS and PanelOLS classes
"""

# pylint: disable-msg=W0212

from __future__ import division

from datetime import datetime
from pandas import compat
from distutils.version import LooseVersion
import nose
import numpy as np
from numpy.testing.decorators import slow

from pandas import date_range, bdate_range
from pandas.core.panel import Panel
from pandas import DataFrame, Index, Series, notnull, datetools
from pandas.stats.api import ols
from pandas.stats.ols import _filter_data
from pandas.stats.plm import NonPooledPanelOLS, PanelOLS
from pandas.util.testing import (assert_almost_equal, assert_series_equal,
                                 assert_frame_equal, assertRaisesRegexp)
import pandas.util.testing as tm
import pandas.compat as compat
from .common import BaseTest

_have_statsmodels = True
try:
    import statsmodels.api as sm
except ImportError:
    try:
        import scikits.statsmodels.api as sm
    except ImportError:
        _have_statsmodels = False


def _check_repr(obj):
    repr(obj)
    str(obj)


def _compare_ols_results(model1, model2):
    tm.assert_isinstance(model1, type(model2))

    if hasattr(model1, '_window_type'):
        _compare_moving_ols(model1, model2)
    else:
        _compare_fullsample_ols(model1, model2)


def _compare_fullsample_ols(model1, model2):
    assert_series_equal(model1.beta, model2.beta)


def _compare_moving_ols(model1, model2):
    assert_frame_equal(model1.beta, model2.beta)


class TestOLS(BaseTest):

    _multiprocess_can_split_ = True

    # TODO: Add tests for OLS y predict
    # TODO: Right now we just check for consistency between full-sample and
    # rolling/expanding results of the panel OLS.  We should also cross-check
    # with trusted implementations of panel OLS (e.g. R).
    # TODO: Add tests for non pooled OLS.

    @classmethod
    def setUpClass(cls):
        super(TestOLS, cls).setUpClass()
        try:
            import matplotlib as mpl
            mpl.use('Agg', warn=False)
        except ImportError:
            pass

        if not _have_statsmodels:
            raise nose.SkipTest("no statsmodels")

    def testOLSWithDatasets_ccard(self):
        self.checkDataSet(sm.datasets.ccard.load(), skip_moving=True)
        self.checkDataSet(sm.datasets.cpunish.load(), skip_moving=True)
        self.checkDataSet(sm.datasets.longley.load(), skip_moving=True)
        self.checkDataSet(sm.datasets.stackloss.load(), skip_moving=True)

    @slow
    def testOLSWithDatasets_copper(self):
        self.checkDataSet(sm.datasets.copper.load())

    @slow
    def testOLSWithDatasets_scotland(self):
        self.checkDataSet(sm.datasets.scotland.load())

        # degenerate case fails on some platforms
        # self.checkDataSet(datasets.ccard.load(), 39, 49) # one col in X all
        # 0s

    def testWLS(self):
        # WLS centered SS changed (fixed) in 0.5.0
        sm_version = sm.version.version
        if sm_version < LooseVersion('0.5.0'):
            raise nose.SkipTest("WLS centered SS not fixed in statsmodels"
                                " version {0}".format(sm_version))

        X = DataFrame(np.random.randn(30, 4), columns=['A', 'B', 'C', 'D'])
        Y = Series(np.random.randn(30))
        weights = X.std(1)

        self._check_wls(X, Y, weights)

        weights.ix[[5, 15]] = np.nan
        Y[[2, 21]] = np.nan
        self._check_wls(X, Y, weights)

    def _check_wls(self, x, y, weights):
        result = ols(y=y, x=x, weights=1 / weights)

        combined = x.copy()
        combined['__y__'] = y
        combined['__weights__'] = weights
        combined = combined.dropna()

        endog = combined.pop('__y__').values
        aweights = combined.pop('__weights__').values
        exog = sm.add_constant(combined.values, prepend=False)

        sm_result = sm.WLS(endog, exog, weights=1 / aweights).fit()

        assert_almost_equal(sm_result.params, result._beta_raw)
        assert_almost_equal(sm_result.resid, result._resid_raw)

        self.checkMovingOLS('rolling', x, y, weights=weights)
        self.checkMovingOLS('expanding', x, y, weights=weights)

    def checkDataSet(self, dataset, start=None, end=None, skip_moving=False):
        exog = dataset.exog[start: end]
        endog = dataset.endog[start: end]
        x = DataFrame(exog, index=np.arange(exog.shape[0]),
                      columns=np.arange(exog.shape[1]))
        y = Series(endog, index=np.arange(len(endog)))

        self.checkOLS(exog, endog, x, y)

        if not skip_moving:
            self.checkMovingOLS('rolling', x, y)
            self.checkMovingOLS('rolling', x, y, nw_lags=0)
            self.checkMovingOLS('expanding', x, y, nw_lags=0)
            self.checkMovingOLS('rolling', x, y, nw_lags=1)
            self.checkMovingOLS('expanding', x, y, nw_lags=1)
            self.checkMovingOLS('expanding', x, y, nw_lags=1, nw_overlap=True)

    def checkOLS(self, exog, endog, x, y):
        reference = sm.OLS(endog, sm.add_constant(exog, prepend=False)).fit()
        result = ols(y=y, x=x)

        # check that sparse version is the same
        sparse_result = ols(y=y.to_sparse(), x=x.to_sparse())
        _compare_ols_results(result, sparse_result)

        assert_almost_equal(reference.params, result._beta_raw)
        assert_almost_equal(reference.df_model, result._df_model_raw)
        assert_almost_equal(reference.df_resid, result._df_resid_raw)
        assert_almost_equal(reference.fvalue, result._f_stat_raw[0])
        assert_almost_equal(reference.pvalues, result._p_value_raw)
        assert_almost_equal(reference.rsquared, result._r2_raw)
        assert_almost_equal(reference.rsquared_adj, result._r2_adj_raw)
        assert_almost_equal(reference.resid, result._resid_raw)
        assert_almost_equal(reference.bse, result._std_err_raw)
        assert_almost_equal(reference.tvalues, result._t_stat_raw)
        assert_almost_equal(reference.cov_params(), result._var_beta_raw)
        assert_almost_equal(reference.fittedvalues, result._y_fitted_raw)

        _check_non_raw_results(result)

    def checkMovingOLS(self, window_type, x, y, weights=None, **kwds):
        window = sm.tools.tools.rank(x.values) * 2

        moving = ols(y=y, x=x, weights=weights, window_type=window_type,
                     window=window, **kwds)

        # check that sparse version is the same
        sparse_moving = ols(y=y.to_sparse(), x=x.to_sparse(),
                            weights=weights,
                            window_type=window_type,
                            window=window, **kwds)
        _compare_ols_results(moving, sparse_moving)

        index = moving._index

        for n, i in enumerate(moving._valid_indices):
            if window_type == 'rolling' and i >= window:
                prior_date = index[i - window + 1]
            else:
                prior_date = index[0]

            date = index[i]

            x_iter = {}
            for k, v in compat.iteritems(x):
                x_iter[k] = v.truncate(before=prior_date, after=date)
            y_iter = y.truncate(before=prior_date, after=date)

            static = ols(y=y_iter, x=x_iter, weights=weights, **kwds)

            self.compare(static, moving, event_index=i,
                         result_index=n)

        _check_non_raw_results(moving)

    FIELDS = ['beta', 'df', 'df_model', 'df_resid', 'f_stat', 'p_value',
              'r2', 'r2_adj', 'rmse', 'std_err', 't_stat',
              'var_beta']

    def compare(self, static, moving, event_index=None,
                result_index=None):

        index = moving._index

        # Check resid if we have a time index specified
        if event_index is not None:
            ref = static._resid_raw[-1]

            label = index[event_index]

            res = moving.resid[label]

            assert_almost_equal(ref, res)

            ref = static._y_fitted_raw[-1]
            res = moving.y_fitted[label]

            assert_almost_equal(ref, res)

        # Check y_fitted

        for field in self.FIELDS:
            attr = '_%s_raw' % field

            ref = getattr(static, attr)
            res = getattr(moving, attr)

            if result_index is not None:
                res = res[result_index]

            assert_almost_equal(ref, res)

    def test_ols_object_dtype(self):
        df = DataFrame(np.random.randn(20, 2), dtype=object)
        model = ols(y=df[0], x=df[1])
        summary = repr(model)


class TestOLSMisc(tm.TestCase):

    _multiprocess_can_split_ = True

    '''
    For test coverage with faux data
    '''
    @classmethod
    def setUpClass(cls):
        super(TestOLSMisc, cls).setUpClass()
        if not _have_statsmodels:
            raise nose.SkipTest("no statsmodels")

    def test_f_test(self):
        x = tm.makeTimeDataFrame()
        y = x.pop('A')

        model = ols(y=y, x=x)

        hyp = '1*B+1*C+1*D=0'
        result = model.f_test(hyp)

        hyp = ['1*B=0',
               '1*C=0',
               '1*D=0']
        result = model.f_test(hyp)
        assert_almost_equal(result['f-stat'], model.f_stat['f-stat'])

        self.assertRaises(Exception, model.f_test, '1*A=0')

    def test_r2_no_intercept(self):
        y = tm.makeTimeSeries()
        x = tm.makeTimeDataFrame()

        x_with = x.copy()
        x_with['intercept'] = 1.

        model1 = ols(y=y, x=x)
        model2 = ols(y=y, x=x_with, intercept=False)
        assert_series_equal(model1.beta, model2.beta)

        # TODO: can we infer whether the intercept is there...
        self.assert_(model1.r2 != model2.r2)

        # rolling

        model1 = ols(y=y, x=x, window=20)
        model2 = ols(y=y, x=x_with, window=20, intercept=False)
        assert_frame_equal(model1.beta, model2.beta)
        self.assert_((model1.r2 != model2.r2).all())

    def test_summary_many_terms(self):
        x = DataFrame(np.random.randn(100, 20))
        y = np.random.randn(100)
        model = ols(y=y, x=x)
        model.summary

    def test_y_predict(self):
        y = tm.makeTimeSeries()
        x = tm.makeTimeDataFrame()
        model1 = ols(y=y, x=x)
        assert_series_equal(model1.y_predict, model1.y_fitted)
        assert_almost_equal(model1._y_predict_raw, model1._y_fitted_raw)

    def test_predict(self):
        y = tm.makeTimeSeries()
        x = tm.makeTimeDataFrame()
        model1 = ols(y=y, x=x)
        assert_series_equal(model1.predict(), model1.y_predict)
        assert_series_equal(model1.predict(x=x), model1.y_predict)
        assert_series_equal(model1.predict(beta=model1.beta), model1.y_predict)

        exog = x.copy()
        exog['intercept'] = 1.
        rs = Series(np.dot(exog.values, model1.beta.values), x.index)
        assert_series_equal(model1.y_predict, rs)

        x2 = x.reindex(columns=x.columns[::-1])
        assert_series_equal(model1.predict(x=x2), model1.y_predict)

        x3 = x2 + 10
        pred3 = model1.predict(x=x3)
        x3['intercept'] = 1.
        x3 = x3.reindex(columns=model1.beta.index)
        expected = Series(np.dot(x3.values, model1.beta.values), x3.index)
        assert_series_equal(expected, pred3)

        beta = Series(0., model1.beta.index)
        pred4 = model1.predict(beta=beta)
        assert_series_equal(Series(0., pred4.index), pred4)

    def test_predict_longer_exog(self):
        exogenous = {"1998": "4760", "1999": "5904", "2000": "4504",
                     "2001": "9808", "2002": "4241", "2003": "4086",
                     "2004": "4687", "2005": "7686", "2006": "3740",
                     "2007": "3075", "2008": "3753", "2009": "4679",
                     "2010": "5468", "2011": "7154", "2012": "4292",
                     "2013": "4283", "2014": "4595", "2015": "9194",
                     "2016": "4221", "2017": "4520"}
        endogenous = {"1998": "691", "1999": "1580", "2000": "80",
                      "2001": "1450", "2002": "555", "2003": "956",
                      "2004": "877", "2005": "614", "2006": "468",
                      "2007": "191"}

        endog = Series(endogenous)
        exog = Series(exogenous)
        model = ols(y=endog, x=exog)

        pred = model.y_predict
        self.assert_(pred.index.equals(exog.index))

    def test_longpanel_series_combo(self):
        wp = tm.makePanel()
        lp = wp.to_frame()

        y = lp.pop('ItemA')
        model = ols(y=y, x=lp, entity_effects=True, window=20)
        self.assert_(notnull(model.beta.values).all())
        tm.assert_isinstance(model, PanelOLS)
        model.summary

    def test_series_rhs(self):
        y = tm.makeTimeSeries()
        x = tm.makeTimeSeries()
        model = ols(y=y, x=x)
        expected = ols(y=y, x={'x': x})
        assert_series_equal(model.beta, expected.beta)

        # GH 5233/5250
        assert_series_equal(model.y_predict, model.predict(x=x))

    def test_various_attributes(self):
        # just make sure everything "works". test correctness elsewhere

        x = DataFrame(np.random.randn(100, 5))
        y = np.random.randn(100)
        model = ols(y=y, x=x, window=20)

        series_attrs = ['rank', 'df', 'forecast_mean', 'forecast_vol']

        for attr in series_attrs:
            value = getattr(model, attr)
            tm.assert_isinstance(value, Series)

        # works
        model._results

    def test_catch_regressor_overlap(self):
        df1 = tm.makeTimeDataFrame().ix[:, ['A', 'B']]
        df2 = tm.makeTimeDataFrame().ix[:, ['B', 'C', 'D']]
        y = tm.makeTimeSeries()

        data = {'foo': df1, 'bar': df2}
        self.assertRaises(Exception, ols, y=y, x=data)

    def test_plm_ctor(self):
        y = tm.makeTimeDataFrame()
        x = {'a': tm.makeTimeDataFrame(),
             'b': tm.makeTimeDataFrame()}

        model = ols(y=y, x=x, intercept=False)
        model.summary

        model = ols(y=y, x=Panel(x))
        model.summary

    def test_plm_attrs(self):
        y = tm.makeTimeDataFrame()
        x = {'a': tm.makeTimeDataFrame(),
             'b': tm.makeTimeDataFrame()}

        rmodel = ols(y=y, x=x, window=10)
        model = ols(y=y, x=x)
        model.resid
        rmodel.resid

    def test_plm_lagged_y_predict(self):
        y = tm.makeTimeDataFrame()
        x = {'a': tm.makeTimeDataFrame(),
             'b': tm.makeTimeDataFrame()}

        model = ols(y=y, x=x, window=10)
        result = model.lagged_y_predict(2)

    def test_plm_f_test(self):
        y = tm.makeTimeDataFrame()
        x = {'a': tm.makeTimeDataFrame(),
             'b': tm.makeTimeDataFrame()}

        model = ols(y=y, x=x)

        hyp = '1*a+1*b=0'
        result = model.f_test(hyp)

        hyp = ['1*a=0',
               '1*b=0']
        result = model.f_test(hyp)
        assert_almost_equal(result['f-stat'], model.f_stat['f-stat'])

    def test_plm_exclude_dummy_corner(self):
        y = tm.makeTimeDataFrame()
        x = {'a': tm.makeTimeDataFrame(),
             'b': tm.makeTimeDataFrame()}

        model = ols(
            y=y, x=x, entity_effects=True, dropped_dummies={'entity': 'D'})
        model.summary

        self.assertRaises(Exception, ols, y=y, x=x, entity_effects=True,
                          dropped_dummies={'entity': 'E'})

    def test_columns_tuples_summary(self):
        # #1837
        X = DataFrame(np.random.randn(10, 2), columns=[('a', 'b'), ('c', 'd')])
        Y = Series(np.random.randn(10))

        # it works!
        model = ols(y=Y, x=X)
        model.summary


class TestPanelOLS(BaseTest):

    _multiprocess_can_split_ = True

    FIELDS = ['beta', 'df', 'df_model', 'df_resid', 'f_stat',
              'p_value', 'r2', 'r2_adj', 'rmse', 'std_err',
              't_stat', 'var_beta']

    _other_fields = ['resid', 'y_fitted']

    def testFiltering(self):
        result = ols(y=self.panel_y2, x=self.panel_x2)

        x = result._x
        index = x.index.get_level_values(0)
        index = Index(sorted(set(index)))
        exp_index = Index([datetime(2000, 1, 1), datetime(2000, 1, 3)])
        self.assertTrue
        (exp_index.equals(index))

        index = x.index.get_level_values(1)
        index = Index(sorted(set(index)))
        exp_index = Index(['A', 'B'])
        self.assertTrue(exp_index.equals(index))

        x = result._x_filtered
        index = x.index.get_level_values(0)
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

        exp_x_filtered = [[6, 14, 1],
                          [9, 17, 1],
                          [30, 48, 1],
                          [11, 20, 1],
                          [12, 21, 1]]
        assert_almost_equal(exp_x_filtered, result._x_filtered.values)

        self.assertTrue(result._x_filtered.index.levels[0].equals(
            result.y_fitted.index))

    def test_wls_panel(self):
        y = tm.makeTimeDataFrame()
        x = Panel({'x1': tm.makeTimeDataFrame(),
                   'x2': tm.makeTimeDataFrame()})

        y.ix[[1, 7], 'A'] = np.nan
        y.ix[[6, 15], 'B'] = np.nan
        y.ix[[3, 20], 'C'] = np.nan
        y.ix[[5, 11], 'D'] = np.nan

        stack_y = y.stack()
        stack_x = DataFrame(dict((k, v.stack())
                                 for k, v in compat.iteritems(x)))

        weights = x.std('items')
        stack_weights = weights.stack()

        stack_y.index = stack_y.index._tuple_index
        stack_x.index = stack_x.index._tuple_index
        stack_weights.index = stack_weights.index._tuple_index

        result = ols(y=y, x=x, weights=1 / weights)
        expected = ols(y=stack_y, x=stack_x, weights=1 / stack_weights)

        assert_almost_equal(result.beta, expected.beta)

        for attr in ['resid', 'y_fitted']:
            rvals = getattr(result, attr).stack().values
            evals = getattr(expected, attr).values
            assert_almost_equal(rvals, evals)

    def testWithTimeEffects(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, time_effects=True)

        assert_almost_equal(result._y_trans.values.flat, [0, -0.5, 0.5])

        exp_x = [[0, 0], [-10.5, -15.5], [10.5, 15.5]]
        assert_almost_equal(result._x_trans.values, exp_x)

        # _check_non_raw_results(result)

    def testWithEntityEffects(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, entity_effects=True)

        assert_almost_equal(result._y.values.flat, [1, 4, 5])

        exp_x = DataFrame([[0., 6., 14., 1.], [0, 9, 17, 1], [1, 30, 48, 1]],
                          index=result._x.index, columns=['FE_B', 'x1', 'x2',
                                                          'intercept'],
                          dtype=float)
        tm.assert_frame_equal(result._x, exp_x.ix[:, result._x.columns])
        # _check_non_raw_results(result)

    def testWithEntityEffectsAndDroppedDummies(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, entity_effects=True,
                     dropped_dummies={'entity': 'B'})

        assert_almost_equal(result._y.values.flat, [1, 4, 5])
        exp_x = DataFrame([[1., 6., 14., 1.], [1, 9, 17, 1], [0, 30, 48, 1]],
                          index=result._x.index, columns=['FE_A', 'x1', 'x2',
                                                          'intercept'],
                          dtype=float)
        tm.assert_frame_equal(result._x, exp_x.ix[:, result._x.columns])
        # _check_non_raw_results(result)

    def testWithXEffects(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, x_effects=['x1'])

        assert_almost_equal(result._y.values.flat, [1, 4, 5])

        res = result._x
        exp_x = DataFrame([[0., 0., 14., 1.], [0, 1, 17, 1], [1, 0, 48, 1]],
                          columns=['x1_30', 'x1_9', 'x2', 'intercept'],
                          index=res.index, dtype=float)
        assert_frame_equal(res, exp_x.reindex(columns=res.columns))

    def testWithXEffectsAndDroppedDummies(self):
        result = ols(y=self.panel_y2, x=self.panel_x2, x_effects=['x1'],
                     dropped_dummies={'x1': 30})

        res = result._x
        assert_almost_equal(result._y.values.flat, [1, 4, 5])
        exp_x = DataFrame([[1., 0., 14., 1.], [0, 1, 17, 1], [0, 0, 48, 1]],
                          columns=['x1_6', 'x1_9', 'x2', 'intercept'],
                          index=res.index, dtype=float)

        assert_frame_equal(res, exp_x.reindex(columns=res.columns))

    def testWithXEffectsAndConversion(self):
        result = ols(y=self.panel_y3, x=self.panel_x3, x_effects=['x1', 'x2'])

        assert_almost_equal(result._y.values.flat, [1, 2, 3, 4])
        exp_x = [[0, 0, 0, 1, 1], [1, 0, 0, 0, 1], [0, 1, 1, 0, 1],
                 [0, 0, 0, 1, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1_B', 'x1_C', 'x2_baz', 'x2_foo', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.columns))

        # _check_non_raw_results(result)

    def testWithXEffectsAndConversionAndDroppedDummies(self):
        result = ols(y=self.panel_y3, x=self.panel_x3, x_effects=['x1', 'x2'],
                     dropped_dummies={'x2': 'foo'})

        assert_almost_equal(result._y.values.flat, [1, 2, 3, 4])
        exp_x = [[0, 0, 0, 0, 1], [1, 0, 1, 0, 1], [0, 1, 0, 1, 1],
                 [0, 0, 0, 0, 1]]
        assert_almost_equal(result._x.values, exp_x)

        exp_index = Index(['x1_B', 'x1_C', 'x2_bar', 'x2_baz', 'intercept'])
        self.assertTrue(exp_index.equals(result._x.columns))

        # _check_non_raw_results(result)

    def testForSeries(self):
        self.checkForSeries(self.series_panel_x, self.series_panel_y,
                            self.series_x, self.series_y)

        self.checkForSeries(self.series_panel_x, self.series_panel_y,
                            self.series_x, self.series_y, nw_lags=0)

        self.checkForSeries(self.series_panel_x, self.series_panel_y,
                            self.series_x, self.series_y, nw_lags=1,
                            nw_overlap=True)

    def testRolling(self):
        self.checkMovingOLS(self.panel_x, self.panel_y)

    def testRollingWithFixedEffects(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            entity_effects=True)
        self.checkMovingOLS(self.panel_x, self.panel_y, intercept=False,
                            entity_effects=True)

    def testRollingWithTimeEffects(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            time_effects=True)

    def testRollingWithNeweyWest(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            nw_lags=1)

    def testRollingWithEntityCluster(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            cluster='entity')
    def testUnknownClusterRaisesValueError(self):
        assertRaisesRegexp(ValueError, "Unrecognized cluster.*ridiculous",
                           self.checkMovingOLS, self.panel_x, self.panel_y,
                                               cluster='ridiculous')
    def testRollingWithTimeEffectsAndEntityCluster(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            time_effects=True, cluster='entity')

    def testRollingWithTimeCluster(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            cluster='time')

    def testRollingWithNeweyWestAndEntityCluster(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            nw_lags=1, cluster='entity')

    def testRollingWithNeweyWestAndTimeEffectsAndEntityCluster(self):
        self.checkMovingOLS(self.panel_x, self.panel_y,
                            nw_lags=1, cluster='entity',
                            time_effects=True)

    def testExpanding(self):
        self.checkMovingOLS(
            self.panel_x, self.panel_y, window_type='expanding')

    def testNonPooled(self):
        self.checkNonPooled(y=self.panel_y, x=self.panel_x)
        self.checkNonPooled(y=self.panel_y, x=self.panel_x,
                            window_type='rolling', window=25, min_periods=10)
    def testUnknownWindowType(self):
        assertRaisesRegexp(ValueError, "window.*ridiculous",
                           self.checkNonPooled, y=self.panel_y, x=self.panel_x,
                           window_type='ridiculous', window=25, min_periods=10)

    def checkNonPooled(self, x, y, **kwds):
        # For now, just check that it doesn't crash
        result = ols(y=y, x=x, pool=False, **kwds)

        _check_repr(result)
        for attr in NonPooledPanelOLS.ATTRIBUTES:
            _check_repr(getattr(result, attr))

    def checkMovingOLS(self, x, y, window_type='rolling', **kwds):
        window = 25  # must be larger than rank of x

        moving = ols(y=y, x=x, window_type=window_type,
                     window=window, **kwds)

        index = moving._index

        for n, i in enumerate(moving._valid_indices):
            if window_type == 'rolling' and i >= window:
                prior_date = index[i - window + 1]
            else:
                prior_date = index[0]

            date = index[i]

            x_iter = {}
            for k, v in compat.iteritems(x):
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

    def test_auto_rolling_window_type(self):
        data = tm.makeTimeDataFrame()
        y = data.pop('A')

        window_model = ols(y=y, x=data, window=20, min_periods=10)
        rolling_model = ols(y=y, x=data, window=20, min_periods=10,
                            window_type='rolling')

        assert_frame_equal(window_model.beta, rolling_model.beta)


def _check_non_raw_results(model):
    _check_repr(model)
    _check_repr(model.resid)
    _check_repr(model.summary_as_matrix)
    _check_repr(model.y_fitted)
    _check_repr(model.y_predict)


def _period_slice(panelModel, i):
    index = panelModel._x_trans.index
    period = index.levels[0][i]

    L, R = index.get_major_bounds(period, period)

    return slice(L, R)


class TestOLSFilter(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        date_index = date_range(datetime(2009, 12, 11), periods=3,
                                freq=datetools.bday)
        ts = Series([3, 1, 4], index=date_index)
        self.TS1 = ts

        date_index = date_range(datetime(2009, 12, 11), periods=5,
                                freq=datetools.bday)
        ts = Series([1, 5, 9, 2, 6], index=date_index)
        self.TS2 = ts

        date_index = date_range(datetime(2009, 12, 11), periods=3,
                                freq=datetools.bday)
        ts = Series([5, np.nan, 3], index=date_index)
        self.TS3 = ts

        date_index = date_range(datetime(2009, 12, 11), periods=5,
                                freq=datetools.bday)
        ts = Series([np.nan, 5, 8, 9, 7], index=date_index)
        self.TS4 = ts

        data = {'x1': self.TS2, 'x2': self.TS4}
        self.DF1 = DataFrame(data=data)

        data = {'x1': self.TS2, 'x2': self.TS4}
        self.DICT1 = data

    def testFilterWithSeriesRHS(self):
        (lhs, rhs, weights, rhs_pre,
         index, valid) = _filter_data(self.TS1, {'x1': self.TS2}, None)
        self.tsAssertEqual(self.TS1, lhs)
        self.tsAssertEqual(self.TS2[:3], rhs['x1'])
        self.tsAssertEqual(self.TS2, rhs_pre['x1'])

    def testFilterWithSeriesRHS2(self):
        (lhs, rhs, weights, rhs_pre,
         index, valid) = _filter_data(self.TS2, {'x1': self.TS1}, None)
        self.tsAssertEqual(self.TS2[:3], lhs)
        self.tsAssertEqual(self.TS1, rhs['x1'])
        self.tsAssertEqual(self.TS1, rhs_pre['x1'])

    def testFilterWithSeriesRHS3(self):
        (lhs, rhs, weights, rhs_pre,
         index, valid) = _filter_data(self.TS3, {'x1': self.TS4}, None)
        exp_lhs = self.TS3[2:3]
        exp_rhs = self.TS4[2:3]
        exp_rhs_pre = self.TS4[1:]
        self.tsAssertEqual(exp_lhs, lhs)
        self.tsAssertEqual(exp_rhs, rhs['x1'])
        self.tsAssertEqual(exp_rhs_pre, rhs_pre['x1'])

    def testFilterWithDataFrameRHS(self):
        (lhs, rhs, weights, rhs_pre,
         index, valid) = _filter_data(self.TS1, self.DF1, None)
        exp_lhs = self.TS1[1:]
        exp_rhs1 = self.TS2[1:3]
        exp_rhs2 = self.TS4[1:3]
        self.tsAssertEqual(exp_lhs, lhs)
        self.tsAssertEqual(exp_rhs1, rhs['x1'])
        self.tsAssertEqual(exp_rhs2, rhs['x2'])

    def testFilterWithDictRHS(self):
        (lhs, rhs, weights, rhs_pre,
         index, valid) = _filter_data(self.TS1, self.DICT1, None)
        exp_lhs = self.TS1[1:]
        exp_rhs1 = self.TS2[1:3]
        exp_rhs2 = self.TS4[1:3]
        self.tsAssertEqual(exp_lhs, lhs)
        self.tsAssertEqual(exp_rhs1, rhs['x1'])
        self.tsAssertEqual(exp_rhs2, rhs['x2'])

    def tsAssertEqual(self, ts1, ts2):
        self.assert_(np.array_equal(ts1, ts2))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
