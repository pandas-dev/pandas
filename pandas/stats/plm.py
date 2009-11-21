"""
Linear regression objects for panel data
"""

# pylint: disable-msg=W0231
# pylint: disable-msg=E1103

from __future__ import division
from itertools import izip, starmap

import numpy as np
from scipy import stats

from pandas.core.panel import WidePanel, LongPanel
from pandas.core.matrix import DataMatrix
from pandas.core.series import Series
from pandas.stats.ols import OLS, MovingOLS
from pandas.util.decorators import cache_readonly

import pandas.stats.common as common
import pandas.stats.math as math
import pandas.stats.moments as moments

class PanelOLS(OLS):
    """Implements panel OLS.

    Parameters
    ----------
    y: DataFrame
    x: Dict of DataFrame or WidePanel
    intercept: bool
        True if you want an intercept.  True by default.
    nw_lags: None or int
        Number of Newey-West lags.  None by default.
    nw_overlap: bool
        Whether there are overlaps in the NW lags.  Defaults to False.
    window_type: int
        FULL_SAMPLE, ROLLING, EXPANDING.  FULL_SAMPLE by default.
    window: int
        size of window (for rolling/expanding OLS)
    weights: DataFrame
        Weight for each observation.  The weights are not normalized;
        they're multiplied directly by each observation.
    pool: bool, default True
        Whether to run pooled panel regression
    entity_effects: bool, deafult False
        Whether to account for entity fixed effects
    time_effects: bool, default False
        Whether to account for time fixed effects
    x_effects: list, default None
        List of x's to account for fixed effects
    dropped_dummies: dict
        Key is the name of the variable for the fixed effect.
        Value is the value of that variable for which we drop the dummy.

        For entity fixed effects, key equals 'entity', e.g. {'entity' : 'US'}

        By default, the first item is dropped if one is not specified.
    cluster: int
        ENTITY or TIME, indicating entity/time clustering
    """
    def __init__(self, y, x, weights=None,
                 intercept=True, nw_lags=None, entity_effects=False,
                 time_effects=False, x_effects=None, cluster=None,
                 dropped_dummies=None, verbose=False, nw_overlap=False):
        self._x_orig = x
        self._y_orig = y
        self._weights = weights
        self._intercept = intercept
        self._nw_lags = nw_lags
        self._nw_overlap = nw_overlap
        self._entity_effects = entity_effects
        self._time_effects = time_effects
        self._x_effects = x_effects
        self._dropped_dummies = dropped_dummies or {}
        self._cluster = cluster
        self._verbose = verbose

        (self._x, self._x_trans,
         self._x_filtered, self._y,
         self._y_trans) = self._prepare_data()

        self._x_raw = self._x.values
        self._x_trans_raw = self._x_trans.values
        self._y_trans_raw = self._y_trans.values.squeeze()

        self._index = self._y.major_axis
        self._T = len(self._index)

        self._nobs = len(self._y_trans_raw)

    def log(self, msg):
        if self._verbose:
            print msg

    def _prepare_data(self):
        """Cleans and converts input data into LongPanel classes.

        If time effects is True, then we turn off intercepts and omit an item
        from every (entity and x) fixed effect.

        Otherwise:
           - If we have an intercept, we omit an item from every fixed effect.
           - Else, we omit an item from every fixed effect except one of them.

        The categorical variables will get dropped from x.
        """
        (x, x_filtered, y, weights,
         weights_filt, cat_mapping) = self._filter_data()

        self.log('Adding dummies to X variables')
        x = self._add_dummies(x, cat_mapping)

        self.log('Adding dummies to filtered X variables')
        x_filtered = self._add_dummies(x_filtered, cat_mapping)

        if self._x_effects:
            x = x.filterItems(x.items - self._x_effects)
            x_filtered = x_filtered.filterItems(x_filtered.items
                                                - self._x_effects)

        if self._time_effects:
            x_regressor = x.subtract(x.mean(broadcast=True))
            y_regressor = y.subtract(y.mean(broadcast=True))

        elif self._intercept:
            # only add intercept when no time effects
            self.log('Adding intercept')
            x = x_regressor = add_intercept(x)
            x_filtered = add_intercept(x_filtered)
            y_regressor = y
        else:
            self.log('No intercept added')

            x_regressor = x
            y_regressor = y

        if weights is not None:
            x = x.multiply(weights)
            x_regressor = x_regressor.multiply(weights)
            x_filtered = x_filtered.multiply(weights_filt)
            y = y.multiply(weights)
            y_regressor = y_regressor.multiply(weights)

        return x, x_regressor, x_filtered, y, y_regressor

    def _filter_data(self):
        """

        """
        data, cat_mapping = self._convert_x()
        x_names = data.keys()

        if isinstance(data, LongPanel):
            data = data.toWide()

        elif not isinstance(data, WidePanel):
            data = WidePanel.fromDict(data)

        if self._weights is not None:
            data['__weights__'] = self._weights

        # Filter x's without y (so we can make a prediction)
        filtered = data.toLong()

        # Filter all data together using toLong
        data['__y__'] = self._y_orig
        data_long = data.toLong()

        x_filt = filtered.filterItems(x_names)
        weights_filt = filtered['__weights__'] if self._weights else None

        x = data_long.filterItems(x_names)
        y = data_long['__y__']
        weights = data_long['__weights__'] if self._weights else None

        return x, x_filt, y, weights, weights_filt, cat_mapping

    def _convert_x(self):

        # Converts non-numeric data in x to floats. x_converted is the
        # DataMatrix with converted values, and x_conversion is a dict that
        # provides the reverse mapping.  For example, if 'A' was converted to 0
        # for x named 'variety', then x_conversion['variety'][0] is 'A'.
        x_converted = {}
        x_conversion = {}
        for key, value in self._x_orig.iteritems():
            df = value
            if _is_numeric(df):
                x_converted[key] = df
            else:
                values = df.values
                distinct_values = sorted(set(values.flat))
                x_conversion[key] = dict(enumerate(distinct_values))
                new_values = np.searchsorted(distinct_values, values)
                x_converted[key] = DataMatrix(new_values, index=df.index,
                                              columns=df.columns)

        data = x_converted.copy()

        return data, x_conversion

    def _add_dummies(self, panel, mapping):
        """
        Add entity and / or categorical dummies to input X LongPanel

        Returns
        -------
        LongPanel
        """
        panel = self._add_entity_effects(panel)
        panel = self._add_categorical_dummies(panel, mapping)

        return panel

    def _add_entity_effects(self, panel):
        """
        Add entity dummies to panel

        Returns
        -------
        LongPanel
        """
        if not self._entity_effects:
            return panel

        self.log('-- Adding entity fixed effect dummies')

        dummies = panel.getAxisDummies(axis='minor')

        if not self._use_all_dummies:
            if 'entity' in self._dropped_dummies:
                to_exclude = self._dropped_dummies.get('entity')
            else:
                to_exclude = panel.minor_axis[0]

            if to_exclude not in dummies.items:
                raise Exception('%s not in %s' % (to_exclude,
                                                  dummies.items))

            self.log('-- Excluding dummy for entity: %s' % to_exclude)

            dummies = dummies.filterItems(dummies.items - [to_exclude])

        dummies = dummies.addPrefix('fe_')
        panel = panel.merge(dummies)

        return panel

    def _add_categorical_dummies(self, panel, cat_mappings):
        """
        Add categorical dummies to panel

        Returns
        -------
        LongPanel
        """
        if not self._x_effects:
            return panel

        dropped_dummy = (self._entity_effects and not self._use_all_dummies)

        for effect in self._x_effects:
            self.log('-- Adding fixed effect dummies for %s' % effect)

            dummies = panel.getItemDummies(effect)

            val_map = cat_mappings.get(effect)
            if val_map:
                val_map = dict((v, k) for k, v in val_map.iteritems())

            if dropped_dummy or not self._use_all_dummies:
                if effect in self._dropped_dummies:
                    to_exclude = self._dropped_dummies.get(effect)
                    mapped_name = val_map[to_exclude] if val_map else to_exclude
                else:
                    to_exclude = mapped_name = dummies.items[0]

                if mapped_name not in dummies.items:
                    raise Exception('%s not in %s' % (to_exclude,
                                                      dummies.items))

                self.log('-- Excluding dummy for %s: %s' % (effect, to_exclude))

                dummies = dummies.filterItems(dummies.items - [mapped_name])
                dropped_dummy = True

            dummies = _convertDummies(dummies, cat_mappings.get(effect))
            dummies = dummies.addPrefix('%s_' % effect)
            panel = panel.merge(dummies)

        return panel

    @property
    def _use_all_dummies(self):
        """
        In the case of using an intercept or including time fixed
        effects, completely partitioning the sample would make the X
        not full rank.
        """
        return (not self._intercept and not self._time_effects)

    @cache_readonly
    def _beta_raw(self):
        """Runs the regression and returns the beta."""
        X = self._x_trans_raw
        Y = self._y_trans_raw

        XX = np.dot(X.T, X)
        XY = np.dot(X.T, Y)

        return math.solve(XX, XY)

    @cache_readonly
    def beta(self):
        return Series(self._beta_raw, index=self._x.items)

    @cache_readonly
    def _df_model_raw(self):
        """Returns the raw model degrees of freedom."""
        return self._df_raw - 1

    @cache_readonly
    def _df_resid_raw(self):
        """Returns the raw residual degrees of freedom."""
        return self._nobs - self._df_raw

    @cache_readonly
    def _df_raw(self):
        """Returns the degrees of freedom."""
        df = math.rank(self._x_trans_raw)
        if self._time_effects:
            df += self._total_times

        return df

    @cache_readonly
    def _f_stat_raw(self):
        """Returns the raw f-stat value."""
        return math.calc_f_stat(self._nw_lags, self._r2_raw, self._r2_adj_raw,
                                self._x.items, self._beta_raw,
                                self._var_beta_raw, self._nobs, self.df)

    @cache_readonly
    def _r2_raw(self):
        Y = self._y_trans_raw
        Y_orig = self._y.values
        X = self._x_trans_raw

        resid = Y - np.dot(X, self._beta_raw)
        SS_err = (resid ** 2).sum()

        SS_total = ((Y_orig - np.mean(Y_orig)) ** 2).sum()

        return 1 - SS_err / SS_total

    @cache_readonly
    def _r2_adj_raw(self):
        """Returns the raw r-squared adjusted values."""
        factor = ((self._nobs - 1) / (self._nobs - self._df_raw))
        return 1 - (1 - self._r2_raw) * factor

    @cache_readonly
    def _resid_raw(self):
        Y = self._y_trans.values.squeeze()
        X = self._x_trans.values
        return Y - np.dot(X, self._beta_raw)

    @cache_readonly
    def resid(self):
        return self._unstack_vector(self._resid_raw)

    @cache_readonly
    def _rmse_raw(self):
        """Returns the raw rmse values."""
        X = self._x_trans_raw
        Y = self._y_trans_raw
        resid = Y - np.dot(X, self._beta_raw)
        ss = (resid ** 2).sum()
        return np.sqrt(ss / (self._nobs - self._df_raw))

    @cache_readonly
    def _var_beta_raw(self):
        cluster_axis = None
        if self._cluster == common.TIME:
            cluster_axis = 0
        elif self._cluster == common.ENTITY:
            cluster_axis = 1

        if self._time_effects:
            xx = math.xx_time_effects(self._x, self._y)
        else:
            xx = np.dot(self._x.values.T, self._x.values)

        return math.var_beta_panel(self._y, self._x, self._beta_raw, xx,
                                   self._rmse_raw, cluster_axis, self._nw_lags,
                                   self.nobs, self._df_raw, self._nw_overlap)

    @cache_readonly
    def _y_fitted_raw(self):
        """Returns the raw fitted y values."""
        return np.dot(self._x_filtered.values, self._beta_raw)

    @cache_readonly
    def y_fitted(self):
        return self._unstack_vector(self._y_fitted_raw,
                                    index=self._x_filtered.index)

    def f_test(self, hypothesis):
        """Runs the F test, given a joint hypothesis.  The hypothesis is
        represented by a collection of equations, in the form

        A*x_1+B*x_2=C

        You must provide the coefficients even if they're 1.  No spaces.

        The equations can be passed as either a single string or a
        list of strings.

        Examples:
        o = ols(...)
        o.f_test('1*x1+2*x2=0,1*x3=0')
        o.f_test(['1*x1+2*x2=0','1*x3=0'])
        """

        x_names = self._x.items

        R = []
        r = []

        if isinstance(hypothesis, str):
            eqs = hypothesis.split(',')
        elif isinstance(hypothesis, list):
            eqs = hypothesis
        else:
            raise Exception('hypothesis must be either string or list')
        for equation in eqs:
            row = np.zeros(len(x_names))
            lhs, rhs = equation.split('=')
            for s in lhs.split('+'):
                ss = s.split('*')
                coeff = float(ss[0])
                x_name = ss[1]
                idx = x_names.indexMap[x_name]
                row[idx] = coeff
            rhs = float(rhs)

            R.append(row)
            r.append(rhs)

        R = np.array(R)
        q = len(r)
        r = np.array(r).reshape(q, 1)

        result = math.calc_F(R, r, self._beta_raw, self._var_beta_raw,
                             self.nobs, self.df)

        return common.f_stat_to_dict(result)

    def _unstack_vector(self, vec, index=None):
        if index is None:
            index = self._y_trans.index
        panel = LongPanel(vec.reshape((len(vec), 1)), ['dummy'],
                          index=index)

        return panel.toWide()['dummy']

    def _unstack_y(self, vec):
        unstacked = self._unstack_vector(vec)
        return unstacked.reindex(self.beta.index)

    @cache_readonly
    def _time_obs_count(self):
        # XXX
        return self._y.count()

    @cache_readonly
    def _time_has_obs(self):
        return self._time_obs_count > 0

    @property
    def _total_times(self):
        return self._time_has_obs.sum()

def _convertDummies(dummies, mapping):
    # cleans up the names of the generated dummies
    new_items = []
    for item in dummies.items:
        if not mapping:
            var = '%g' % item if isinstance(item, float) else '%s' % item
            new_items.append(var)
        else:
            # renames the dummies if a conversion dict is provided
            new_items.append(mapping[int(item)])

    dummies = LongPanel(dummies.values, new_items, dummies.index)

    return dummies

def _is_numeric(df):
    for col in df:
        if df[col].dtype.name == 'object':
            return False

    return True

def add_intercept(panel, name='intercept'):
    """
    Add column of ones to input panel

    Parameters
    ----------
    panel: Panel (Long or Wide)
    name: string, default 'intercept']

    Returns
    -------
    New object (same type as input)
    """
    panel = panel.copy()
    panel[name] = 1

    return panel

class MovingPanelOLS(PanelOLS, MovingOLS):
    """Implements rolling/expanding panel OLS.

    Parameters
    ----------
    y: DataFrame
    x: Dict of DataFrame
    intercept: bool
        True if you want an intercept.
    nw_lags: None or int
        Number of Newey-West lags.
    window_type: int
        FULL_SAMPLE, ROLLING, EXPANDING.  FULL_SAMPLE by default.
    window: int
        size of window (for rolling/expanding OLS)
    weights: DataFrame
        Weight for each observation.  The weights are not normalized;
        they're multiplied directly by each observation.
    pool: bool
        Whether to run pooled panel regression.  Defaults to true.
    entity_effects: bool
        Whether to account for entity fixed effects.  Defaults to false.
    time_effects: bool
        Whether to account for time fixed effects.  Defaults to false.
    x_effects: list
        List of x's to account for fixed effects.  Defaults to none.
    dropped_dummies: dict
        Key is the name of the variable for the fixed effect.
        Value is the value of that variable for which we drop the dummy.

        For entity fixed effects, key equals 'entity'.

        By default, the first dummy is dropped if no dummy is specified.
    cluster: int
        ENTITY or TIME, indicating entity/time clustering
    """
    def __init__(self, y, x, weights=None,
                 window_type=common.ROLLING, window=10, min_periods=0,
                 intercept=True,
                 nw_lags=None, nw_overlap=False,
                 entity_effects=False,
                 time_effects=False,
                 x_effects=None,
                 cluster=None,
                 dropped_dummies=None,
                 verbose=False):

        self._args = dict(weights=weights,
                          intercept=intercept,
                          nw_lags=nw_lags,
                          nw_overlap=nw_overlap,
                          entity_effects=entity_effects,
                          time_effects=time_effects,
                          x_effects=x_effects,
                          cluster=cluster,
                          dropped_dummies=dropped_dummies,
                          verbose=verbose)

        PanelOLS.__init__(self, y=y, x=x, **self._args)

        self._window_type = common._get_window_type(window_type)
        self._window = window
        self._min_periods = min_periods

    @cache_readonly
    def beta(self):
        return DataMatrix(self._beta_raw,
                          index=self._result_index,
                          columns=self._x.items)

    @cache_readonly
    def _beta_raw(self):
        """Runs the regression and returns the beta."""
        beta, indices = self._rolling_ols_call

        return beta[indices]

    @cache_readonly
    def rank(self):
        return Series(self._rank_raw, index=self._result_index)

    @cache_readonly
    def _rank_raw(self):
        rank = self._rolling_rank

        return rank[self._valid_indices]

    @cache_readonly
    def _result_index(self):
        return self._index[self._valid_indices]

    @property
    def _valid_indices(self):
        return self._rolling_ols_call[1]

    @cache_readonly
    def _rolling_ols_call(self):
        return self._calc_betas()

    def _calc_betas(self):
        N = len(self._index)
        K = len(self._x.items)

        betas = np.empty((N, K), dtype=float)
        betas[:] = np.NaN

        valid = self._time_has_obs
        enough = self._enough_obs
        window_obs = self._window_time_obs

        # Use transformed (demeaned) Y, X variables
        cum_xx = self._cum_xx(self._x_trans)
        cum_xy = self._cum_xy(self._x_trans, self._y_trans)

        for i in xrange(N):
            # XXX
            if not valid[i] or not enough[i]:
                continue

            xx = cum_xx[i]
            xy = cum_xy[i]
            obs = window_obs[i]
            if self._is_rolling and i >= obs:
                xx = xx - cum_xx[i - obs]
                xy = xy - cum_xy[i - obs]

            betas[i] = math.solve(xx, xy).squeeze()

        have_betas = np.arange(N)[-np.isnan(betas).any(axis=1)]

        return betas, have_betas

    def _cum_xx(self, x):
        dates = x.index.major_axis
        K = len(x.items)
        valid = self._time_has_obs
        cum_xx = []

        last = np.zeros((K, K))
        for i in xrange(len(dates)):
            if not valid[i]:
                cum_xx.append(last)
                continue

            date = dates[i]
            x_slice = x.getValueSlice(date, date)
            xx = last = last + np.dot(x_slice.T, x_slice)
            cum_xx.append(xx)

        return cum_xx

    def _cum_xy(self, x, y):
        dates = x.index.major_axis
        valid = self._time_has_obs
        cum_xy = []

        last = np.zeros((len(x.items), 1))
        for i in xrange(len(dates)):
            if not valid[i]:
                cum_xy.append(last)
                continue

            date = dates[i]
            x_slice = x.getValueSlice(date, date)
            y_slice = y.getValueSlice(date, date)

            xy = last = last + np.dot(x_slice.T, y_slice)
            cum_xy.append(xy)

        return cum_xy

    @property
    def _is_rolling(self):
        return self._window_type == common.ROLLING

    @cache_readonly
    def _rolling_rank(self):
        dates = self._x.index.major_axis

        N = len(dates)
        ranks = np.empty(N, dtype=float)
        ranks[:] = np.NaN

        enough = self._enough_obs
        time_periods = self._window_time_obs

        for i in xrange(N):
            if not enough[i]:
                continue

            if self._is_rolling:
                prior_date = dates[i - time_periods[i] + 1]
            else:
                prior_date = dates[0]

            date = dates[i]
            x_slice = self._x.getValueSlice(prior_date, date)
            ranks[i] = math.rank(x_slice)

        return ranks

    @cache_readonly
    def _df_raw(self):
        """Returns the degrees of freedom."""
        df = self._rolling_rank

        if self._time_effects:
            df += self._window_time_obs

        return df[self._valid_indices]

    @cache_readonly
    def _df_resid_raw(self):
        """Returns the raw residual degrees of freedom."""
        return self._window_nobs - self._df_raw

    @cache_readonly
    def _var_beta_raw(self):
        """Returns the raw covariance of beta."""
        x = self._x
        y = self._y
        dates = x.index.major_axis

        cluster_axis = None
        if self._cluster == common.TIME:
            cluster_axis = 0
        elif self._cluster == common.ENTITY:
            cluster_axis = 1

        time_periods = self._window_time_obs
        window_nobs = self._window_nobs
        rmse = self._rmse_raw
        beta = self._beta_raw
        df = self._df_raw

        if not self._time_effects:
            # Non-transformed X

            cum_xx = self._cum_xx(self._x)

        results = []
        for n, i in enumerate(self._valid_indices):
            obs = time_periods[i]

            if self._is_rolling:
                prior_date = dates[i - obs + 1]
            else:
                prior_date = dates[0]

            date = dates[i]

            x_slice = x.getSlice(prior_date, date)
            y_slice = y.getSlice(prior_date, date)

            if self._time_effects:
                xx = math.xx_time_effects(x_slice, y_slice)
            else:
                xx = cum_xx[i]
                if self._is_rolling and i >= obs:
                    xx = xx - cum_xx[i - time_periods[i]]

            result = math.var_beta_panel(y_slice, x_slice, beta[n], xx, rmse[n],
                                         cluster_axis, self._nw_lags,
                                         window_nobs[n], df[n],
                                         self._nw_overlap)

            results.append(result)

        return np.array(results)

    def f_test(self, hypothesis):
        raise Exception('f_test not supported for rolling/expanding OLS')

    @cache_readonly
    def _f_stat_raw(self):
        """Returns the raw f-stat value."""
        items = self._x.items
        nobs = self._window_nobs
        df = self._df_raw
        df_resid = nobs - df

        # var_beta has not been newey-west adjusted
        if self._nw_lags is None:
            F = self._r2_raw / (self._r2_raw - self._r2_adj_raw)

            q = len(items)
            if 'intercept' in items:
                q -= 1

            def get_result_simple(Fst, d):
                return Fst, (q, d), 1 - stats.f.cdf(Fst, q, d)

            # Compute the P-value for each pair
            result = starmap(get_result_simple, izip(F, df_resid))

            return list(result)

        K = len(items)
        R = np.eye(K)
        r = np.zeros((K, 1))

        intercept = items.indexMap.get('intercept')

        if intercept is not None:
            R = np.concatenate((R[0 : intercept], R[intercept + 1:]))
            r = np.concatenate((r[0 : intercept], r[intercept + 1:]))

        def get_result(beta, vcov, n, d):
            return math.calc_F(R, r, beta, vcov, n, d)

        results = starmap(get_result,
                          izip(self._beta_raw, self._var_beta_raw, nobs, df))

        return list(results)

    @cache_readonly
    def _resid_stats(self):
        Y = self._y_trans
        Y_orig = self._y
        X = self._x_trans
        dates = self._index
        time_periods = self._window_time_obs

        sst = []
        sse = []

        for n, index in enumerate(self._valid_indices):

            if self._is_rolling:
                prior_date = dates[index - time_periods[index] + 1]
            else:
                prior_date = dates[0]

            date = dates[index]

            X_slice = X.getValueSlice(prior_date, date)
            Y_slice = Y.getValueSlice(prior_date, date).squeeze()
            Y_orig_slice = Y_orig.getValueSlice(prior_date, date).squeeze()

            beta_slice = self._beta_raw[n]

            resid = Y_slice - np.dot(X_slice, beta_slice)
            SS_err = (resid ** 2).sum()

            Y_mean = np.mean(Y_orig_slice)
            SS_total = ((Y_orig_slice - Y_mean) ** 2).sum()

            sse.append(SS_err)
            sst.append(SS_total)

        sse = np.array(sse)
        sst = np.array(sst)

        return {
            'sse' : sse,
            'sst' : sst,
        }

    @cache_readonly
    def _rmse_raw(self):
        """Returns the raw rmse values."""
        return np.sqrt(self._resid_stats['sse'] / self._df_resid_raw)

    @cache_readonly
    def _r2_raw(self):
        rs = self._resid_stats
        return 1 - rs['sse'] / rs['sst']

    @cache_readonly
    def _r2_adj_raw(self):
        """Returns the raw r-squared adjusted values."""
        nobs = self._window_nobs
        factors = (nobs - 1) / (nobs - self._df_raw)
        return 1 - (1 - self._r2_raw) * factors

    @cache_readonly
    def _t_stat_raw(self):
        """Returns the raw t-stat value."""
        return np.nan_to_num(self._beta_raw / self._std_err_raw)

    @cache_readonly
    def _p_value_raw(self):
        """Returns the raw p values."""
        get_prob = lambda a, b: 2 * (1 - stats.t.cdf(a, b))

        result = starmap(get_prob,
                         izip(np.fabs(self._t_stat_raw),
                              self._window_nobs - self._df_raw))

        result = np.array(list(result))

        return result

    @cache_readonly
    def _resid_raw(self):
        beta_matrix = self._beta_matrix(lag=0)

        Y = self._y_trans.values.squeeze()
        X = self._x_trans.values
        resid = Y - (X * beta_matrix).sum(1)

        return resid

    @cache_readonly
    def _y_fitted_raw(self):
        x = self._x_raw
        betas = self._beta_matrix(lag=0)
        return (betas * x).sum(1)

    @cache_readonly
    def _y_predict_raw(self):
        """Returns the raw predicted y values."""
        x = self._x_raw
        betas = self._beta_matrix(lag=1)
        return (betas * x).sum(1)

    @cache_readonly
    def resid(self):
        return self._unstack_y(self._resid_raw)

    @cache_readonly
    def y_fitted(self):
        return self._unstack_y(self._y_fitted_raw)

    @cache_readonly
    def y_predict(self):
        """Returns the predicted y values."""
        return self._unstack_y(self._y_predict_raw)

    def _beta_matrix(self, lag=0):
        assert(lag >= 0)

        labels = self._y_trans.index.major_labels - lag
        indexer = self._valid_indices.searchsorted(labels, side='left')

        beta_matrix = self._beta_raw[indexer]
        beta_matrix[labels < 0] = np.NaN

        return beta_matrix

    @cache_readonly
    def _window_nobs_raw(self):
        if self._window_type == common.EXPANDING:
            window = len(self._index)
        else:
            window = self._window

        result = moments.rollingSum(self._time_obs_count, window,
                                    minPeriods=1)

        return result.astype(int)

    @cache_readonly
    def _window_nobs(self):
        return self._window_nobs_raw[self._valid_indices]

    @cache_readonly
    def _window_time_obs(self):
        window_obs = moments.rollingSum(self._time_obs_count > 0,
                                        self._window,
                                        minPeriods=1)

        window_obs[np.isnan(window_obs)] = 0
        return window_obs.astype(int)

    @cache_readonly
    def _enough_obs(self):
        return self._window_nobs_raw >= max(self._min_periods,
                                            len(self._x.items) * 2)

def create_ols_dict(attr):
    def attr_getter(self):
        d = {}
        for k, v in self.results.iteritems():
            result = getattr(v, attr)
            d[k] = result

        return d

    return attr_getter

def create_ols_attr(attr):
    return property(create_ols_dict(attr))

class NonPooledPanelOLS(object):
    """Implements non-pooled panel OLS.

    Parameters
    ----------
    y: DataFrame
    x: Series, DataFrame, or dict of Series
    intercept: bool
        True if you want an intercept.
    nw_lags: None or int
        Number of Newey-West lags.
    window_type: int
        FULL_SAMPLE, ROLLING, EXPANDING.  FULL_SAMPLE by default.
    window: int
        size of window (for rolling/expanding OLS)
    """

    ATTRIBUTES = [
        'beta',
        'df',
        'df_model',
        'df_resid',
        'f_stat',
        'p_value',
        'r2',
        'r2_adj',
        'resid',
        'rmse',
        'std_err',
        'summary_as_matrix',
        't_stat',
        'var_beta',
        'x',
        'y',
        'y_fitted',
        'y_predict'
    ]

    def __init__(self, y, x, window_type=common.FULL_SAMPLE, window=None,
                 intercept=True, nw_lags=None, nw_overlap=False):

        for attr in self.ATTRIBUTES:
            setattr(self.__class__, attr, create_ols_attr(attr))

        results = {}

        for entity in y:
            entity_y = y[entity]

            entity_x = {}
            for x_var in x:
                entity_x[x_var] = x[x_var][entity]

            from pandas.stats.interface import ols
            results[entity] = ols(y=entity_y,
                                  x=entity_x,
                                  window_type=window_type,
                                  window=window,
                                  intercept=intercept,
                                  nw_lags=nw_lags,
                                  nw_overlap=nw_overlap)

        self.results = results
