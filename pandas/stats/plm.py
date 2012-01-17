"""
Linear regression objects for panel data
"""

# pylint: disable-msg=W0231
# pylint: disable-msg=E1101,E1103

from __future__ import division
import warnings

import numpy as np

from pandas.core.panel import Panel
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from pandas.core.sparse import SparsePanel
from pandas.stats.ols import OLS, MovingOLS
import pandas.stats.common as com
import pandas.stats.math as math
from pandas.util.decorators import cache_readonly

class PanelOLS(OLS):
    """Implements panel OLS.

    See ols function docs
    """
    _panel_model = True

    def __init__(self, y, x, weights=None, intercept=True, nw_lags=None,
                 entity_effects=False, time_effects=False, x_effects=None,
                 cluster=None, dropped_dummies=None, verbose=False,
                 nw_overlap=False):
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
        self._cluster = com._get_cluster_type(cluster)
        self._verbose = verbose

        (self._x, self._x_trans,
         self._x_filtered, self._y,
         self._y_trans) = self._prepare_data()

        self._index = self._x.index.levels[0]

        self._T = len(self._index)

    def log(self, msg):
        if self._verbose: # pragma: no cover
            print msg

    def _prepare_data(self):
        """Cleans and stacks input data into DataFrame objects

        If time effects is True, then we turn off intercepts and omit an item
        from every (entity and x) fixed effect.

        Otherwise:
           - If we have an intercept, we omit an item from every fixed effect.
           - Else, we omit an item from every fixed effect except one of them.

        The categorical variables will get dropped from x.
        """
        (x, x_filtered, y, weights, cat_mapping) = self._filter_data()

        self.log('Adding dummies to X variables')
        x = self._add_dummies(x, cat_mapping)

        self.log('Adding dummies to filtered X variables')
        x_filtered = self._add_dummies(x_filtered, cat_mapping)

        if self._x_effects:
            x = x.drop(self._x_effects, axis=1)
            x_filtered = x_filtered.drop(self._x_effects, axis=1)

        if self._time_effects:
            x_regressor = x.sub(x.mean(level=0), level=0)

            unstacked_y = y.unstack()
            y_regressor = unstacked_y.sub(unstacked_y.mean(1), axis=0).stack()
            y_regressor.index = y.index

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
            assert(y_regressor.index.equals(weights.index))
            assert(x_regressor.index.equals(weights.index))

            rt_weights = np.sqrt(weights)
            y_regressor = y_regressor * rt_weights
            x_regressor = x_regressor.mul(rt_weights, axis=0)

        return x, x_regressor, x_filtered, y, y_regressor

    def _filter_data(self):
        """

        """
        data = self._x_orig
        cat_mapping = {}

        if isinstance(data, DataFrame):
            data = data.to_panel()
        else:
            if isinstance(data, Panel):
                data = data.copy()

            if not isinstance(data, SparsePanel):
                data, cat_mapping = self._convert_x(data)

            if not isinstance(data, Panel):
                data = Panel.from_dict(data, intersect=True)

        x_names = data.items

        if self._weights is not None:
            data['__weights__'] = self._weights

        # Filter x's without y (so we can make a prediction)
        filtered = data.to_frame()

        # Filter all data together using to_frame

        # convert to DataFrame
        y = self._y_orig
        if isinstance(y, Series):
            y = y.unstack()

        data['__y__'] = y
        data_long = data.to_frame()

        x_filt = filtered.filter(x_names)
        x = data_long.filter(x_names)
        y = data_long['__y__']

        if self._weights:
            weights = data_long['__weights__']
        else:
            weights = None

        return x, x_filt, y, weights, cat_mapping

    def _convert_x(self, x):
        # Converts non-numeric data in x to floats. x_converted is the
        # DataFrame with converted values, and x_conversion is a dict that
        # provides the reverse mapping.  For example, if 'A' was converted to 0
        # for x named 'variety', then x_conversion['variety'][0] is 'A'.
        x_converted = {}
        cat_mapping = {}
        # x can be either a dict or a Panel, but in Python 3, dicts don't have
        # .iteritems
        iteritems = getattr(x, 'iteritems', x.items)
        for key, df in iteritems():
            assert(isinstance(df, DataFrame))

            if _is_numeric(df):
                x_converted[key] = df
            else:
                try:
                    df = df.astype(float)
                except (TypeError, ValueError):
                    values = df.values
                    distinct_values = sorted(set(values.flat))
                    cat_mapping[key] = dict(enumerate(distinct_values))
                    new_values = np.searchsorted(distinct_values, values)
                    x_converted[key] = DataFrame(new_values, index=df.index,
                                                  columns=df.columns)

        if len(cat_mapping) == 0:
            x_converted = x

        return x_converted, cat_mapping

    def _add_dummies(self, panel, mapping):
        """
        Add entity and / or categorical dummies to input X DataFrame

        Returns
        -------
        DataFrame
        """
        panel = self._add_entity_effects(panel)
        panel = self._add_categorical_dummies(panel, mapping)

        return panel

    def _add_entity_effects(self, panel):
        """
        Add entity dummies to panel

        Returns
        -------
        DataFrame
        """
        from pandas.core.reshape import make_axis_dummies

        if not self._entity_effects:
            return panel

        self.log('-- Adding entity fixed effect dummies')

        dummies = make_axis_dummies(panel, 'minor')

        if not self._use_all_dummies:
            if 'entity' in self._dropped_dummies:
                to_exclude = str(self._dropped_dummies.get('entity'))
            else:
                to_exclude = dummies.columns[0]

            if to_exclude not in dummies.columns:
                raise Exception('%s not in %s' % (to_exclude,
                                                  dummies.columns))

            self.log('-- Excluding dummy for entity: %s' % to_exclude)

            dummies = dummies.filter(dummies.columns - [to_exclude])

        dummies = dummies.add_prefix('FE_')
        panel = panel.join(dummies)

        return panel

    def _add_categorical_dummies(self, panel, cat_mappings):
        """
        Add categorical dummies to panel

        Returns
        -------
        DataFrame
        """
        from pandas.core.reshape import make_column_dummies

        if not self._x_effects:
            return panel

        dropped_dummy = (self._entity_effects and not self._use_all_dummies)

        for effect in self._x_effects:
            self.log('-- Adding fixed effect dummies for %s' % effect)

            dummies = make_column_dummies(panel, effect, prefix=False)

            val_map = cat_mappings.get(effect)
            if val_map:
                val_map = dict((v, k) for k, v in val_map.iteritems())

            if dropped_dummy or not self._use_all_dummies:
                if effect in self._dropped_dummies:
                    to_exclude = mapped_name = self._dropped_dummies.get(effect)

                    if val_map:
                        mapped_name = val_map[to_exclude]
                else:
                    to_exclude = mapped_name = dummies.columns[0]

                if mapped_name not in dummies.columns: # pragma: no cover
                    raise Exception('%s not in %s' % (to_exclude,
                                                      dummies.columns))

                self.log('-- Excluding dummy for %s: %s' % (effect, to_exclude))

                dummies = dummies.filter(dummies.columns - [mapped_name])
                dropped_dummy = True

            dummies = _convertDummies(dummies, cat_mappings.get(effect))
            dummies = dummies.add_prefix('%s_' % effect)
            panel = panel.join(dummies)

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
        X = self._x_trans.values
        Y = self._y_trans.values.squeeze()

        beta, _, _, _ = np.linalg.lstsq(X, Y)

        return beta

    @cache_readonly
    def beta(self):
        return Series(self._beta_raw, index=self._x.columns)

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
        df = math.rank(self._x_trans.values)
        if self._time_effects:
            df += self._total_times

        return df

    @cache_readonly
    def _r2_raw(self):
        Y = self._y_trans.values.squeeze()
        X = self._x_trans.values

        resid = Y - np.dot(X, self._beta_raw)

        SSE = (resid ** 2).sum()

        if self._use_centered_tss:
            SST = ((Y - np.mean(Y)) ** 2).sum()
        else:
            SST = (Y**2).sum()

        return 1 - SSE / SST

    @property
    def _use_centered_tss(self):
        # has_intercept = np.abs(self._resid_raw.sum()) < _FP_ERR
        return self._intercept or self._entity_effects or self._time_effects

    @cache_readonly
    def _r2_adj_raw(self):
        """Returns the raw r-squared adjusted values."""
        nobs = self._nobs
        factors = (nobs - 1) / (nobs - self._df_raw)
        return 1 - (1 - self._r2_raw) * factors

    @cache_readonly
    def _resid_raw(self):
        Y = self._y.values.squeeze()
        X = self._x.values
        return Y - np.dot(X, self._beta_raw)

    @cache_readonly
    def resid(self):
        return self._unstack_vector(self._resid_raw)

    @cache_readonly
    def _rmse_raw(self):
        """Returns the raw rmse values."""
        # X = self._x.values
        # Y = self._y.values.squeeze()

        X = self._x_trans.values
        Y = self._y_trans.values.squeeze()

        resid = Y - np.dot(X, self._beta_raw)
        ss = (resid ** 2).sum()
        return np.sqrt(ss / (self._nobs - self._df_raw))

    @cache_readonly
    def _var_beta_raw(self):
        cluster_axis = None
        if self._cluster == 'time':
            cluster_axis = 0
        elif self._cluster == 'entity':
            cluster_axis = 1

        x = self._x
        y = self._y

        if self._time_effects:
            xx = _xx_time_effects(x, y)
        else:
            xx = np.dot(x.values.T, x.values)

        return _var_beta_panel(y, x, self._beta_raw, xx,
                               self._rmse_raw, cluster_axis, self._nw_lags,
                               self._nobs, self._df_raw, self._nw_overlap)

    @cache_readonly
    def _y_fitted_raw(self):
        """Returns the raw fitted y values."""
        return np.dot(self._x.values, self._beta_raw)

    @cache_readonly
    def y_fitted(self):
        return self._unstack_vector(self._y_fitted_raw, index=self._x.index)

    def _unstack_vector(self, vec, index=None):
        if index is None:
            index = self._y_trans.index
        panel = DataFrame(vec, index=index, columns=['dummy'])
        return panel.to_panel()['dummy']

    def _unstack_y(self, vec):
        unstacked = self._unstack_vector(vec)
        return unstacked.reindex(self.beta.index)

    @cache_readonly
    def _time_obs_count(self):
        return self._y_trans.count(level=0).values

    @cache_readonly
    def _time_has_obs(self):
        return self._time_obs_count > 0

    @property
    def _nobs(self):
        return len(self._y)

def _convertDummies(dummies, mapping):
    # cleans up the names of the generated dummies
    new_items = []
    for item in dummies.columns:
        if not mapping:
            var = str(item)
            if isinstance(item, float):
                var = '%g' % item

            new_items.append(var)
        else:
            # renames the dummies if a conversion dict is provided
            new_items.append(mapping[int(item)])

    dummies = DataFrame(dummies.values, index=dummies.index,
                        columns=new_items)

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
    panel: Panel / DataFrame
    name: string, default 'intercept']

    Returns
    -------
    New object (same type as input)
    """
    panel = panel.copy()
    panel[name] = 1.

    return panel.consolidate()

class MovingPanelOLS(MovingOLS, PanelOLS):
    """Implements rolling/expanding panel OLS.

    See ols function docs
    """
    _panel_model = True

    def __init__(self, y, x, weights=None,
                 window_type='expanding', window=None,
                 min_periods=None,
                 min_obs=None,
                 intercept=True,
                 nw_lags=None, nw_overlap=False,
                 entity_effects=False,
                 time_effects=False,
                 x_effects=None,
                 cluster=None,
                 dropped_dummies=None,
                 verbose=False):

        self._args = dict(intercept=intercept,
                          nw_lags=nw_lags,
                          nw_overlap=nw_overlap,
                          entity_effects=entity_effects,
                          time_effects=time_effects,
                          x_effects=x_effects,
                          cluster=cluster,
                          dropped_dummies=dropped_dummies,
                          verbose=verbose)

        PanelOLS.__init__(self, y=y, x=x, weights=weights,
                          **self._args)

        self._set_window(window_type, window, min_periods)

        if min_obs is None:
            min_obs = len(self._x.columns) + 1

        self._min_obs = min_obs

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

    def lagged_y_predict(self, lag=1):
        """
        Compute forecast Y value lagging coefficient by input number
        of time periods

        Parameters
        ----------
        lag : int

        Returns
        -------
        DataFrame
        """
        x = self._x.values
        betas = self._beta_matrix(lag=lag)
        return self._unstack_y((betas * x).sum(1))

    @cache_readonly
    def _rolling_ols_call(self):
        return self._calc_betas(self._x_trans, self._y_trans)

    @cache_readonly
    def _df_raw(self):
        """Returns the degrees of freedom."""
        df = self._rolling_rank()

        if self._time_effects:
            df += self._window_time_obs

        return df[self._valid_indices]

    @cache_readonly
    def _var_beta_raw(self):
        """Returns the raw covariance of beta."""
        x = self._x
        y = self._y

        dates = x.index.levels[0]

        cluster_axis = None
        if self._cluster == 'time':
            cluster_axis = 0
        elif self._cluster == 'entity':
            cluster_axis = 1

        nobs = self._nobs
        rmse = self._rmse_raw
        beta = self._beta_raw
        df = self._df_raw
        window = self._window

        if not self._time_effects:
            # Non-transformed X
            cum_xx = self._cum_xx(x)

        results = []
        for n, i in enumerate(self._valid_indices):
            if self._is_rolling and i >= window:
                prior_date = dates[i - window + 1]
            else:
                prior_date = dates[0]

            date = dates[i]

            x_slice = x.truncate(prior_date, date)
            y_slice = y.truncate(prior_date, date)

            if self._time_effects:
                xx = _xx_time_effects(x_slice, y_slice)
            else:
                xx = cum_xx[i]
                if self._is_rolling and i >= window:
                    xx = xx - cum_xx[i - window]

            result = _var_beta_panel(y_slice, x_slice, beta[n], xx, rmse[n],
                                    cluster_axis, self._nw_lags,
                                    nobs[n], df[n], self._nw_overlap)

            results.append(result)

        return np.array(results)

    @cache_readonly
    def _resid_raw(self):
        beta_matrix = self._beta_matrix(lag=0)

        Y = self._y.values.squeeze()
        X = self._x.values
        resid = Y - (X * beta_matrix).sum(1)

        return resid

    @cache_readonly
    def _y_fitted_raw(self):
        x = self._x.values
        betas = self._beta_matrix(lag=0)
        return (betas * x).sum(1)

    @cache_readonly
    def _y_predict_raw(self):
        """Returns the raw predicted y values."""
        x = self._x.values
        betas = self._beta_matrix(lag=1)
        return (betas * x).sum(1)

    def _beta_matrix(self, lag=0):
        assert(lag >= 0)

        index = self._y_trans.index
        major_labels = index.labels[0]
        labels = major_labels - lag
        indexer = self._valid_indices.searchsorted(labels, side='left')

        beta_matrix = self._beta_raw[indexer]
        beta_matrix[labels < self._valid_indices[0]] = np.NaN

        return beta_matrix

    @cache_readonly
    def _enough_obs(self):
        # XXX: what's the best way to determine where to start?
        # TODO: write unit tests for this

        rank_threshold = len(self._x.columns) + 1
        if self._min_obs < rank_threshold: # pragma: no cover
            warnings.warn('min_obs is smaller than rank of X matrix')

        enough_observations = self._nobs_raw >= self._min_obs
        enough_time_periods = self._window_time_obs >= self._min_periods
        return enough_time_periods & enough_observations

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
    y : DataFrame
    x : Series, DataFrame, or dict of Series
    intercept : bool
        True if you want an intercept.
    nw_lags : None or int
        Number of Newey-West lags.
    window_type : {'full_sample', 'rolling', 'expanding'}
        'full_sample' by default
    window : int
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

    def __init__(self, y, x, window_type='full_sample', window=None,
                 min_periods=None, intercept=True, nw_lags=None,
                 nw_overlap=False):

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
                                  min_periods=min_periods,
                                  intercept=intercept,
                                  nw_lags=nw_lags,
                                  nw_overlap=nw_overlap)

        self.results = results


def _var_beta_panel(y, x, beta, xx, rmse, cluster_axis,
                   nw_lags, nobs, df, nw_overlap):
    from pandas.core.frame import group_agg
    xx_inv = math.inv(xx)

    yv = y.values

    if cluster_axis is None:
        if nw_lags is None:
            return xx_inv * (rmse ** 2)
        else:
            resid = yv - np.dot(x.values, beta)
            m = (x.values.T * resid).T

            xeps = math.newey_west(m, nw_lags, nobs, df, nw_overlap)

            return np.dot(xx_inv, np.dot(xeps, xx_inv))
    else:
        Xb = np.dot(x.values, beta).reshape((len(x.values), 1))
        resid = DataFrame(yv[:, None] - Xb, index=y.index, columns=['resid'])

        if cluster_axis == 1:
            x = x.swaplevel(0, 1).sortlevel(0)
            resid = resid.swaplevel(0, 1).sortlevel(0)

        m = group_agg(x.values * resid.values, x.index._bounds,
                      lambda x: np.sum(x, axis=0))

        if nw_lags is None:
            nw_lags = 0

        xox = 0
        for i in range(len(x.index.levels[0])):
            xox += math.newey_west(m[i : i + 1], nw_lags,
                                   nobs, df, nw_overlap)

        return np.dot(xx_inv, np.dot(xox, xx_inv))

def _xx_time_effects(x, y):
    """
    Returns X'X - (X'T) (T'T)^-1 (T'X)
    """
    # X'X
    xx = np.dot(x.values.T, x.values)
    xt = x.sum(level=0).values

    count = y.unstack().count(1).values
    selector = count > 0

    # X'X - (T'T)^-1 (T'X)
    xt = xt[selector]
    count = count[selector]

    return xx - np.dot(xt.T / count, xt)


