# -*- coding: utf-8 -*-

import numpy as np

import pandas as pd
from pandas.core import common as com
from pandas import api
from pandas.api import types
from pandas.util import testing as tm

_multiprocess_can_split_ = True


class Base(object):

    def check(self, namespace, expected, ignored=None):
        # see which names are in the namespace, minus optional
        # ignored ones
        # compare vs the expected

        result = sorted([f for f in dir(namespace) if not f.startswith('_')])
        if ignored is not None:
            result = sorted(list(set(result) - set(ignored)))

        expected = sorted(expected)
        tm.assert_almost_equal(result, expected)


class TestPDApi(Base, tm.TestCase):

    # these are optionally imported based on testing
    # & need to be ignored
    ignored = ['tests', 'rpy', 'locale']

    # top-level sub-packages
    lib = ['api', 'compat', 'computation', 'core',
           'indexes', 'formats', 'pandas',
           'test', 'tools', 'tseries',
           'types', 'util', 'options', 'io']

    # top-level packages that are c-imports, should rename to _*
    # to avoid naming conflicts
    lib_to_rename = ['algos', 'hashtable', 'tslib', 'msgpack', 'sparse',
                     'json', 'lib', 'index', 'parser']

    # these are already deprecated; awaiting removal
    deprecated_modules = ['ols', 'stats', 'datetools']

    # misc
    misc = ['IndexSlice', 'NaT']

    # top-level classes
    classes = ['Categorical', 'CategoricalIndex', 'DataFrame', 'DateOffset',
               'DatetimeIndex', 'ExcelFile', 'ExcelWriter', 'Float64Index',
               'Grouper', 'HDFStore', 'Index', 'Int64Index', 'MultiIndex',
               'Period', 'PeriodIndex', 'RangeIndex',
               'Series', 'SparseArray', 'SparseDataFrame',
               'SparseSeries', 'TimeGrouper', 'Timedelta',
               'TimedeltaIndex', 'Timestamp']

    # these are already deprecated; awaiting removal
    deprecated_classes = ['TimeSeries', 'WidePanel',
                          'SparseTimeSeries', 'Panel4D',
                          'SparseList']

    # these should be deprecated in the future
    deprecated_classes_in_future = ['Term', 'Panel']

    # these should be removed from top-level namespace
    remove_classes_from_top_level_namespace = ['Expr']

    # external modules exposed in pandas namespace
    modules = ['np', 'datetime']

    # top-level functions
    funcs = ['bdate_range', 'concat', 'crosstab', 'cut',
             'date_range', 'eval',
             'factorize', 'get_dummies', 'get_store',
             'infer_freq', 'isnull', 'lreshape',
             'match', 'melt', 'notnull', 'offsets',
             'merge', 'merge_ordered', 'merge_asof',
             'period_range',
             'pivot', 'pivot_table', 'plot_params', 'qcut',
             'scatter_matrix',
             'show_versions', 'timedelta_range', 'unique',
             'value_counts', 'wide_to_long']

    # top-level option funcs
    funcs_option = ['reset_option', 'describe_option', 'get_option',
                    'option_context', 'set_option',
                    'set_eng_float_format']

    # top-level read_* funcs
    funcs_read = ['read_clipboard', 'read_csv', 'read_excel', 'read_fwf',
                  'read_gbq', 'read_hdf', 'read_html', 'read_json',
                  'read_msgpack', 'read_pickle', 'read_sas', 'read_sql',
                  'read_sql_query', 'read_sql_table', 'read_stata',
                  'read_table']

    # top-level to_* funcs
    funcs_to = ['to_datetime', 'to_msgpack',
                'to_numeric', 'to_pickle', 'to_timedelta']

    # these should be deprecated in the future
    deprecated_funcs_in_future = ['pnow', 'groupby', 'info']

    # these are already deprecated; awaiting removal
    deprecated_funcs = ['ewma', 'ewmcorr', 'ewmcov', 'ewmstd', 'ewmvar',
                        'ewmvol', 'expanding_apply', 'expanding_corr',
                        'expanding_count', 'expanding_cov', 'expanding_kurt',
                        'expanding_max', 'expanding_mean', 'expanding_median',
                        'expanding_min', 'expanding_quantile',
                        'expanding_skew', 'expanding_std', 'expanding_sum',
                        'expanding_var', 'fama_macbeth', 'rolling_apply',
                        'rolling_corr', 'rolling_count', 'rolling_cov',
                        'rolling_kurt', 'rolling_max', 'rolling_mean',
                        'rolling_median', 'rolling_min', 'rolling_quantile',
                        'rolling_skew', 'rolling_std', 'rolling_sum',
                        'rolling_var', 'rolling_window', 'ordered_merge']

    def test_api(self):

        self.check(pd,
                   self.lib + self.lib_to_rename + self.misc +
                   self.modules + self.deprecated_modules +
                   self.classes + self.deprecated_classes +
                   self.deprecated_classes_in_future +
                   self.remove_classes_from_top_level_namespace +
                   self.funcs + self.funcs_option +
                   self.funcs_read + self.funcs_to +
                   self.deprecated_funcs +
                   self.deprecated_funcs_in_future,
                   self.ignored)


class TestApi(Base, tm.TestCase):

    allowed = ['tests', 'types']

    def test_api(self):

        self.check(api, self.allowed)


class TestTypes(Base, tm.TestCase):

    allowed = ['is_any_int_dtype', 'is_bool', 'is_bool_dtype',
               'is_categorical', 'is_categorical_dtype', 'is_complex',
               'is_complex_dtype', 'is_datetime64_any_dtype',
               'is_datetime64_dtype', 'is_datetime64_ns_dtype',
               'is_datetime64tz_dtype', 'is_datetimetz', 'is_dtype_equal',
               'is_extension_type', 'is_float', 'is_float_dtype',
               'is_floating_dtype', 'is_int64_dtype', 'is_integer',
               'is_integer_dtype', 'is_number', 'is_numeric_dtype',
               'is_object_dtype', 'is_scalar', 'is_sparse',
               'is_string_dtype',
               'is_timedelta64_dtype', 'is_timedelta64_ns_dtype',
               'is_period', 'is_period_dtype',
               'is_re', 'is_re_compilable',
               'is_dict_like', 'is_iterator',
               'is_list_like', 'is_hashable',
               'is_named_tuple', 'is_sequence',
               'pandas_dtype']

    def test_types(self):

        self.check(types, self.allowed)

    def check_deprecation(self, fold, fnew):
        with tm.assert_produces_warning(DeprecationWarning):
            try:
                result = fold('foo')
                expected = fnew('foo')
                self.assertEqual(result, expected)
            except TypeError:
                self.assertRaises(TypeError,
                                  lambda: fnew('foo'))
            except AttributeError:
                self.assertRaises(AttributeError,
                                  lambda: fnew('foo'))

    def test_deprecation_core_common(self):

        # test that we are in fact deprecating
        # the pandas.core.common introspectors
        for t in self.allowed:
            self.check_deprecation(getattr(com, t), getattr(types, t))

    def test_deprecation_core_common_array_equivalent(self):

        with tm.assert_produces_warning(DeprecationWarning):
            com.array_equivalent(np.array([1, 2]), np.array([1, 2]))

    def test_deprecation_core_common_moved(self):

        # these are in pandas.types.common
        l = ['is_datetime_arraylike',
             'is_datetime_or_timedelta_dtype',
             'is_datetimelike',
             'is_datetimelike_v_numeric',
             'is_datetimelike_v_object',
             'is_datetimetz',
             'is_int_or_datetime_dtype',
             'is_period_arraylike',
             'is_string_like',
             'is_string_like_dtype']

        from pandas.types import common as c
        for t in l:
            self.check_deprecation(getattr(com, t), getattr(c, t))

    def test_removed_from_core_common(self):

        for t in ['is_null_datelike_scalar',
                  'ensure_float']:
            self.assertRaises(AttributeError, lambda: getattr(com, t))


class TestDatetools(tm.TestCase):

    def test_deprecation_access_func(self):
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            pd.datetools.to_datetime('2016-01-01')

    def test_deprecation_access_obj(self):
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            pd.datetools.monthEnd

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
