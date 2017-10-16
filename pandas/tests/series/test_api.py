# coding=utf-8
# pylint: disable-msg=E1101,W0612
from collections import OrderedDict

import pytest

import numpy as np
import pandas as pd

from pandas import Index, Series, DataFrame, date_range
from pandas.core.indexes.datetimes import Timestamp

from pandas.compat import range
from pandas import compat
import pandas.io.formats.printing as printing
from pandas.util.testing import (assert_series_equal,
                                 ensure_clean)
import pandas.util.testing as tm

from .common import TestData


class SharedWithSparse(object):
    """
    A collection of tests Series and SparseSeries can share.

    In generic tests on this class, use ``self._assert_series_equal()``
    which is implemented in sub-classes.
    """
    def _assert_series_equal(self, left, right):
        """Dispatch to series class dependent assertion"""
        raise NotImplementedError

    def test_scalarop_preserve_name(self):
        result = self.ts * 2
        assert result.name == self.ts.name

    def test_copy_name(self):
        result = self.ts.copy()
        assert result.name == self.ts.name

    def test_copy_index_name_checking(self):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy

        self.ts.index.name = None
        assert self.ts.index.name is None
        assert self.ts is self.ts

        cp = self.ts.copy()
        cp.index.name = 'foo'
        printing.pprint_thing(self.ts.index.name)
        assert self.ts.index.name is None

    def test_append_preserve_name(self):
        result = self.ts[:5].append(self.ts[5:])
        assert result.name == self.ts.name

    def test_binop_maybe_preserve_name(self):
        # names match, preserve
        result = self.ts * self.ts
        assert result.name == self.ts.name
        result = self.ts.mul(self.ts)
        assert result.name == self.ts.name

        result = self.ts * self.ts[:-2]
        assert result.name == self.ts.name

        # names don't match, don't preserve
        cp = self.ts.copy()
        cp.name = 'something else'
        result = self.ts + cp
        assert result.name is None
        result = self.ts.add(cp)
        assert result.name is None

        ops = ['add', 'sub', 'mul', 'div', 'truediv', 'floordiv', 'mod', 'pow']
        ops = ops + ['r' + op for op in ops]
        for op in ops:
            # names match, preserve
            s = self.ts.copy()
            result = getattr(s, op)(s)
            assert result.name == self.ts.name

            # names don't match, don't preserve
            cp = self.ts.copy()
            cp.name = 'changed'
            result = getattr(s, op)(cp)
            assert result.name is None

    def test_combine_first_name(self):
        result = self.ts.combine_first(self.ts[:5])
        assert result.name == self.ts.name

    def test_getitem_preserve_name(self):
        result = self.ts[self.ts > 0]
        assert result.name == self.ts.name

        result = self.ts[[0, 2, 4]]
        assert result.name == self.ts.name

        result = self.ts[5:10]
        assert result.name == self.ts.name

    def test_pickle(self):
        unp_series = self._pickle_roundtrip(self.series)
        unp_ts = self._pickle_roundtrip(self.ts)
        assert_series_equal(unp_series, self.series)
        assert_series_equal(unp_ts, self.ts)

    def _pickle_roundtrip(self, obj):

        with ensure_clean() as path:
            obj.to_pickle(path)
            unpickled = pd.read_pickle(path)
            return unpickled

    def test_argsort_preserve_name(self):
        result = self.ts.argsort()
        assert result.name == self.ts.name

    def test_sort_index_name(self):
        result = self.ts.sort_index(ascending=False)
        assert result.name == self.ts.name

    def test_to_sparse_pass_name(self):
        result = self.ts.to_sparse()
        assert result.name == self.ts.name

    def test_constructor_dict(self):
        d = {'a': 0., 'b': 1., 'c': 2.}
        result = self.series_klass(d)
        expected = self.series_klass(d, index=sorted(d.keys()))
        self._assert_series_equal(result, expected)

        result = self.series_klass(d, index=['b', 'c', 'd', 'a'])
        expected = self.series_klass([1, 2, np.nan, 0],
                                     index=['b', 'c', 'd', 'a'])
        self._assert_series_equal(result, expected)

    def test_constructor_subclass_dict(self):
        data = tm.TestSubDict((x, 10.0 * x) for x in range(10))
        series = self.series_klass(data)
        expected = self.series_klass(dict(compat.iteritems(data)))
        self._assert_series_equal(series, expected)

    def test_constructor_ordereddict(self):
        # GH3283
        data = OrderedDict(
            ('col%s' % i, np.random.random()) for i in range(12))

        series = self.series_klass(data)
        expected = self.series_klass(list(data.values()), list(data.keys()))
        self._assert_series_equal(series, expected)

        # Test with subclass
        class A(OrderedDict):
            pass

        series = self.series_klass(A(data))
        self._assert_series_equal(series, expected)

    def test_constructor_dict_multiindex(self):
        d = {('a', 'a'): 0., ('b', 'a'): 1., ('b', 'c'): 2.}
        _d = sorted(d.items())
        result = self.series_klass(d)
        expected = self.series_klass(
            [x[1] for x in _d],
            index=pd.MultiIndex.from_tuples([x[0] for x in _d]))
        self._assert_series_equal(result, expected)

        d['z'] = 111.
        _d.insert(0, ('z', d['z']))
        result = self.series_klass(d)
        expected = self.series_klass([x[1] for x in _d],
                                     index=pd.Index([x[0] for x in _d],
                                                    tupleize_cols=False))
        result = result.reindex(index=expected.index)
        self._assert_series_equal(result, expected)

    def test_constructor_dict_timedelta_index(self):
        # GH #12169 : Resample category data with timedelta index
        # construct Series from dict as data and TimedeltaIndex as index
        # will result NaN in result Series data
        expected = self.series_klass(
            data=['A', 'B', 'C'],
            index=pd.to_timedelta([0, 10, 20], unit='s')
        )

        result = self.series_klass(
            data={pd.to_timedelta(0, unit='s'): 'A',
                  pd.to_timedelta(10, unit='s'): 'B',
                  pd.to_timedelta(20, unit='s'): 'C'},
            index=pd.to_timedelta([0, 10, 20], unit='s')
        )
        self._assert_series_equal(result, expected)


class TestSeriesMisc(TestData, SharedWithSparse):

    series_klass = Series
    # SharedWithSparse tests use generic, series_klass-agnostic assertion
    _assert_series_equal = staticmethod(tm.assert_series_equal)

    def test_tab_completion(self):
        # GH 9910
        s = Series(list('abcd'))
        # Series of str values should have .str but not .dt/.cat in __dir__
        assert 'str' in dir(s)
        assert 'dt' not in dir(s)
        assert 'cat' not in dir(s)

        # similiarly for .dt
        s = Series(date_range('1/1/2015', periods=5))
        assert 'dt' in dir(s)
        assert 'str' not in dir(s)
        assert 'cat' not in dir(s)

        # Similarly for .cat, but with the twist that str and dt should be
        # there if the categories are of that type first cat and str.
        s = Series(list('abbcd'), dtype="category")
        assert 'cat' in dir(s)
        assert 'str' in dir(s)  # as it is a string categorical
        assert 'dt' not in dir(s)

        # similar to cat and str
        s = Series(date_range('1/1/2015', periods=5)).astype("category")
        assert 'cat' in dir(s)
        assert 'str' not in dir(s)
        assert 'dt' in dir(s)  # as it is a datetime categorical

    def test_not_hashable(self):
        s_empty = Series()
        s = Series([1])
        pytest.raises(TypeError, hash, s_empty)
        pytest.raises(TypeError, hash, s)

    def test_contains(self):
        tm.assert_contains_all(self.ts.index, self.ts)

    def test_iter(self):
        for i, val in enumerate(self.series):
            assert val == self.series[i]

        for i, val in enumerate(self.ts):
            assert val == self.ts[i]

    def test_keys(self):
        # HACK: By doing this in two stages, we avoid 2to3 wrapping the call
        # to .keys() in a list()
        getkeys = self.ts.keys
        assert getkeys() is self.ts.index

    def test_values(self):
        tm.assert_almost_equal(self.ts.values, self.ts, check_dtype=False)

    def test_iteritems(self):
        for idx, val in compat.iteritems(self.series):
            assert val == self.series[idx]

        for idx, val in compat.iteritems(self.ts):
            assert val == self.ts[idx]

        # assert is lazy (genrators don't define reverse, lists do)
        assert not hasattr(self.series.iteritems(), 'reverse')

    def test_items(self):
        for idx, val in self.series.items():
            assert val == self.series[idx]

        for idx, val in self.ts.items():
            assert val == self.ts[idx]

        # assert is lazy (genrators don't define reverse, lists do)
        assert not hasattr(self.series.items(), 'reverse')

    def test_raise_on_info(self):
        s = Series(np.random.randn(10))
        with pytest.raises(AttributeError):
            s.info()

    def test_copy(self):

        for deep in [None, False, True]:
            s = Series(np.arange(10), dtype='float64')

            # default deep is True
            if deep is None:
                s2 = s.copy()
            else:
                s2 = s.copy(deep=deep)

            s2[::2] = np.NaN

            if deep is None or deep is True:
                # Did not modify original Series
                assert np.isnan(s2[0])
                assert not np.isnan(s[0])
            else:
                # we DID modify the original Series
                assert np.isnan(s2[0])
                assert np.isnan(s[0])

        # GH 11794
        # copy of tz-aware
        expected = Series([Timestamp('2012/01/01', tz='UTC')])
        expected2 = Series([Timestamp('1999/01/01', tz='UTC')])

        for deep in [None, False, True]:

            s = Series([Timestamp('2012/01/01', tz='UTC')])

            if deep is None:
                s2 = s.copy()
            else:
                s2 = s.copy(deep=deep)

            s2[0] = pd.Timestamp('1999/01/01', tz='UTC')

            # default deep is True
            if deep is None or deep is True:
                # Did not modify original Series
                assert_series_equal(s2, expected2)
                assert_series_equal(s, expected)
            else:
                # we DID modify the original Series
                assert_series_equal(s2, expected2)
                assert_series_equal(s, expected2)

    def test_axis_alias(self):
        s = Series([1, 2, np.nan])
        assert_series_equal(s.dropna(axis='rows'), s.dropna(axis='index'))
        assert s.dropna().sum('rows') == 3
        assert s._get_axis_number('rows') == 0
        assert s._get_axis_name('rows') == 'index'

    def test_numpy_unique(self):
        # it works!
        np.unique(self.ts)

    def test_ndarray_compat(self):

        # test numpy compat with Series as sub-class of NDFrame
        tsdf = DataFrame(np.random.randn(1000, 3), columns=['A', 'B', 'C'],
                         index=date_range('1/1/2000', periods=1000))

        def f(x):
            return x[x.idxmax()]

        result = tsdf.apply(f)
        expected = tsdf.max()
        tm.assert_series_equal(result, expected)

        # .item()
        s = Series([1])
        result = s.item()
        assert result == 1
        assert s.item() == s.iloc[0]

        # using an ndarray like function
        s = Series(np.random.randn(10))
        result = Series(np.ones_like(s))
        expected = Series(1, index=range(10), dtype='float64')
        tm.assert_series_equal(result, expected)

        # ravel
        s = Series(np.random.randn(10))
        tm.assert_almost_equal(s.ravel(order='F'), s.values.ravel(order='F'))

        # compress
        # GH 6658
        s = Series([0, 1., -1], index=list('abc'))
        result = np.compress(s > 0, s)
        tm.assert_series_equal(result, Series([1.], index=['b']))

        result = np.compress(s < -1, s)
        # result empty Index(dtype=object) as the same as original
        exp = Series([], dtype='float64', index=Index([], dtype='object'))
        tm.assert_series_equal(result, exp)

        s = Series([0, 1., -1], index=[.1, .2, .3])
        result = np.compress(s > 0, s)
        tm.assert_series_equal(result, Series([1.], index=[.2]))

        result = np.compress(s < -1, s)
        # result empty Float64Index as the same as original
        exp = Series([], dtype='float64', index=Index([], dtype='float64'))
        tm.assert_series_equal(result, exp)

    def test_str_attribute(self):
        # GH9068
        methods = ['strip', 'rstrip', 'lstrip']
        s = Series([' jack', 'jill ', ' jesse ', 'frank'])
        for method in methods:
            expected = Series([getattr(str, method)(x) for x in s.values])
            assert_series_equal(getattr(Series.str, method)(s.str), expected)

        # str accessor only valid with string values
        s = Series(range(5))
        with tm.assert_raises_regex(AttributeError,
                                    'only use .str accessor'):
            s.str.repeat(2)

    def test_empty_method(self):
        s_empty = pd.Series()
        assert s_empty.empty

        for full_series in [pd.Series([1]), pd.Series(index=[1])]:
            assert not full_series.empty

    def test_tab_complete_warning(self, ip):
        # https://github.com/pandas-dev/pandas/issues/16409
        pytest.importorskip('IPython', minversion="6.0.0")
        from IPython.core.completer import provisionalcompleter

        code = "import pandas as pd; s = pd.Series()"
        ip.run_code(code)
        with tm.assert_produces_warning(None):
            with provisionalcompleter('ignore'):
                list(ip.Completer.completions('s.', 1))
