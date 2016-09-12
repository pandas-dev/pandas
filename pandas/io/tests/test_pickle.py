# pylint: disable=E1101,E1103,W0232

""" manage legacy pickle tests """

import nose
import os

from distutils.version import LooseVersion

import pandas as pd
from pandas import Index
from pandas.compat import u, is_platform_little_endian
import pandas
import pandas.util.testing as tm
from pandas.tseries.offsets import Day, MonthEnd


class TestPickle():
    """
    How to add pickle tests:

    1. Install pandas version intended to output the pickle.

    2. Execute "generate_legacy_storage_files.py" to create the pickle.
    $ python generate_legacy_storage_files.py <output_dir> pickle

    3. Move the created pickle to "data/legacy_pickle/<version>" directory.

    NOTE: TestPickle can't be a subclass of tm.Testcase to use test generator.
    http://stackoverflow.com/questions/6689537/
    nose-test-generators-inside-class
    """
    _multiprocess_can_split_ = True

    def setUp(self):
        from pandas.io.tests.generate_legacy_storage_files import (
            create_pickle_data)
        self.data = create_pickle_data()
        self.path = u('__%s__.pickle' % tm.rands(10))

    def compare_element(self, result, expected, typ, version=None):
        if isinstance(expected, Index):
            tm.assert_index_equal(expected, result)
            return

        if typ.startswith('sp_'):
            comparator = getattr(tm, "assert_%s_equal" % typ)
            comparator(result, expected, exact_indices=False)
        elif typ == 'timestamp':
            if expected is pd.NaT:
                assert result is pd.NaT
            else:
                tm.assert_equal(result, expected)
                tm.assert_equal(result.freq, expected.freq)
        else:
            comparator = getattr(tm, "assert_%s_equal" %
                                 typ, tm.assert_almost_equal)
            comparator(result, expected)

    def compare(self, vf, version):

        # py3 compat when reading py2 pickle
        try:
            data = pandas.read_pickle(vf)
        except (ValueError) as e:
            if 'unsupported pickle protocol:' in str(e):
                # trying to read a py3 pickle in py2
                return
            else:
                raise

        for typ, dv in data.items():
            for dt, result in dv.items():
                try:
                    expected = self.data[typ][dt]
                except (KeyError):
                    if version in ('0.10.1', '0.11.0') and dt == 'reg':
                        break
                    else:
                        raise

                # use a specific comparator
                # if available
                comparator = "compare_{typ}_{dt}".format(typ=typ, dt=dt)
                comparator = getattr(self, comparator, self.compare_element)
                comparator(result, expected, typ, version)
        return data

    def compare_sp_series_ts(self, res, exp, typ, version):
        # SparseTimeSeries integrated into SparseSeries in 0.12.0
        # and deprecated in 0.17.0
        if version and LooseVersion(version) <= "0.12.0":
            tm.assert_sp_series_equal(res, exp, check_series_type=False)
        else:
            tm.assert_sp_series_equal(res, exp)

    def compare_series_ts(self, result, expected, typ, version):
        # GH 7748
        tm.assert_series_equal(result, expected)
        tm.assert_equal(result.index.freq, expected.index.freq)
        tm.assert_equal(result.index.freq.normalize, False)
        tm.assert_series_equal(result > 0, expected > 0)

        # GH 9291
        freq = result.index.freq
        tm.assert_equal(freq + Day(1), Day(2))

        res = freq + pandas.Timedelta(hours=1)
        tm.assert_equal(isinstance(res, pandas.Timedelta), True)
        tm.assert_equal(res, pandas.Timedelta(days=1, hours=1))

        res = freq + pandas.Timedelta(nanoseconds=1)
        tm.assert_equal(isinstance(res, pandas.Timedelta), True)
        tm.assert_equal(res, pandas.Timedelta(days=1, nanoseconds=1))

    def compare_series_dt_tz(self, result, expected, typ, version):
        # 8260
        # dtype is object < 0.17.0
        if LooseVersion(version) < '0.17.0':
            expected = expected.astype(object)
            tm.assert_series_equal(result, expected)
        else:
            tm.assert_series_equal(result, expected)

    def compare_series_cat(self, result, expected, typ, version):
        # Categorical dtype is added in 0.15.0
        # ordered is changed in 0.16.0
        if LooseVersion(version) < '0.15.0':
            tm.assert_series_equal(result, expected, check_dtype=False,
                                   check_categorical=False)
        elif LooseVersion(version) < '0.16.0':
            tm.assert_series_equal(result, expected, check_categorical=False)
        else:
            tm.assert_series_equal(result, expected)

    def compare_frame_dt_mixed_tzs(self, result, expected, typ, version):
        # 8260
        # dtype is object < 0.17.0
        if LooseVersion(version) < '0.17.0':
            expected = expected.astype(object)
            tm.assert_frame_equal(result, expected)
        else:
            tm.assert_frame_equal(result, expected)

    def compare_frame_cat_onecol(self, result, expected, typ, version):
        # Categorical dtype is added in 0.15.0
        # ordered is changed in 0.16.0
        if LooseVersion(version) < '0.15.0':
            tm.assert_frame_equal(result, expected, check_dtype=False,
                                  check_categorical=False)
        elif LooseVersion(version) < '0.16.0':
            tm.assert_frame_equal(result, expected, check_categorical=False)
        else:
            tm.assert_frame_equal(result, expected)

    def compare_frame_cat_and_float(self, result, expected, typ, version):
        self.compare_frame_cat_onecol(result, expected, typ, version)

    def compare_index_period(self, result, expected, typ, version):
        tm.assert_index_equal(result, expected)
        tm.assertIsInstance(result.freq, MonthEnd)
        tm.assert_equal(result.freq, MonthEnd())
        tm.assert_equal(result.freqstr, 'M')
        tm.assert_index_equal(result.shift(2), expected.shift(2))

    def compare_sp_frame_float(self, result, expected, typ, version):
        if LooseVersion(version) <= '0.18.1':
            tm.assert_sp_frame_equal(result, expected, exact_indices=False,
                                     check_dtype=False)
        else:
            tm.assert_sp_frame_equal(result, expected)

    def read_pickles(self, version):
        if not is_platform_little_endian():
            raise nose.SkipTest("known failure on non-little endian")

        pth = tm.get_data_path('legacy_pickle/{0}'.format(str(version)))
        n = 0
        for f in os.listdir(pth):
            vf = os.path.join(pth, f)
            data = self.compare(vf, version)

            if data is None:
                continue
            n += 1
        assert n > 0, 'Pickle files are not tested'

    def test_pickles(self):
        pickle_path = tm.get_data_path('legacy_pickle')
        n = 0
        for v in os.listdir(pickle_path):
            pth = os.path.join(pickle_path, v)
            if os.path.isdir(pth):
                yield self.read_pickles, v
            n += 1
        assert n > 0, 'Pickle files are not tested'

    def test_round_trip_current(self):

        try:
            import cPickle as c_pickle

            def c_pickler(obj, path):
                with open(path, 'wb') as fh:
                    c_pickle.dump(obj, fh, protocol=-1)

            def c_unpickler(path):
                with open(path, 'rb') as fh:
                    fh.seek(0)
                    return c_pickle.load(fh)
        except:
            c_pickler = None
            c_unpickler = None

        import pickle as python_pickle

        def python_pickler(obj, path):
            with open(path, 'wb') as fh:
                python_pickle.dump(obj, fh, protocol=-1)

        def python_unpickler(path):
            with open(path, 'rb') as fh:
                fh.seek(0)
                return python_pickle.load(fh)

        for typ, dv in self.data.items():
            for dt, expected in dv.items():

                for writer in [pd.to_pickle, c_pickler, python_pickler]:
                    if writer is None:
                        continue

                    with tm.ensure_clean(self.path) as path:

                        # test writing with each pickler
                        writer(expected, path)

                        # test reading with each unpickler
                        result = pd.read_pickle(path)
                        self.compare_element(result, expected, typ)

                        if c_unpickler is not None:
                            result = c_unpickler(path)
                            self.compare_element(result, expected, typ)

                        result = python_unpickler(path)
                        self.compare_element(result, expected, typ)

    def test_pickle_v0_14_1(self):

        # we have the name warning
        # 10482
        with tm.assert_produces_warning(UserWarning):
            cat = pd.Categorical(values=['a', 'b', 'c'],
                                 categories=['a', 'b', 'c', 'd'],
                                 name='foobar', ordered=False)
        pickle_path = os.path.join(tm.get_data_path(),
                                   'categorical_0_14_1.pickle')
        # This code was executed once on v0.14.1 to generate the pickle:
        #
        # cat = Categorical(labels=np.arange(3), levels=['a', 'b', 'c', 'd'],
        #                   name='foobar')
        # with open(pickle_path, 'wb') as f: pickle.dump(cat, f)
        #
        tm.assert_categorical_equal(cat, pd.read_pickle(pickle_path))

    def test_pickle_v0_15_2(self):
        # ordered -> _ordered
        # GH 9347

        # we have the name warning
        # 10482
        with tm.assert_produces_warning(UserWarning):
            cat = pd.Categorical(values=['a', 'b', 'c'],
                                 categories=['a', 'b', 'c', 'd'],
                                 name='foobar', ordered=False)
        pickle_path = os.path.join(tm.get_data_path(),
                                   'categorical_0_15_2.pickle')
        # This code was executed once on v0.15.2 to generate the pickle:
        #
        # cat = Categorical(labels=np.arange(3), levels=['a', 'b', 'c', 'd'],
        #                   name='foobar')
        # with open(pickle_path, 'wb') as f: pickle.dump(cat, f)
        #
        tm.assert_categorical_equal(cat, pd.read_pickle(pickle_path))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
