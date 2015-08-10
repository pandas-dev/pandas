# pylint: disable=E1101,E1103,W0232

""" manage legacy pickle tests """

from datetime import datetime, timedelta
import operator
import pickle as pkl
import nose
import os

from distutils.version import LooseVersion

import numpy as np
import pandas.util.testing as tm
import pandas as pd
from pandas import Index
from pandas.sparse.tests import test_sparse
from pandas import compat
from pandas.compat import u
from pandas.util.misc import is_little_endian
import pandas
from pandas.tseries.offsets import Day, MonthEnd


class TestPickle():
    """
    How to add pickle tests:

    1. Install pandas version intended to output the pickle.

    2. Execute "generate_legacy_storage_files.py" to create the pickle.
    $ python generate_legacy_storage_files.py <output_dir> pickle

    3. Move the created pickle to "data/legacy_pickle/<version>" directory.

    NOTE: TestPickle can't be a subclass of tm.Testcase to use test generator.
    http://stackoverflow.com/questions/6689537/nose-test-generators-inside-class
    """
    _multiprocess_can_split_ = True

    def setUp(self):
        from pandas.io.tests.generate_legacy_storage_files import create_pickle_data
        self.data = create_pickle_data()
        self.path = u('__%s__.pickle' % tm.rands(10))

    def compare_element(self, result, expected, typ, version=None):
        if isinstance(expected,Index):
            tm.assert_index_equal(expected, result)
            return

        if typ.startswith('sp_'):
            comparator = getattr(test_sparse,"assert_%s_equal" % typ)
            comparator(result,expected,exact_indices=False)
        else:
            comparator = getattr(tm,"assert_%s_equal" % typ,tm.assert_almost_equal)
            comparator(result,expected)

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
                    continue

                # use a specific comparator
                # if available
                comparator = getattr(self,"compare_{typ}_{dt}".format(typ=typ,dt=dt), self.compare_element)
                comparator(result, expected, typ, version)
        return data

    def compare_series_dt_tz(self, result, expected, typ, version):
        # 8260
        # dtype is object < 0.17.0
        if LooseVersion(version) < '0.17.0':
            expected = expected.astype(object)
            tm.assert_series_equal(result, expected)
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

    def read_pickles(self, version):
        if not is_little_endian():
            raise nose.SkipTest("known failure on non-little endian")

        pth = tm.get_data_path('legacy_pickle/{0}'.format(str(version)))
        n = 0
        for f in os.listdir(pth):
            vf = os.path.join(pth, f)
            data = self.compare(vf, version)

            if data is None:
                continue

            if 'series' in data:
                if 'ts' in data['series']:
                    self._validate_timeseries(data['series']['ts'], self.data['series']['ts'])
                    self._validate_frequency(data['series']['ts'])
            if 'index' in data:
                if 'period' in data['index']:
                    self._validate_periodindex(data['index']['period'],
                                               self.data['index']['period'])
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
            def c_pickler(obj,path):
                with open(path,'wb') as fh:
                    c_pickle.dump(obj,fh,protocol=-1)

            def c_unpickler(path):
                with open(path,'rb') as fh:
                    fh.seek(0)
                    return c_pickle.load(fh)
        except:
            c_pickler = None
            c_unpickler = None

        import pickle as python_pickle

        def python_pickler(obj,path):
            with open(path,'wb') as fh:
                python_pickle.dump(obj,fh,protocol=-1)

        def python_unpickler(path):
            with open(path,'rb') as fh:
                fh.seek(0)
                return python_pickle.load(fh)

        for typ, dv in self.data.items():
            for dt, expected in dv.items():

                for writer in [pd.to_pickle, c_pickler, python_pickler ]:
                    if writer is None:
                        continue

                    with tm.ensure_clean(self.path) as path:

                        # test writing with each pickler
                        writer(expected,path)

                        # test reading with each unpickler
                        result = pd.read_pickle(path)
                        self.compare_element(result, expected, typ)

                        if c_unpickler is not None:
                            result = c_unpickler(path)
                            self.compare_element(result, expected, typ)

                        result = python_unpickler(path)
                        self.compare_element(result, expected, typ)

    def _validate_timeseries(self, pickled, current):
        # GH 7748
        tm.assert_series_equal(pickled, current)
        tm.assert_equal(pickled.index.freq, current.index.freq)
        tm.assert_equal(pickled.index.freq.normalize, False)
        tm.assert_numpy_array_equal(pickled > 0, current > 0)

    def _validate_frequency(self, pickled):
        # GH 9291
        freq = pickled.index.freq
        result = freq + Day(1)
        tm.assert_equal(result, Day(2))

        result = freq + pandas.Timedelta(hours=1)
        tm.assert_equal(isinstance(result, pandas.Timedelta), True)
        tm.assert_equal(result, pandas.Timedelta(days=1, hours=1))

        result = freq + pandas.Timedelta(nanoseconds=1)
        tm.assert_equal(isinstance(result, pandas.Timedelta), True)
        tm.assert_equal(result, pandas.Timedelta(days=1, nanoseconds=1))

    def _validate_periodindex(self, pickled, current):
        tm.assert_index_equal(pickled, current)
        tm.assertIsInstance(pickled.freq, MonthEnd)
        tm.assert_equal(pickled.freq, MonthEnd())
        tm.assert_equal(pickled.freqstr, 'M')
        tm.assert_index_equal(pickled.shift(2), current.shift(2))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
