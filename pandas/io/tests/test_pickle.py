# pylint: disable=E1101,E1103,W0232

""" manage legacy pickle tests """

from datetime import datetime, timedelta
import operator
import pickle as pkl
import nose
import os

import numpy as np
import pandas.util.testing as tm
import pandas as pd
from pandas import Index
from pandas.sparse.tests import test_sparse
from pandas import compat
from pandas.compat import u
from pandas.util.misc import is_little_endian
import pandas

class TestPickle():
    """
    How to add pickle tests:

    1. Install pandas version intended to output the pickle.

    2. Execute "generate_legacy_pkcles.py" to create the pickle.
    $ python generate_legacy_pickles.py <version> <output_dir>

    3. Move the created pickle to "data/legacy_pickle/<version>" directory.

    NOTE: TestPickle can't be a subclass of tm.Testcase to use test generator.
    http://stackoverflow.com/questions/6689537/nose-test-generators-inside-class
    """
    _multiprocess_can_split_ = True

    def setUp(self):
        from pandas.io.tests.generate_legacy_pickles import create_data
        self.data = create_data()
        self.path = u('__%s__.pickle' % tm.rands(10))

    def compare_element(self, typ, result, expected):
        if isinstance(expected,Index):
            tm.assert_index_equal(expected, result)
            return

        if typ.startswith('sp_'):
            comparator = getattr(test_sparse,"assert_%s_equal" % typ)
            comparator(result,expected,exact_indices=False)
        else:
            comparator = getattr(tm,"assert_%s_equal" % typ)
            comparator(result,expected)

    def compare(self, vf):

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

                self.compare_element(typ, result, expected)
        return data

    def read_pickles(self, version):
        if not is_little_endian():
            raise nose.SkipTest("known failure on non-little endian")

        pth = tm.get_data_path('legacy_pickle/{0}'.format(str(version)))
        n = 0
        for f in os.listdir(pth):
            vf = os.path.join(pth, f)
            data = self.compare(vf)

            if data is None:
                continue

            if 'series' in data:
                if 'ts' in data['series']:
                    self._validate_timeseries(data['series']['ts'], self.data['series']['ts'])
                    self._validate_frequency(data['series']['ts'])
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
                        self.compare_element(typ, result, expected)

                        if c_unpickler is not None:
                            result = c_unpickler(path)
                            self.compare_element(typ, result, expected)

                        result = python_unpickler(path)
                        self.compare_element(typ, result, expected)

    def _validate_timeseries(self, pickled, current):
        # GH 7748
        tm.assert_series_equal(pickled, current)
        tm.assert_equal(pickled.index.freq, current.index.freq)
        tm.assert_equal(pickled.index.freq.normalize, False)
        tm.assert_numpy_array_equal(pickled > 0, current > 0)

    def _validate_frequency(self, pickled):
        # GH 9291
        from pandas.tseries.offsets import Day
        freq = pickled.index.freq
        result = freq + Day(1)
        tm.assert_equal(result, Day(2))

        result = freq + pandas.Timedelta(hours=1)
        tm.assert_equal(isinstance(result, pandas.Timedelta), True)
        tm.assert_equal(result, pandas.Timedelta(days=1, hours=1))

        result = freq + pandas.Timedelta(nanoseconds=1)
        tm.assert_equal(isinstance(result, pandas.Timedelta), True)
        tm.assert_equal(result, pandas.Timedelta(days=1, nanoseconds=1))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
