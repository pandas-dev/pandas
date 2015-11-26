import nose

import os
import datetime
import numpy as np
import sys
from distutils.version import LooseVersion

from pandas import compat
from pandas.compat import u
from pandas import (Series, DataFrame, Panel, MultiIndex, bdate_range,
                    date_range, period_range, Index, SparseSeries, SparseDataFrame,
                    SparsePanel)
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean, assert_index_equal
from pandas.tests.test_series import assert_series_equal
from pandas.tests.test_frame import assert_frame_equal
from pandas.tests.test_panel import assert_panel_equal

import pandas
from pandas.sparse.tests.test_sparse import assert_sp_series_equal, assert_sp_frame_equal
from pandas import Timestamp, tslib

nan = np.nan

from pandas.io.packers import to_msgpack, read_msgpack

_multiprocess_can_split_ = False


def check_arbitrary(a, b):

    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        assert(len(a) == len(b))
        for a_, b_ in zip(a, b):
            check_arbitrary(a_, b_)
    elif isinstance(a, Panel):
        assert_panel_equal(a, b)
    elif isinstance(a, DataFrame):
        assert_frame_equal(a, b)
    elif isinstance(a, Series):
        assert_series_equal(a, b)
    elif isinstance(a, Index):
        assert_index_equal(a, b)
    else:
        assert(a == b)

class TestPackers(tm.TestCase):

    def setUp(self):
        self.path = '__%s__.msg' % tm.rands(10)

    def tearDown(self):
        pass

    def encode_decode(self, x, compress=None, **kwargs):
        with ensure_clean(self.path) as p:
            to_msgpack(p, x, compress=compress, **kwargs)
            return read_msgpack(p, **kwargs)

class TestAPI(TestPackers):

    def test_string_io(self):

        df = DataFrame(np.random.randn(10,2))
        s = df.to_msgpack(None)
        result = read_msgpack(s)
        tm.assert_frame_equal(result,df)

        s = df.to_msgpack()
        result = read_msgpack(s)
        tm.assert_frame_equal(result,df)

        s = df.to_msgpack()
        result = read_msgpack(compat.BytesIO(s))
        tm.assert_frame_equal(result,df)

        s = to_msgpack(None,df)
        result = read_msgpack(s)
        tm.assert_frame_equal(result, df)

        with ensure_clean(self.path) as p:

            s = df.to_msgpack()
            fh = open(p,'wb')
            fh.write(s)
            fh.close()
            result = read_msgpack(p)
            tm.assert_frame_equal(result, df)

    def test_iterator_with_string_io(self):

        dfs = [ DataFrame(np.random.randn(10,2)) for i in range(5) ]
        s = to_msgpack(None,*dfs)
        for i, result in enumerate(read_msgpack(s,iterator=True)):
            tm.assert_frame_equal(result,dfs[i])

    def test_invalid_arg(self):
        #GH10369
        class A(object):
            def __init__(self):
                self.read = 0

        tm.assertRaises(ValueError, read_msgpack, path_or_buf=None)
        tm.assertRaises(ValueError, read_msgpack, path_or_buf={})
        tm.assertRaises(ValueError, read_msgpack, path_or_buf=A())


class TestNumpy(TestPackers):

    def test_numpy_scalar_float(self):
        x = np.float32(np.random.rand())
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

    def test_numpy_scalar_complex(self):
        x = np.complex64(np.random.rand() + 1j * np.random.rand())
        x_rec = self.encode_decode(x)
        self.assertTrue(np.allclose(x, x_rec))

    def test_scalar_float(self):
        x = np.random.rand()
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

    def test_scalar_complex(self):
        x = np.random.rand() + 1j * np.random.rand()
        x_rec = self.encode_decode(x)
        self.assertTrue(np.allclose(x, x_rec))

    def test_list_numpy_float(self):
        x = [np.float32(np.random.rand()) for i in range(5)]
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

    def test_list_numpy_float_complex(self):
        if not hasattr(np, 'complex128'):
            raise nose.SkipTest('numpy cant handle complex128')

        x = [np.float32(np.random.rand()) for i in range(5)] + \
            [np.complex128(np.random.rand() + 1j * np.random.rand())
             for i in range(5)]
        x_rec = self.encode_decode(x)
        self.assertTrue(np.allclose(x, x_rec))

    def test_list_float(self):
        x = [np.random.rand() for i in range(5)]
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

    def test_list_float_complex(self):
        x = [np.random.rand() for i in range(5)] + \
            [(np.random.rand() + 1j * np.random.rand()) for i in range(5)]
        x_rec = self.encode_decode(x)
        self.assertTrue(np.allclose(x, x_rec))

    def test_dict_float(self):
        x = {'foo': 1.0, 'bar': 2.0}
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

    def test_dict_complex(self):
        x = {'foo': 1.0 + 1.0j, 'bar': 2.0 + 2.0j}
        x_rec = self.encode_decode(x)
        self.assertEqual(x, x_rec)
        for key in x:
            self.assertEqual(type(x[key]), type(x_rec[key]))

    def test_dict_numpy_float(self):
        x = {'foo': np.float32(1.0), 'bar': np.float32(2.0)}
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

    def test_dict_numpy_complex(self):
        x = {'foo': np.complex128(1.0 + 1.0j),
             'bar': np.complex128(2.0 + 2.0j)}
        x_rec = self.encode_decode(x)
        self.assertEqual(x, x_rec)
        for key in x:
            self.assertEqual(type(x[key]), type(x_rec[key]))

    def test_numpy_array_float(self):

        # run multiple times
        for n in range(10):
            x = np.random.rand(10)
            for dtype in ['float32','float64']:
                x = x.astype(dtype)
                x_rec = self.encode_decode(x)
                tm.assert_almost_equal(x,x_rec)

    def test_numpy_array_complex(self):
        x = (np.random.rand(5) + 1j * np.random.rand(5)).astype(np.complex128)
        x_rec = self.encode_decode(x)
        self.assertTrue(all(map(lambda x, y: x == y, x, x_rec)) and
                        x.dtype == x_rec.dtype)

    def test_list_mixed(self):
        x = [1.0, np.float32(3.5), np.complex128(4.25), u('foo')]
        x_rec = self.encode_decode(x)
        tm.assert_almost_equal(x,x_rec)

class TestBasic(TestPackers):

    def test_timestamp(self):

        for i in [Timestamp(
            '20130101'), Timestamp('20130101', tz='US/Eastern'),
                Timestamp('201301010501')]:
            i_rec = self.encode_decode(i)
            self.assertEqual(i, i_rec)

    def test_datetimes(self):

        # fails under 2.6/win32 (np.datetime64 seems broken)

        if LooseVersion(sys.version) < '2.7':
            raise nose.SkipTest('2.6 with np.datetime64 is broken')

        for i in [datetime.datetime(
            2013, 1, 1), datetime.datetime(2013, 1, 1, 5, 1),
                datetime.date(2013, 1, 1), np.datetime64(datetime.datetime(2013, 1, 5, 2, 15))]:
            i_rec = self.encode_decode(i)
            self.assertEqual(i, i_rec)

    def test_timedeltas(self):

        for i in [datetime.timedelta(days=1),
                  datetime.timedelta(days=1, seconds=10),
                  np.timedelta64(1000000)]:
            i_rec = self.encode_decode(i)
            self.assertEqual(i, i_rec)


class TestIndex(TestPackers):

    def setUp(self):
        super(TestIndex, self).setUp()

        self.d = {
            'string': tm.makeStringIndex(100),
            'date': tm.makeDateIndex(100),
            'int': tm.makeIntIndex(100),
            'float': tm.makeFloatIndex(100),
            'empty': Index([]),
            'tuple': Index(zip(['foo', 'bar', 'baz'], [1, 2, 3])),
            'period': Index(period_range('2012-1-1', freq='M', periods=3)),
            'date2': Index(date_range('2013-01-1', periods=10)),
            'bdate': Index(bdate_range('2013-01-02', periods=10)),
        }

        self.mi = {
            'reg': MultiIndex.from_tuples([('bar', 'one'), ('baz', 'two'), ('foo', 'two'),
                                           ('qux', 'one'), ('qux', 'two')], names=['first', 'second']),
        }

    def test_basic_index(self):

        for s, i in self.d.items():
            i_rec = self.encode_decode(i)
            self.assertTrue(i.equals(i_rec))

        # datetime with no freq (GH5506)
        i = Index([Timestamp('20130101'),Timestamp('20130103')])
        i_rec = self.encode_decode(i)
        self.assertTrue(i.equals(i_rec))

        # datetime with timezone
        i = Index([Timestamp('20130101 9:00:00'),Timestamp('20130103 11:00:00')]).tz_localize('US/Eastern')
        i_rec = self.encode_decode(i)
        self.assertTrue(i.equals(i_rec))

    def test_multi_index(self):

        for s, i in self.mi.items():
            i_rec = self.encode_decode(i)
            self.assertTrue(i.equals(i_rec))

    def test_unicode(self):
        i = tm.makeUnicodeIndex(100)

        # this currently fails
        self.assertRaises(UnicodeEncodeError, self.encode_decode, i)

        #i_rec = self.encode_decode(i)
        #self.assertTrue(i.equals(i_rec))


class TestSeries(TestPackers):

    def setUp(self):
        super(TestSeries, self).setUp()

        self.d = {}

        s = tm.makeStringSeries()
        s.name = 'string'
        self.d['string'] = s

        s = tm.makeObjectSeries()
        s.name = 'object'
        self.d['object'] = s

        s = Series(tslib.iNaT, dtype='M8[ns]', index=range(5))
        self.d['date'] = s

        data = {
            'A': [0., 1., 2., 3., np.nan],
            'B': [0, 1, 0, 1, 0],
            'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
            'D': date_range('1/1/2009', periods=5),
            'E': [0., 1, Timestamp('20100101'), 'foo', 2.],
        }

        self.d['float'] = Series(data['A'])
        self.d['int'] = Series(data['B'])
        self.d['mixed'] = Series(data['E'])

    def test_basic(self):

        # run multiple times here
        for n in range(10):
            for s, i in self.d.items():
                i_rec = self.encode_decode(i)
                assert_series_equal(i, i_rec)


class TestNDFrame(TestPackers):

    def setUp(self):
        super(TestNDFrame, self).setUp()

        data = {
            'A': [0., 1., 2., 3., np.nan],
            'B': [0, 1, 0, 1, 0],
            'C': ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
            'D': date_range('1/1/2009', periods=5),
            'E': [0., 1, Timestamp('20100101'), 'foo', 2.],
        }

        self.frame = {
            'float': DataFrame(dict(A=data['A'], B=Series(data['A']) + 1)),
            'int': DataFrame(dict(A=data['B'], B=Series(data['B']) + 1)),
            'mixed': DataFrame(dict([(k, data[k]) for k in ['A', 'B', 'C', 'D']]))}

        self.panel = {
            'float': Panel(dict(ItemA=self.frame['float'], ItemB=self.frame['float'] + 1))}

    def test_basic_frame(self):

        for s, i in self.frame.items():
            i_rec = self.encode_decode(i)
            assert_frame_equal(i, i_rec)

    def test_basic_panel(self):

        for s, i in self.panel.items():
            i_rec = self.encode_decode(i)
            assert_panel_equal(i, i_rec)

    def test_multi(self):

        i_rec = self.encode_decode(self.frame)
        for k in self.frame.keys():
            assert_frame_equal(self.frame[k], i_rec[k])

        l = tuple(
            [self.frame['float'], self.frame['float'].A, self.frame['float'].B, None])
        l_rec = self.encode_decode(l)
        check_arbitrary(l, l_rec)

        # this is an oddity in that packed lists will be returned as tuples
        l = [self.frame['float'], self.frame['float']
             .A, self.frame['float'].B, None]
        l_rec = self.encode_decode(l)
        self.assertIsInstance(l_rec, tuple)
        check_arbitrary(l, l_rec)

    def test_iterator(self):

        l = [self.frame['float'], self.frame['float']
             .A, self.frame['float'].B, None]

        with ensure_clean(self.path) as path:
            to_msgpack(path, *l)
            for i, packed in enumerate(read_msgpack(path, iterator=True)):
                check_arbitrary(packed, l[i])

    def tests_datetimeindex_freq_issue(self):

        # GH 5947
        # inferring freq on the datetimeindex
        df = DataFrame([1, 2, 3], index=date_range('1/1/2013', '1/3/2013'))
        result = self.encode_decode(df)
        assert_frame_equal(result, df)

        df = DataFrame([1, 2], index=date_range('1/1/2013', '1/2/2013'))
        result = self.encode_decode(df)
        assert_frame_equal(result, df)

    def test_dataframe_duplicate_column_names(self):

        # GH 9618
        expected_1 = DataFrame(columns=['a', 'a'])
        expected_2 = DataFrame(columns=[1]*100)
        expected_2.loc[0] = np.random.randn(100)
        expected_3 = DataFrame(columns=[1, 1])
        expected_3.loc[0] = ['abc', np.nan]

        result_1 = self.encode_decode(expected_1)
        result_2 = self.encode_decode(expected_2)
        result_3 = self.encode_decode(expected_3)

        assert_frame_equal(result_1, expected_1)
        assert_frame_equal(result_2, expected_2)
        assert_frame_equal(result_3, expected_3)


class TestSparse(TestPackers):

    def _check_roundtrip(self, obj, comparator, **kwargs):

        # currently these are not implemetned
        #i_rec = self.encode_decode(obj)
        #comparator(obj, i_rec, **kwargs)
        self.assertRaises(NotImplementedError, self.encode_decode, obj)

    def test_sparse_series(self):

        s = tm.makeStringSeries()
        s[3:5] = np.nan
        ss = s.to_sparse()
        self._check_roundtrip(ss, tm.assert_series_equal,
                              check_series_type=True)

        ss2 = s.to_sparse(kind='integer')
        self._check_roundtrip(ss2, tm.assert_series_equal,
                              check_series_type=True)

        ss3 = s.to_sparse(fill_value=0)
        self._check_roundtrip(ss3, tm.assert_series_equal,
                              check_series_type=True)

    def test_sparse_frame(self):

        s = tm.makeDataFrame()
        s.ix[3:5, 1:3] = np.nan
        s.ix[8:10, -2] = np.nan
        ss = s.to_sparse()

        self._check_roundtrip(ss, tm.assert_frame_equal,
                              check_frame_type=True)

        ss2 = s.to_sparse(kind='integer')
        self._check_roundtrip(ss2, tm.assert_frame_equal,
                              check_frame_type=True)

        ss3 = s.to_sparse(fill_value=0)
        self._check_roundtrip(ss3, tm.assert_frame_equal,
                              check_frame_type=True)

    def test_sparse_panel(self):

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            items = ['x', 'y', 'z']
            p = Panel(dict((i, tm.makeDataFrame().ix[:2, :2]) for i in items))
            sp = p.to_sparse()

            self._check_roundtrip(sp, tm.assert_panel_equal,
                                  check_panel_type=True)

            sp2 = p.to_sparse(kind='integer')
            self._check_roundtrip(sp2, tm.assert_panel_equal,
                                  check_panel_type=True)

            sp3 = p.to_sparse(fill_value=0)
            self._check_roundtrip(sp3, tm.assert_panel_equal,
                                  check_panel_type=True)


class TestCompression(TestPackers):
    """See https://github.com/pydata/pandas/pull/9783
    """

    def setUp(self):
        super(TestCompression, self).setUp()
        data = {
            'A': np.arange(1000, dtype=np.float64),
            'B': np.arange(1000, dtype=np.int32),
            'C': list(100 * 'abcdefghij'),
            'D': date_range(datetime.datetime(2015, 4, 1), periods=1000),
            'E': [datetime.timedelta(days=x) for x in range(1000)],
        }
        self.frame = {
            'float': DataFrame(dict((k, data[k]) for k in ['A', 'A'])),
            'int': DataFrame(dict((k, data[k]) for k in ['B', 'B'])),
            'mixed': DataFrame(data),
        }

    def test_plain(self):
        i_rec = self.encode_decode(self.frame)
        for k in self.frame.keys():
            assert_frame_equal(self.frame[k], i_rec[k])

    def test_compression_zlib(self):
        i_rec = self.encode_decode(self.frame, compress='zlib')
        for k in self.frame.keys():
            assert_frame_equal(self.frame[k], i_rec[k])

    def test_compression_blosc(self):
        try:
            import blosc
        except ImportError:
            raise nose.SkipTest('no blosc')
        i_rec = self.encode_decode(self.frame, compress='blosc')
        for k in self.frame.keys():
            assert_frame_equal(self.frame[k], i_rec[k])


class TestEncoding(TestPackers):
        def setUp(self):
            super(TestEncoding, self).setUp()
            data = {
                'A': [compat.u('\u2019')] * 1000,
                'B': np.arange(1000, dtype=np.int32),
                'C': list(100 * 'abcdefghij'),
                'D': date_range(datetime.datetime(2015, 4, 1), periods=1000),
                'E': [datetime.timedelta(days=x) for x in range(1000)],
                'G': [400] * 1000
            }
            self.frame = {
                'float': DataFrame(dict((k, data[k]) for k in ['A', 'A'])),
                'int': DataFrame(dict((k, data[k]) for k in ['B', 'B'])),
                'mixed': DataFrame(data),
            }
            self.utf_encodings = ['utf8', 'utf16', 'utf32']

        def test_utf(self):
            # GH10581
            for encoding in self.utf_encodings:
                for frame in compat.itervalues(self.frame):
                    result = self.encode_decode(frame, encoding=encoding)
                    assert_frame_equal(result, frame)


class TestMsgpack():
    """
    How to add msgpack tests:

    1. Install pandas version intended to output the msgpack.
TestPackers
    2. Execute "generate_legacy_storage_files.py" to create the msgpack.
    $ python generate_legacy_storage_files.py <output_dir> msgpack

    3. Move the created pickle to "data/legacy_msgpack/<version>" directory.

    NOTE: TestMsgpack can't be a subclass of tm.Testcase to use test generator.
    http://stackoverflow.com/questions/6689537/nose-test-generators-inside-class
    """
    def setUp(self):
        from pandas.io.tests.generate_legacy_storage_files import (
            create_msgpack_data, create_data)
        self.data = create_msgpack_data()
        self.all_data = create_data()
        self.path = u('__%s__.msgpack' % tm.rands(10))
        self.minimum_structure = {'series': ['float', 'int', 'mixed', 'ts', 'mi', 'dup'],
                                  'frame': ['float', 'int', 'mixed', 'mi'],
                                  'panel': ['float'],
                                  'index': ['int', 'date', 'period'],
                                  'mi': ['reg2']}

    def check_min_structure(self, data):
        for typ, v in self.minimum_structure.items():
            assert typ in data, '"{0}" not found in unpacked data'.format(typ)
            for kind in v:
                assert kind in data[typ], '"{0}" not found in data["{1}"]'.format(kind, typ)

    def compare(self, vf, version):
        data = read_msgpack(vf)
        self.check_min_structure(data)
        for typ, dv in data.items():
            assert typ in self.all_data, 'unpacked data contains extra key "{0}"'.format(typ)
            for dt, result in dv.items():
                assert dt in self.all_data[typ], 'data["{0}"] contains extra key "{1}"'.format(typ, dt)
                try:
                    expected = self.data[typ][dt]
                except KeyError:
                    continue

                # use a specific comparator
                # if available
                comparator = getattr(self,"compare_{typ}_{dt}".format(typ=typ,dt=dt), None)
                if comparator is not None:
                    comparator(result, expected, typ, version)
                else:
                    check_arbitrary(result, expected)

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

    def read_msgpacks(self, version):

        pth = tm.get_data_path('legacy_msgpack/{0}'.format(str(version)))
        n = 0
        for f in os.listdir(pth):
            vf = os.path.join(pth, f)
            self.compare(vf, version)
            n += 1
        assert n > 0, 'Msgpack files are not tested'

    def test_msgpack(self):
        msgpack_path = tm.get_data_path('legacy_msgpack')
        n = 0
        for v in os.listdir(msgpack_path):
            pth = os.path.join(msgpack_path, v)
            if os.path.isdir(pth):
                yield self.read_msgpacks, v
            n += 1
        assert n > 0, 'Msgpack files are not tested'


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
