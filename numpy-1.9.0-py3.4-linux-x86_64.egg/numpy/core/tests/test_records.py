from __future__ import division, absolute_import, print_function

from os import path
import numpy as np
from numpy.testing import *
from numpy.compat import asbytes, asunicode

import warnings
import collections
import pickle


class TestFromrecords(TestCase):
    def test_fromrecords(self):
        r = np.rec.fromrecords([[456, 'dbe', 1.2], [2, 'de', 1.3]],
                            names='col1,col2,col3')
        assert_equal(r[0].item(), (456, 'dbe', 1.2))

    def test_method_array(self):
        r = np.rec.array(asbytes('abcdefg') * 100, formats='i2,a3,i4', shape=3, byteorder='big')
        assert_equal(r[1].item(), (25444, asbytes('efg'), 1633837924))

    def test_method_array2(self):
        r = np.rec.array([(1, 11, 'a'), (2, 22, 'b'), (3, 33, 'c'), (4, 44, 'd'), (5, 55, 'ex'),
                     (6, 66, 'f'), (7, 77, 'g')], formats='u1,f4,a1')
        assert_equal(r[1].item(), (2, 22.0, asbytes('b')))

    def test_recarray_slices(self):
        r = np.rec.array([(1, 11, 'a'), (2, 22, 'b'), (3, 33, 'c'), (4, 44, 'd'), (5, 55, 'ex'),
                     (6, 66, 'f'), (7, 77, 'g')], formats='u1,f4,a1')
        assert_equal(r[1::2][1].item(), (4, 44.0, asbytes('d')))

    def test_recarray_fromarrays(self):
        x1 = np.array([1, 2, 3, 4])
        x2 = np.array(['a', 'dd', 'xyz', '12'])
        x3 = np.array([1.1, 2, 3, 4])
        r = np.rec.fromarrays([x1, x2, x3], names='a,b,c')
        assert_equal(r[1].item(), (2, 'dd', 2.0))
        x1[1] = 34
        assert_equal(r.a, np.array([1, 2, 3, 4]))

    def test_recarray_fromfile(self):
        data_dir = path.join(path.dirname(__file__), 'data')
        filename = path.join(data_dir, 'recarray_from_file.fits')
        fd = open(filename, 'rb')
        fd.seek(2880 * 2)
        r = np.rec.fromfile(fd, formats='f8,i4,a5', shape=3, byteorder='big')
        fd.seek(2880 * 2)
        r = np.rec.array(fd, formats='f8,i4,a5', shape=3, byteorder='big')
        fd.close()

    def test_recarray_from_obj(self):
        count = 10
        a = np.zeros(count, dtype='O')
        b = np.zeros(count, dtype='f8')
        c = np.zeros(count, dtype='f8')
        for i in range(len(a)):
            a[i] = list(range(1, 10))

        mine = np.rec.fromarrays([a, b, c], names='date,data1,data2')
        for i in range(len(a)):
            assert_((mine.date[i] == list(range(1, 10))))
            assert_((mine.data1[i] == 0.0))
            assert_((mine.data2[i] == 0.0))

    def test_recarray_from_repr(self):
        x = np.rec.array([ (1, 2)], dtype=[('a', np.int8), ('b', np.int8)])
        y = eval("np." + repr(x))
        assert_(isinstance(y, np.recarray))
        assert_equal(y, x)

    def test_recarray_from_names(self):
        ra = np.rec.array([
            (1, 'abc', 3.7000002861022949, 0),
            (2, 'xy', 6.6999998092651367, 1),
            (0, ' ', 0.40000000596046448, 0)],
                       names='c1, c2, c3, c4')
        pa = np.rec.fromrecords([
            (1, 'abc', 3.7000002861022949, 0),
            (2, 'xy', 6.6999998092651367, 1),
            (0, ' ', 0.40000000596046448, 0)],
                       names='c1, c2, c3, c4')
        assert_(ra.dtype == pa.dtype)
        assert_(ra.shape == pa.shape)
        for k in range(len(ra)):
            assert_(ra[k].item() == pa[k].item())

    def test_recarray_conflict_fields(self):
        ra = np.rec.array([(1, 'abc', 2.3), (2, 'xyz', 4.2),
                        (3, 'wrs', 1.3)],
                       names='field, shape, mean')
        ra.mean = [1.1, 2.2, 3.3]
        assert_array_almost_equal(ra['mean'], [1.1, 2.2, 3.3])
        assert_(type(ra.mean) is type(ra.var))
        ra.shape = (1, 3)
        assert_(ra.shape == (1, 3))
        ra.shape = ['A', 'B', 'C']
        assert_array_equal(ra['shape'], [['A', 'B', 'C']])
        ra.field = 5
        assert_array_equal(ra['field'], [[5, 5, 5]])
        assert_(isinstance(ra.field, collections.Callable))

    def test_fromrecords_with_explicit_dtype(self):
        a = np.rec.fromrecords([(1, 'a'), (2, 'bbb')],
                                dtype=[('a', int), ('b', np.object)])
        assert_equal(a.a, [1, 2])
        assert_equal(a[0].a, 1)
        assert_equal(a.b, ['a', 'bbb'])
        assert_equal(a[-1].b, 'bbb')
        #
        ndtype = np.dtype([('a', int), ('b', np.object)])
        a = np.rec.fromrecords([(1, 'a'), (2, 'bbb')], dtype=ndtype)
        assert_equal(a.a, [1, 2])
        assert_equal(a[0].a, 1)
        assert_equal(a.b, ['a', 'bbb'])
        assert_equal(a[-1].b, 'bbb')


class TestRecord(TestCase):
    def setUp(self):
        self.data = np.rec.fromrecords([(1, 2, 3), (4, 5, 6)],
                            dtype=[("col1", "<i4"),
                                   ("col2", "<i4"),
                                   ("col3", "<i4")])

    def test_assignment1(self):
        a = self.data
        assert_equal(a.col1[0], 1)
        a[0].col1 = 0
        assert_equal(a.col1[0], 0)

    def test_assignment2(self):
        a = self.data
        assert_equal(a.col1[0], 1)
        a.col1[0] = 0
        assert_equal(a.col1[0], 0)

    def test_invalid_assignment(self):
        a = self.data
        def assign_invalid_column(x):
            x[0].col5 = 1
        self.assertRaises(AttributeError, assign_invalid_column, a)

    def test_out_of_order_fields(self):
        """Ticket #1431."""
        x = self.data[['col1', 'col2']]
        y = self.data[['col2', 'col1']]
        assert_equal(x[0][0], y[0][1])

    def test_pickle_1(self):
        # Issue #1529
        a = np.array([(1, [])], dtype=[('a', np.int32), ('b', np.int32, 0)])
        assert_equal(a, pickle.loads(pickle.dumps(a)))
        assert_equal(a[0], pickle.loads(pickle.dumps(a[0])))

    def test_pickle_2(self):
        a = self.data
        assert_equal(a, pickle.loads(pickle.dumps(a)))
        assert_equal(a[0], pickle.loads(pickle.dumps(a[0])))


def test_find_duplicate():
    l1 = [1, 2, 3, 4, 5, 6]
    assert_(np.rec.find_duplicate(l1) == [])

    l2 = [1, 2, 1, 4, 5, 6]
    assert_(np.rec.find_duplicate(l2) == [1])

    l3 = [1, 2, 1, 4, 1, 6, 2, 3]
    assert_(np.rec.find_duplicate(l3) == [1, 2])

    l3 = [2, 2, 1, 4, 1, 6, 2, 3]
    assert_(np.rec.find_duplicate(l3) == [2, 1])

if __name__ == "__main__":
    run_module_suite()
