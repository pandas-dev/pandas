# pylint: disable-msg=E1101,W0612

import nose  # noqa
from numpy import nan
import pandas as pd

from pandas import DataFrame, bdate_range, Panel
from pandas.core.index import Index
import pandas.util.testing as tm
from pandas.sparse.api import SparseSeries, SparsePanel
import pandas.tests.test_panel as test_panel


def panel_data1():
    index = bdate_range('1/1/2011', periods=8)

    return DataFrame({
        'A': [nan, nan, nan, 0, 1, 2, 3, 4],
        'B': [0, 1, 2, 3, 4, nan, nan, nan],
        'C': [0, 1, 2, nan, nan, nan, 3, 4],
        'D': [nan, 0, 1, nan, 2, 3, 4, nan]
    }, index=index)


def panel_data2():
    index = bdate_range('1/1/2011', periods=9)

    return DataFrame({
        'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5],
        'B': [0, 1, 2, 3, 4, 5, nan, nan, nan],
        'C': [0, 1, 2, nan, nan, nan, 3, 4, 5],
        'D': [nan, 0, 1, nan, 2, 3, 4, 5, nan]
    }, index=index)


def panel_data3():
    index = bdate_range('1/1/2011', periods=10).shift(-2)

    return DataFrame({
        'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
        'B': [0, 1, 2, 3, 4, 5, 6, nan, nan, nan],
        'C': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
        'D': [nan, 0, 1, nan, 2, 3, 4, 5, 6, nan]
    }, index=index)


class TestSparsePanel(tm.TestCase, test_panel.SafeForLongAndSparse,
                      test_panel.SafeForSparse):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.data_dict = {
            'ItemA': panel_data1(),
            'ItemB': panel_data2(),
            'ItemC': panel_data3(),
            'ItemD': panel_data1(),
        }
        with tm.assert_produces_warning(FutureWarning):
            self.panel = SparsePanel(self.data_dict)

    @staticmethod
    def _test_op(panel, op):
        # arithmetic tests
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = op(panel, 1)
        tm.assert_sp_frame_equal(result['ItemA'], op(panel['ItemA'], 1))

    def test_constructor(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            self.assertRaises(ValueError, SparsePanel, self.data_dict,
                              items=['Item0', 'ItemA', 'ItemB'])
            with tm.assertRaisesRegexp(TypeError,
                                       "input must be a dict, a 'list' was "
                                       "passed"):
                SparsePanel(['a', 'b', 'c'])

    # deprecation GH11157
    def test_deprecation(self):
        with tm.assert_produces_warning(FutureWarning):
            SparsePanel()

    # GH 9272
    def test_constructor_empty(self):
        with tm.assert_produces_warning(FutureWarning):
            sp = SparsePanel()
        self.assertEqual(len(sp.items), 0)
        self.assertEqual(len(sp.major_axis), 0)
        self.assertEqual(len(sp.minor_axis), 0)

    def test_from_dict(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            fd = SparsePanel.from_dict(self.data_dict)
        tm.assert_sp_panel_equal(fd, self.panel)

    def test_pickle(self):
        def _test_roundtrip(panel):
            result = self.round_trip_pickle(panel)
            tm.assertIsInstance(result.items, Index)
            tm.assertIsInstance(result.major_axis, Index)
            tm.assertIsInstance(result.minor_axis, Index)
            tm.assert_sp_panel_equal(panel, result)

        _test_roundtrip(self.panel)

    def test_dense_to_sparse(self):
        wp = Panel.from_dict(self.data_dict)
        dwp = wp.to_sparse()
        tm.assertIsInstance(dwp['ItemA']['A'], SparseSeries)

    def test_to_dense(self):
        dwp = self.panel.to_dense()
        dwp2 = Panel.from_dict(self.data_dict)
        tm.assert_panel_equal(dwp, dwp2)

    def test_to_frame(self):

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):

            def _compare_with_dense(panel):
                slp = panel.to_frame()
                dlp = panel.to_dense().to_frame()

                self.assert_numpy_array_equal(slp.values, dlp.values)
                self.assertTrue(slp.index.equals(dlp.index))

            _compare_with_dense(self.panel)
            _compare_with_dense(self.panel.reindex(items=['ItemA']))

            with tm.assert_produces_warning(FutureWarning):
                zero_panel = SparsePanel(self.data_dict, default_fill_value=0)
            self.assertRaises(Exception, zero_panel.to_frame)

            self.assertRaises(Exception, self.panel.to_frame,
                              filter_observations=False)

    def test_long_to_wide_sparse(self):
        pass

    def test_values(self):
        pass

    def test_setitem(self):
        self.panel['ItemE'] = self.panel['ItemC']
        self.panel['ItemF'] = self.panel['ItemC'].to_dense()

        tm.assert_sp_frame_equal(self.panel['ItemE'], self.panel['ItemC'])
        tm.assert_sp_frame_equal(self.panel['ItemF'], self.panel['ItemC'])

        expected = pd.Index(['ItemA', 'ItemB', 'ItemC',
                             'ItemD', 'ItemE', 'ItemF'])
        tm.assert_index_equal(self.panel.items, expected)

        self.assertRaises(Exception, self.panel.__setitem__, 'item6', 1)

    def test_set_value(self):
        def _check_loc(item, major, minor, val=1.5):
            res = self.panel.set_value(item, major, minor, val)
            self.assertIsNot(res, self.panel)
            self.assertEqual(res.get_value(item, major, minor), val)

        _check_loc('ItemA', self.panel.major_axis[4], self.panel.minor_axis[3])
        _check_loc('ItemF', self.panel.major_axis[4], self.panel.minor_axis[3])
        _check_loc('ItemF', 'foo', self.panel.minor_axis[3])
        _check_loc('ItemE', 'foo', 'bar')

    def test_delitem_pop(self):
        del self.panel['ItemB']
        tm.assert_index_equal(self.panel.items,
                              pd.Index(['ItemA', 'ItemC', 'ItemD']))
        crackle = self.panel['ItemC']
        pop = self.panel.pop('ItemC')
        self.assertIs(pop, crackle)
        tm.assert_almost_equal(self.panel.items, pd.Index(['ItemA', 'ItemD']))

        self.assertRaises(KeyError, self.panel.__delitem__, 'ItemC')

    def test_copy(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            cop = self.panel.copy()
        tm.assert_sp_panel_equal(cop, self.panel)

    def test_reindex(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):

            def _compare_with_dense(swp, items, major, minor):
                swp_re = swp.reindex(items=items, major=major, minor=minor)
                dwp_re = swp.to_dense().reindex(items=items, major=major,
                                                minor=minor)
                tm.assert_panel_equal(swp_re.to_dense(), dwp_re)

            _compare_with_dense(self.panel, self.panel.items[:2],
                                self.panel.major_axis[::2],
                                self.panel.minor_axis[::2])
            _compare_with_dense(self.panel, None, self.panel.major_axis[::2],
                                self.panel.minor_axis[::2])

            self.assertRaises(ValueError, self.panel.reindex)

            # TODO: do something about this later...
            self.assertRaises(Exception, self.panel.reindex,
                              items=['item0', 'ItemA', 'ItemB'])

            # test copying
            cp = self.panel.reindex(self.panel.major_axis, copy=True)
            cp['ItemA']['E'] = cp['ItemA']['A']
            self.assertNotIn('E', self.panel['ItemA'])

    def test_operators(self):
        def _check_ops(panel):
            def _dense_comp(op):
                with tm.assert_produces_warning(FutureWarning,
                                                check_stacklevel=False):
                    dense = panel.to_dense()
                    sparse_result = op(panel)
                    dense_result = op(dense)
                    tm.assert_panel_equal(sparse_result.to_dense(),
                                          dense_result)

            def _mixed_comp(op):
                with tm.assert_produces_warning(FutureWarning,
                                                check_stacklevel=False):
                    result = op(panel, panel.to_dense())
                    expected = op(panel.to_dense(), panel.to_dense())
                    tm.assert_panel_equal(result, expected)

            op1 = lambda x: x + 2

            _dense_comp(op1)
            op2 = lambda x: x.add(x.reindex(major=x.major_axis[::2]))
            _dense_comp(op2)
            op3 = lambda x: x.subtract(x.mean(0), axis=0)
            _dense_comp(op3)
            op4 = lambda x: x.subtract(x.mean(1), axis=1)
            _dense_comp(op4)
            op5 = lambda x: x.subtract(x.mean(2), axis=2)
            _dense_comp(op5)

            _mixed_comp(Panel.multiply)
            _mixed_comp(Panel.subtract)

            # TODO: this case not yet supported!
            # op6 = lambda x: x.add(x.to_frame())
            # _dense_comp(op6)

        _check_ops(self.panel)

    def test_major_xs(self):
        def _dense_comp(sparse):
            dense = sparse.to_dense()

            for idx in sparse.major_axis:
                dslice = dense.major_xs(idx)
                sslice = sparse.major_xs(idx)
                tm.assert_frame_equal(dslice, sslice)

        _dense_comp(self.panel)

    def test_minor_xs(self):
        def _dense_comp(sparse):
            dense = sparse.to_dense()

            for idx in sparse.minor_axis:
                dslice = dense.minor_xs(idx)
                sslice = sparse.minor_xs(idx).to_dense()
                tm.assert_frame_equal(dslice, sslice)

        _dense_comp(self.panel)


if __name__ == '__main__':
    import nose  # noqa
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
