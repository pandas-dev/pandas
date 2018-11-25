from warnings import catch_warnings

import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, Series
from pandas.util import testing as tm


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
class TestMultiIndexIloc(object):

    def test_iloc_getitem_multiindex2(self):
        # TODO(wesm): fix this
        pytest.skip('this test was being suppressed, '
                    'needs to be fixed')

        arr = np.random.randn(3, 3)
        df = DataFrame(arr, columns=[[2, 2, 4], [6, 8, 10]],
                       index=[[4, 4, 8], [8, 10, 12]])

        rs = df.iloc[2]
        xp = Series(arr[2], index=df.columns)
        tm.assert_series_equal(rs, xp)

        rs = df.iloc[:, 2]
        xp = Series(arr[:, 2], index=df.index)
        tm.assert_series_equal(rs, xp)

        rs = df.iloc[2, 2]
        xp = df.values[2, 2]
        assert rs == xp

        # for multiple items
        # GH 5528
        rs = df.iloc[[0, 1]]
        xp = df.xs(4, drop_level=False)
        tm.assert_frame_equal(rs, xp)

        tup = zip(*[['a', 'a', 'b', 'b'], ['x', 'y', 'x', 'y']])
        index = MultiIndex.from_tuples(tup)
        df = DataFrame(np.random.randn(4, 4), index=index)
        rs = df.iloc[[2, 3]]
        xp = df.xs('b', drop_level=False)
        tm.assert_frame_equal(rs, xp)

    def test_iloc_getitem_multiindex(self):
        mi_labels = DataFrame(np.random.randn(4, 3),
                              columns=[['i', 'i', 'j'], ['A', 'A', 'B']],
                              index=[['i', 'i', 'j', 'k'],
                                     ['X', 'X', 'Y', 'Y']])

        mi_int = DataFrame(np.random.randn(3, 3),
                           columns=[[2, 2, 4], [6, 8, 10]],
                           index=[[4, 4, 8], [8, 10, 12]])

        # the first row
        rs = mi_int.iloc[0]
        with catch_warnings(record=True):
            xp = mi_int.ix[4].ix[8]
        tm.assert_series_equal(rs, xp, check_names=False)
        assert rs.name == (4, 8)
        assert xp.name == 8

        # 2nd (last) columns
        rs = mi_int.iloc[:, 2]
        with catch_warnings(record=True):
            xp = mi_int.ix[:, 2]
        tm.assert_series_equal(rs, xp)

        # corner column
        rs = mi_int.iloc[2, 2]
        with catch_warnings(record=True):
            # First level is int - so use .loc rather than .ix (GH 21593)
            xp = mi_int.loc[(8, 12), (4, 10)]
        assert rs == xp

        # this is basically regular indexing
        rs = mi_labels.iloc[2, 2]
        with catch_warnings(record=True):
            xp = mi_labels.ix['j'].ix[:, 'j'].ix[0, 0]
        assert rs == xp

    def test_frame_getitem_setitem_slice(
            self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data
        # getitem
        result = frame.iloc[:4]
        expected = frame[:4]
        tm.assert_frame_equal(result, expected)

        # setitem
        cp = frame.copy()
        cp.iloc[:4] = 0

        assert (cp.values[:4] == 0).all()
        assert (cp.values[4:] != 0).all()

    def test_indexing_ambiguity_bug_1678(self):
        columns = MultiIndex.from_tuples([('Ohio', 'Green'), ('Ohio', 'Red'), (
            'Colorado', 'Green')])
        index = MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1), ('b', 2)
                                        ])

        frame = DataFrame(np.arange(12).reshape((4, 3)), index=index,
                          columns=columns)

        result = frame.iloc[:, 1]
        exp = frame.loc[:, ('Ohio', 'Red')]
        assert isinstance(result, Series)
        tm.assert_series_equal(result, exp)

    def test_iloc_mi(self):
        # GH 13797
        # Test if iloc can handle integer locations in MultiIndexed DataFrame

        data = [['str00', 'str01'], ['str10', 'str11'], ['str20', 'srt21'],
                ['str30', 'str31'], ['str40', 'str41']]

        mi = MultiIndex.from_tuples(
            [('CC', 'A'), ('CC', 'B'), ('CC', 'B'), ('BB', 'a'), ('BB', 'b')])

        expected = DataFrame(data)
        df_mi = DataFrame(data, index=mi)

        result = DataFrame([[df_mi.iloc[r, c] for c in range(2)]
                            for r in range(5)])

        tm.assert_frame_equal(result, expected)
