import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, Panel, Series
from pandas.util import testing as tm


@pytest.mark.filterwarnings('ignore:\\nPanel:FutureWarning')
class TestMultiIndexPanel(object):

    def test_iloc_getitem_panel_multiindex(self):

        # GH 7199
        # Panel with multi-index
        multi_index = MultiIndex.from_tuples([('ONE', 'one'),
                                              ('TWO', 'two'),
                                              ('THREE', 'three')],
                                             names=['UPPER', 'lower'])

        simple_index = [x[0] for x in multi_index]
        wd1 = Panel(items=['First', 'Second'],
                    major_axis=['a', 'b', 'c', 'd'],
                    minor_axis=multi_index)

        wd2 = Panel(items=['First', 'Second'],
                    major_axis=['a', 'b', 'c', 'd'],
                    minor_axis=simple_index)

        expected1 = wd1['First'].iloc[[True, True, True, False], [0, 2]]
        result1 = wd1.iloc[0, [True, True, True, False], [0, 2]]  # WRONG
        tm.assert_frame_equal(result1, expected1)

        expected2 = wd2['First'].iloc[[True, True, True, False], [0, 2]]
        result2 = wd2.iloc[0, [True, True, True, False], [0, 2]]
        tm.assert_frame_equal(result2, expected2)

        expected1 = DataFrame(index=['a'], columns=multi_index,
                              dtype='float64')
        result1 = wd1.iloc[0, [0], [0, 1, 2]]
        tm.assert_frame_equal(result1, expected1)

        expected2 = DataFrame(index=['a'], columns=simple_index,
                              dtype='float64')
        result2 = wd2.iloc[0, [0], [0, 1, 2]]
        tm.assert_frame_equal(result2, expected2)

        # GH 7516
        mi = MultiIndex.from_tuples([(0, 'x'), (1, 'y'), (2, 'z')])
        p = Panel(np.arange(3 * 3 * 3, dtype='int64').reshape(3, 3, 3),
                  items=['a', 'b', 'c'], major_axis=mi,
                  minor_axis=['u', 'v', 'w'])
        result = p.iloc[:, 1, 0]
        expected = Series([3, 12, 21], index=['a', 'b', 'c'], name='u')
        tm.assert_series_equal(result, expected)

        result = p.loc[:, (1, 'y'), 'u']
        tm.assert_series_equal(result, expected)
