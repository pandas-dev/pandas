from warnings import catch_warnings

import numpy as np
import pytest

from pandas import Panel, date_range
from pandas.util import testing as tm


@pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
class TestPanel(object):

    def test_iloc_getitem_panel(self):

        with catch_warnings(record=True):
            # GH 7189
            p = Panel(np.arange(4 * 3 * 2).reshape(4, 3, 2),
                      items=['A', 'B', 'C', 'D'],
                      major_axis=['a', 'b', 'c'],
                      minor_axis=['one', 'two'])

            result = p.iloc[1]
            expected = p.loc['B']
            tm.assert_frame_equal(result, expected)

            result = p.iloc[1, 1]
            expected = p.loc['B', 'b']
            tm.assert_series_equal(result, expected)

            result = p.iloc[1, 1, 1]
            expected = p.loc['B', 'b', 'two']
            assert result == expected

            # combined
            result = p.iloc[0, [True, True], [0, 1]]
            expected = p.loc['A', ['a', 'b'], ['one', 'two']]
            tm.assert_frame_equal(result, expected)

            # out-of-bounds exception
            with pytest.raises(IndexError):
                p.iloc[tuple([10, 5])]

            with pytest.raises(IndexError):
                p.iloc[0, [True, True], [0, 1, 2]]

            # trying to use a label
            with pytest.raises(ValueError):
                p.iloc[tuple(['j', 'D'])]

            # GH
            p = Panel(
                np.random.rand(4, 3, 2), items=['A', 'B', 'C', 'D'],
                major_axis=['U', 'V', 'W'], minor_axis=['X', 'Y'])
            expected = p['A']

            result = p.iloc[0, :, :]
            tm.assert_frame_equal(result, expected)

            result = p.iloc[0, [True, True, True], :]
            tm.assert_frame_equal(result, expected)

            result = p.iloc[0, [True, True, True], [0, 1]]
            tm.assert_frame_equal(result, expected)

            with pytest.raises(IndexError):
                p.iloc[0, [True, True, True], [0, 1, 2]]

            with pytest.raises(IndexError):
                p.iloc[0, [True, True, True], [2]]

    def test_iloc_panel_issue(self):

        with catch_warnings(record=True):
            # see gh-3617
            p = Panel(np.random.randn(4, 4, 4))

            assert p.iloc[:3, :3, :3].shape == (3, 3, 3)
            assert p.iloc[1, :3, :3].shape == (3, 3)
            assert p.iloc[:3, 1, :3].shape == (3, 3)
            assert p.iloc[:3, :3, 1].shape == (3, 3)
            assert p.iloc[1, 1, :3].shape == (3, )
            assert p.iloc[1, :3, 1].shape == (3, )
            assert p.iloc[:3, 1, 1].shape == (3, )

    @pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
    def test_panel_getitem(self):

        with catch_warnings(record=True):
            # with an object-like
            # GH 9140
            class TestObject(object):

                def __str__(self):
                    return "TestObject"

            obj = TestObject()

            p = Panel(np.random.randn(1, 5, 4), items=[obj],
                      major_axis=date_range('1/1/2000', periods=5),
                      minor_axis=['A', 'B', 'C', 'D'])

            expected = p.iloc[0]
            result = p[obj]
            tm.assert_frame_equal(result, expected)
