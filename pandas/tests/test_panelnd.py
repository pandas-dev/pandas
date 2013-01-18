from datetime import datetime
import os
import operator
import unittest
import nose

import numpy as np

from pandas.core import panelnd
from pandas.core.panel import Panel
import pandas.core.common as com
from pandas.util import py3compat

from pandas.util.testing import (assert_panel_equal,
                                 assert_panel4d_equal,
                                 assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal)
import pandas.util.testing as tm


class TestPanelnd(unittest.TestCase):

    def setUp(self):
        pass

    def test_4d_construction(self):

        # create a 4D
        Panel4D = panelnd.create_nd_panel_factory(
            klass_name='Panel4D',
            axis_orders=['labels', 'items', 'major_axis', 'minor_axis'],
            axis_slices={'items': 'items', 'major_axis': 'major_axis',
                         'minor_axis': 'minor_axis'},
            slicer=Panel,
            axis_aliases={'major': 'major_axis', 'minor': 'minor_axis'},
            stat_axis=2)

        p4d = Panel4D(dict(L1=tm.makePanel(), L2=tm.makePanel()))

    def test_4d_construction_alt(self):

        # create a 4D
        Panel4D = panelnd.create_nd_panel_factory(
            klass_name='Panel4D',
            axis_orders=['labels', 'items', 'major_axis', 'minor_axis'],
            axis_slices={'items': 'items', 'major_axis': 'major_axis',
                         'minor_axis': 'minor_axis'},
            slicer='Panel',
            axis_aliases={'major': 'major_axis', 'minor': 'minor_axis'},
            stat_axis=2)

        p4d = Panel4D(dict(L1=tm.makePanel(), L2=tm.makePanel()))

    def test_4d_construction_error(self):

        # create a 4D
        self.assertRaises(Exception,
                          panelnd.create_nd_panel_factory,
                          klass_name='Panel4D',
                          axis_orders=['labels', 'items', 'major_axis',
                                       'minor_axis'],
                          axis_slices={'items': 'items',
                                       'major_axis': 'major_axis',
                                       'minor_axis': 'minor_axis'},
                          slicer='foo',
                          axis_aliases={'major': 'major_axis',
                                        'minor': 'minor_axis'},
                          stat_axis=2)

    def test_5d_construction(self):

        # create a 4D
        Panel4D = panelnd.create_nd_panel_factory(
            klass_name='Panel4D',
            axis_orders=['labels1', 'items', 'major_axis', 'minor_axis'],
            axis_slices={'items': 'items', 'major_axis': 'major_axis',
                         'minor_axis': 'minor_axis'},
            slicer=Panel,
            axis_aliases={'major': 'major_axis', 'minor': 'minor_axis'},
            stat_axis=2)

        p4d = Panel4D(dict(L1=tm.makePanel(), L2=tm.makePanel()))

        # create a 5D
        Panel5D = panelnd.create_nd_panel_factory(
            klass_name='Panel5D',
            axis_orders=['cool1', 'labels1', 'items', 'major_axis',
                         'minor_axis'],
            axis_slices={'labels1': 'labels1', 'items': 'items',
                         'major_axis': 'major_axis',
                         'minor_axis': 'minor_axis'},
            slicer=Panel4D,
            axis_aliases={'major': 'major_axis', 'minor': 'minor_axis'},
            stat_axis=2)

        p5d = Panel5D(dict(C1=p4d))

        # slice back to 4d
        results = p5d.ix['C1', :, :, 0:3, :]
        expected = p4d.ix[:, :, 0:3, :]
        assert_panel_equal(results['L1'], expected['L1'])

        # test a transpose
        # results  = p5d.transpose(1,2,3,4,0)
        # expected =

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
