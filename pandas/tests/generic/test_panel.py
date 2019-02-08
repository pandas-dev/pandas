# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from warnings import catch_warnings, simplefilter

from pandas import Panel
from pandas.util.testing import assert_panel_equal

from .test_generic import Generic


class TestPanel(Generic):
    _typ = Panel
    _comparator = lambda self, x, y: assert_panel_equal(x, y, by_blocks=True)


# run all the tests, but wrap each in a warning catcher
for t in ['test_rename', 'test_get_numeric_data',
          'test_get_default', 'test_nonzero',
          'test_downcast', 'test_constructor_compound_dtypes',
          'test_head_tail',
          'test_size_compat', 'test_split_compat',
          'test_unexpected_keyword',
          'test_stat_unexpected_keyword', 'test_api_compat',
          'test_stat_non_defaults_args',
          'test_truncate_out_of_bounds',
          'test_metadata_propagation', 'test_copy_and_deepcopy',
          'test_pct_change', 'test_sample']:

    def f():
        def tester(self):
            f = getattr(super(TestPanel, self), t)
            with catch_warnings(record=True):
                simplefilter("ignore", FutureWarning)
                f()
        return tester

    setattr(TestPanel, t, f())
