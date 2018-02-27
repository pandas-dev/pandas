from warnings import catch_warnings

from pandas.tests.io.pytables.common import (_check_roundtrip,
                                             ensure_clean_store)
import pandas.util.testing as tm
from pandas.util.testing import assert_panel_equal


def test_wide():
    with catch_warnings(record=True):
        wp = tm.makePanel()
        _check_roundtrip(wp, assert_panel_equal)


def test_wide_table_dups():
    with ensure_clean_store() as store:
        with catch_warnings(record=True):
            wp = tm.makePanel()
            store.put('panel', wp, format='table')
            store.put('panel', wp, format='table', append=True)

            recons = store['panel']

            assert_panel_equal(recons, wp)


def test_long():
    def _check(left, right):
        assert_panel_equal(left.to_panel(), right.to_panel())

    with catch_warnings(record=True):
        wp = tm.makePanel()
        _check_roundtrip(wp.to_frame(), _check)
