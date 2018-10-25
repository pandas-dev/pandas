"""
Tests for pd.Period behavior that depend on the implementations of Index/Series
"""
import pytest

import pandas.util.testing as tm
from pandas import Index, Period, Series, Timestamp
from pandas.core.indexes.period import pnow

boxes = [lambda x: x, lambda x: Series([x]), lambda x: Index([x])]
ids = ['identity', 'Series', 'Index']


@pytest.mark.parametrize('lbox', boxes, ids=ids)
@pytest.mark.parametrize('rbox', boxes, ids=ids)
def test_add_timestamp_raises(self, rbox, lbox):
    # GH#17983
    ts = Timestamp('2017')
    per = Period('2017', freq='M')

    # We may get a different message depending on which class raises
    # the error.
    msg = (r"cannot add|unsupported operand|"
           r"can only operate on a|incompatible type|"
           r"ufunc add cannot use operands")
    with tm.assert_raises_regex(TypeError, msg):
        lbox(ts) + rbox(per)

    with tm.assert_raises_regex(TypeError, msg):
        lbox(per) + rbox(ts)

    with tm.assert_raises_regex(TypeError, msg):
        lbox(per) + rbox(per)


def test_pnow(self):
    # deprecation, xref GH#13790
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        pnow('D')
