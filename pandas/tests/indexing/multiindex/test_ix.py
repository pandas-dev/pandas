from warnings import catch_warnings, simplefilter

import pytest

from pandas.compat import lrange


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
class TestMultiIndexIx(object):

    def test_frame_setitem_ix(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data
        frame.loc[('bar', 'two'), 'B'] = 5
        assert frame.loc[('bar', 'two'), 'B'] == 5

        # with integer labels
        df = frame.copy()
        df.columns = lrange(3)
        df.loc[('bar', 'two'), 1] = 7
        assert df.loc[('bar', 'two'), 1] == 7

        with catch_warnings(record=True):
            simplefilter("ignore", DeprecationWarning)
            df = frame.copy()
            df.columns = lrange(3)
            df.ix[('bar', 'two'), 1] = 7
        assert df.loc[('bar', 'two'), 1] == 7
