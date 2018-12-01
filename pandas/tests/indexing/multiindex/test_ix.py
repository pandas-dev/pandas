import pytest

from pandas.compat import lrange


@pytest.mark.filterwarnings("ignore:\\n.ix:DeprecationWarning")
def test_frame_setitem_ix(multiindex_dataframe_random_data):
    frame = multiindex_dataframe_random_data
    frame.loc[('bar', 'two'), 'B'] = 5
    assert frame.loc[('bar', 'two'), 'B'] == 5

    # with integer labels
    df = frame.copy()
    df.columns = lrange(3)
    df.loc[('bar', 'two'), 1] = 7
    assert df.loc[('bar', 'two'), 1] == 7

    df = frame.copy()
    df.columns = lrange(3)
    df.ix[('bar', 'two'), 1] = 7
    assert df.loc[('bar', 'two'), 1] == 7
