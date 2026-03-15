import pandas as pd
import datetime as dt

def test_loc_setitem_expansion_preserves_tz_aware_dtype():
    # GH#55423: Test for tz-aware datetime expansion
    df = pd.DataFrame([{'id': 1}, {'id': 2}, {'id': 3}])
    _time = dt.datetime.fromtimestamp(1695887042, tz=dt.timezone.utc)
    df.loc[df.id >= 2, 'time'] = _time
    assert str(df['time'].dtype).startswith('datetime64[ns, UTC]')
