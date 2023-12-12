import pandas as pd
import pytest
# to test and assert that Timestamp.fromisoformat should return a timezone aware timestamp

def test_timestamp_fromisoformat_timezone_handling():
    # Test case for ISO 8601 string with Zulu (UTC) timezone
    timestr_utc = '2023-11-05T08:30:00Z'
    ts_utc = pd.Timestamp.fromisoformat(timestr_utc)
    assert ts_utc == pd.to_datetime(timestr_utc), "UTC timezone not handled correctly"

    # Test case for ISO 8601 string with the offset timezone
    timestr_offset = '2023-11-05T08:30:00+0000'
    ts_offset = pd.Timestamp.fromisoformat(timestr_offset)
    assert ts_offset == pd.to_datetime(timestr_offset), "Offset timezone not handled correctly"

    # Test case for a roundtrip conversion
    ts_now = pd.Timestamp.utcnow()
    assert ts_now == pd.Timestamp.fromisoformat(ts_now.isoformat()), "Roundtrip conversion loses timezone information"
