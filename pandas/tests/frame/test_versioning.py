# pandas/tests/frame/test_versioning.py
import pandas as pd
import pytest


def test_snapshot_and_restore_returns_dataframe():
    df = pd.DataFrame({"x": [1, 2, 3]})
    sid = df.snapshot("t1")
    assert sid in df.list_snapshots()
    df.loc[0, "x"] = 99
    restored = df.restore(sid)
    assert list(restored["x"]) == [1, 2, 3]


def test_restore_inplace_mutates_dataframe():
    df = pd.DataFrame({"x": [1, 2, 3]})
    sid = df.snapshot("t2")
    df.loc[1, "x"] = 999
    df.restore(sid, inplace=True)
    assert list(df["x"]) == [1, 2, 3]


def test_drop_and_clear_behaviour():
    df = pd.DataFrame({"a": [1, 2]})
    sid1 = df.snapshot("s1")
    sid2 = df.snapshot("s2")
    assert set(df.list_snapshots()) == {sid1, sid2}
    df.drop_snapshot(sid1)
    assert sid1 not in df.list_snapshots()
    df.clear_snapshots()
    assert df.list_snapshots() == []


def test_snapshot_on_empty_dataframe():
    df = pd.DataFrame()
    sid = df.snapshot()
    df.loc[0, "a"] = 1
    restored = df.restore(sid)
    assert restored.empty


def test_copy_does_not_inherit_snapshots():
    df = pd.DataFrame({"a": [1, 2, 3]})
    sid = df.snapshot("orig")
    df2 = df.copy()
    # design decision: copies do not copy snapshots
    assert df2.list_snapshots() == []
    assert sid in df.list_snapshots()


def test_missing_snapshot_raises():
    df = pd.DataFrame({"x": [1]})
    with pytest.raises(KeyError):
        df.restore("no-such-snapshot")
