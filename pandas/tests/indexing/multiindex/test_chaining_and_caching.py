import numpy as np

from pandas._libs import index as libindex

from pandas import (
    DataFrame,
    MultiIndex,
    RangeIndex,
    Series,
)
import pandas._testing as tm


def test_detect_chained_assignment():
    # Inplace ops, originally from:
    # https://stackoverflow.com/questions/20508968/series-fillna-in-a-multiindex-dataframe-does-not-fill-is-this-a-bug
    a = [12, 23]
    b = [123, None]
    c = [1234, 2345]
    d = [12345, 23456]
    tuples = [("eyes", "left"), ("eyes", "right"), ("ears", "left"), ("ears", "right")]
    events = {
        ("eyes", "left"): a,
        ("eyes", "right"): b,
        ("ears", "left"): c,
        ("ears", "right"): d,
    }
    multiind = MultiIndex.from_tuples(tuples, names=["part", "side"])
    zed = DataFrame(events, index=["a", "b"], columns=multiind)

    with tm.raises_chained_assignment_error():
        zed["eyes"]["right"].fillna(value=555, inplace=True)


def test_cache_updating():
    # 5216
    # make sure that we don't try to set a dead cache
    a = np.random.default_rng(2).random((10, 3))
    df = DataFrame(a, columns=["x", "y", "z"])
    df_original = df.copy()
    tuples = [(i, j) for i in range(5) for j in range(2)]
    index = MultiIndex.from_tuples(tuples)
    df.index = index

    # setting via chained assignment
    # but actually works, since everything is a view

    with tm.raises_chained_assignment_error():
        df.loc[0]["z"].iloc[0] = 1.0

    assert df.loc[(0, 0), "z"] == df_original.loc[0, "z"]

    # correct setting
    df.loc[(0, 0), "z"] = 2
    result = df.loc[(0, 0), "z"]
    assert result == 2


def test_indexer_caching(monkeypatch):
    # GH5727
    # make sure that indexers are in the _internal_names_set
    size_cutoff = 20
    with monkeypatch.context():
        monkeypatch.setattr(libindex, "_SIZE_CUTOFF", size_cutoff)
        index = MultiIndex.from_arrays([np.arange(size_cutoff), np.arange(size_cutoff)])
        s = Series(np.zeros(size_cutoff), index=index)

        # setitem
        s[s == 0] = 1
    expected = Series(np.ones(size_cutoff), index=index)
    tm.assert_series_equal(s, expected)


def test_set_names_only_clears_level_cache():
    mi = MultiIndex.from_arrays([range(4), range(4)], names=["a", "b"])
    mi.dtypes
    mi.is_monotonic_increasing
    mi._engine
    mi.levels
    old_cache_keys = sorted(mi._cache.keys())
    assert old_cache_keys == ["_engine", "dtypes", "is_monotonic_increasing", "levels"]
    mi.names = ["A", "B"]
    new_cache_keys = sorted(mi._cache.keys())
    assert new_cache_keys == ["_engine", "dtypes", "is_monotonic_increasing"]
    new_levels = mi.levels
    tm.assert_index_equal(new_levels[0], RangeIndex(4, name="A"))
    tm.assert_index_equal(new_levels[1], RangeIndex(4, name="B"))
