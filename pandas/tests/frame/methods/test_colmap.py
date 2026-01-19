import pandas as pd


def test_colmap_lower():
    df = pd.DataFrame(columns=["A", "B"])
    df.colmap(str.lower)
    assert list(df.columns) == ["a", "b"]


def test_colmap_replace():
    df = pd.DataFrame(columns=["col 1", "col 2"])
    df.colmap(lambda c: c.replace(" ", "_"))
    assert list(df.columns) == ["col_1", "col_2"]


def test_colmap_returns_self():
    df = pd.DataFrame(columns=["A"])
    out = df.colmap(str.lower)
    assert out is df
