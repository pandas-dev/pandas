import pytest

def test_read_csv_with_polars(tmp_path):
    pl = pytest.importorskip("polars")
    pd = pytest.importorskip("pandas")

    file = tmp_path / "data.csv"
    file.write_text("a,b\n1,2\n3,4")

    df = pd.read_csv(file, engine="polars")
    assert df.shape == (2, 2)
    assert list(df.columns) == ["a", "b"]
