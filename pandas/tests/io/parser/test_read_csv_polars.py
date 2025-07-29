import pytest


def test_read_csv_with_polars(tmp_path):
    pd = pytest.importorskip("pandas")

    # Create a simple CSV file
    file = tmp_path / "sample.csv"
    file.write_text("a,b\n1,2\n3,4")

    # Read using engine='polars'
    df = pd.read_csv(file, engine="polars")

    assert df.shape == (2, 2)
    assert list(df.columns) == ["a", "b"]
    assert df.iloc[0, 0] == 1
    assert df.iloc[1, 1] == 4
