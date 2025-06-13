import pandas as pd


def test_basemaskedarray_map():
    for dtype, data, expected_data in [
        ("Int32", [1, 2, None, 4], [2, 3, pd.NA, 5]),
    ]:
        s = pd.Series(data, dtype=dtype)

        def transform(x):
            if x is None:
                return x
            return x + 1

        result = s.map(transform)
        expected = pd.Series(expected_data, dtype=result.dtype)

        assert result.tolist() == expected.tolist()
