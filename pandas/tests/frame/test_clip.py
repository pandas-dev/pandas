import pandas as pd

def test_clip_lower_greater_than_upper():
    df = pd.DataFrame({'A': [1, 5, 10]})
    result = df.clip(lower=10, upper=5)
    expected = pd.DataFrame({'A': [5, 5, 10]})
    pd.testing.assert_frame_equal(result, expected)
