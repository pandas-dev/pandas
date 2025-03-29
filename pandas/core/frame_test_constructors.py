import pandas as pd
def test_empty_df_preserve_col():
    rows = []
    df = pd.DataFrame.from_records(iter(rows), columns=['col_1', 'Col_2'], nrows=0)
    assert list(df.columns)==['col_1', 'Col_2']
    assert len(df) == 0
    