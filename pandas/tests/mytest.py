import pandas as pd

def test_me():
    pd.eval(
        """
        A = df.A - df.B
        B = df.A + df.B
        """,
        target=pd.DataFrame(),
    )