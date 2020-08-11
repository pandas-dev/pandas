import pandas as pd

df1 = pd.DataFrame(
    {
        "a": [3, 4, 5, 3, 5, 4, 1, 2, 3],
        "b": [1, 3, 4, 5, 6, 5, 4, 4, 4],
        "c": list("aabaaddce"),
        "d": [3, 4, 5, 3, 5, 4, 1, 2, 3],
        "e": [1, 3, 4, 5, 6, 5, 4, 4, 4],
        "f": list("aabaaddce"),
    }
)

assert len(df1.groupby("b").agg(["cumsum"])) == len(df1)
