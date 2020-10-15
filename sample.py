import pandas as pd

df = pd.DataFrame(
    {"A": ["x", "x", "y"], "B": ["a", "b", "b"], "C": [1, 1, 1]}
).set_index(["A", "B"])
df["D"] = 2
result = df.unstack("B").fillna(0)
assert result.notna().all().all()
