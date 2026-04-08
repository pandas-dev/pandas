import pandas as pd

df1 = pd.DataFrame({"A": [1, 2]}, index=[1, 2])
df2 = pd.DataFrame({"A": [100, 200]}, index=["1", "2"])

print("Before update:")
print(df1)

df1.update(df2)

print("\nAfter update:")
print(df1)

