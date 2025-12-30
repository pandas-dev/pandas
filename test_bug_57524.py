import pandas as pd

s1 = pd.Series(
    [23, 22, 21], index=pd.Index(["a", "b", "c"], name="index a"), dtype=pd.Int64Dtype()
)
s2 = pd.Series(
    [21, 22, 23],
    index=pd.Index(["a", "b", "c"], name="index b", dtype=pd.StringDtype()),
    dtype=pd.Int64Dtype(),
)

print("s1:", s1)
print("s1.index:", s1.index)
print("s1.index dtype:", s1.index.dtype)
print()
print("s2:", s2)
print("s2.index:", s2.index)
print("s2.index dtype:", s2.index.dtype)
print()

try:
    result = s1 - s2
    print("Result:", result)
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
