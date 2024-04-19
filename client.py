import pandas as pd

df1 = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
store = pd.HDFStore("store.h5", "w")  # doctest: +SKIP
store.put("data", df1, format="fixed")  # doctest: +SKIP
df2 = pd.DataFrame([[5, 6], [7, 8]], columns=["A", "B"])
store.append("data", df2)  # doctest: +SKIP
store.close()  # doctest: +SKIP