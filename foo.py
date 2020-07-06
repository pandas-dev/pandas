import pandas as pd

# https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_json.html
#df = pd.DataFrame([['a', 'b'], ['c', 'd']],index=['row 1', 'row 2'],columns=['col 1', 'col 2'])

df = pd.DataFrame([["a", "b"], ["c", "d"]], columns=["x", "x"])
print(df.to_json(orient="split"))


