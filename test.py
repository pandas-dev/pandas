# contributing guide: https://pandas.pydata.org/docs/dev/development/contributing.html#pushing-your-changes
import pandas as pd

df = pd.DataFrame({
    "yes": [5000, 2, 3],
    "Byesyes": [4, 5, 6],
    "no": [7, 8, 9],
    "Byesno": [10, 11, 12],
    "YES": [13, 14, 15],
    "NO": [16, 17, 18],
    "YESYES": [19, 20, 21],
})

# Test the DataFrame creation
print(df.select_by_substr("skibidi"))