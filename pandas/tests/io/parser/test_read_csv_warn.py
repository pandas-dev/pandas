import pandas as pd

def my_bad_line_handler(bad_line):
    print("Bad line encountered:", bad_line)
    return None

df = pd.read_csv(
    "test.csv",  # make sure this file exists in same folder or adjust the path
    on_bad_lines=my_bad_line_handler,
    index_col=0,
    engine="python",  # ✅ add this line
)
print(df)
