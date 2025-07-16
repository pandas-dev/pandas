# mypy: ignore-errors
import numpy as np

import pandas as pd


def print_side_by_side(df1, df2):
    # Convert to string and split into lines
    df1_str = df1.to_string(index=False).split("\n")
    df2_str = df2.to_string(index=False).split("\n")

    # Pad lines to the same length for alignment
    max_len_1 = max(len(line) for line in df1_str)
    padded_df1 = [line.ljust(max_len_1) for line in df1_str]

    # Print side-by-side
    print("Result".ljust(max_len_1) + "  |  Expected")
    for line1, line2 in zip(padded_df1, df2_str):
        print(f"{line1}  |  {line2}")


data = np.arange(50).reshape(10, 5)
fill_val = 5

# data = pd.array(np.random.choice([True, False], size=(10, 5)), dtype="boolean")
# fill_val = True

# data = pd.array([i for i in range(50)], dtype="int")
# fill_val = 5

print(type(data))

columns = list("ABCDE")
df = pd.DataFrame(data, columns=columns)
for i in range(5):
    df.iat[i, i] = np.nan
    df.iat[i + 1, i] = np.nan
    df.iat[i + 4, i] = np.nan

df_base = df.iloc[:, :-1]
df_mult = df.iloc[:, [-1]]

mask = df.isna().values
mask = mask[:, :-1] & mask[:, [-1]]

df_result = df_base.mul(df_mult, axis=0, fill_value=fill_val)
df_expected = (df_base.fillna(fill_val).mul(df_mult.fillna(fill_val), axis=0)).mask(
    mask, np.nan
)

print_side_by_side(df_result, df_expected)
# tm.assert_frame_equal(df_result, df_expected)
