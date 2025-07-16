# mypy: ignore-errors
import numpy as np

import pandas as pd
import pandas._testing as tm


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


# data = np.arange(50).reshape(10, 5)
# fill_val = 5

# data = pd.array(np.random.choice([True, False], size=(10, 5)), dtype="boolean")
# fill_val = True

data = np.arange(50).reshape(10, 5)
# data_mult = pd.array([i for i in range(10)], dtype=tm.SIGNED_INT_NUMPY_DTYPES[0])
data_mult = pd.array(list(range(10)), dtype=tm.SIGNED_INT_EA_DTYPES[0])
fill_val = 5

# print(tm.ALL_INT_DTYPES)
# print(tm.SIGNED_INT_EA_DTYPES)
# tm.SIGNED_INT_NUMPY_DTYPES[0]
print(type(data_mult))

# TODO masking not working with EA with dim > 1
# NOTE currently trying to get EA testing set up

columns = list("ABCDE")
df_base = pd.DataFrame(data, columns=columns)
for i in range(5):
    df_base.iat[i, i] = np.nan
    df_base.iat[i + 1, i] = np.nan
    df_base.iat[i + 4, i] = np.nan

mask = df_base.isna().values

data_mult_re = data_mult.reshape(10, 1)
mask = mask[:, :-1] & data_mult_re

df_result = df_base.mul(data_mult, axis=0, fill_value=fill_val)
print(df_result)
# df_expected = (df_base.fillna(fill_val).mul(data_mult.fillna(fill_val),
# axis=0)).mask(mask, np.nan)

# print_side_by_side(df_result, df_expected)
# # tm.assert_frame_equal(df_result, df_expected)
