import pandas as pd

pd.Series([2**64], dtype=object).convert_dtypes()
# OverflowError: Python int too large to convert to C long

pd.Series([-(2**63) - 1], dtype=object).convert_dtypes()
# OverflowError: Python int too large to convert to C long