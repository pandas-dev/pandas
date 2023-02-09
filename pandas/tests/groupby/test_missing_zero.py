import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.testing import assert_frame_equal, assert_series_equal


np.random.seed(42)

dates = pd.date_range('1/1/2001', '31/12/2001', freq = 'd')
sales = [np.random.randint(100) for _ in range(len(dates))]
product = [['A', 'B', 'C'][np.random.randint(3)] for _ in range(len(dates))]

full_df = pd.DataFrame({'Dates': dates,
                   'Sales': sales,
                   'Product': product
                   })
march = full_df.Dates.dt.month == 3
incomplete_df = full_df[~march]

monthly = pd.Grouper(key='Dates', freq='M')

product_df = incomplete_df[incomplete_df["Product"] == "A"]

# sum sales
display_df = product_df.groupby(monthly)['Sales'].sum()
print("-- OK --")
print(display_df)

print("-- BUG --")
# bug_display_df = product_df.groupby(['Product', monthly])['Sales'].sum()
bug_display_df = product_df.groupby(['Product', monthly])['Sales'].sum().droplevel("Product", axis="index")
print(bug_display_df)

assert assert_series_equal(display_df, bug_display_df)
