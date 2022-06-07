import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

np.random.seed(42)

dates = pd.date_range('1/1/2001', '31/12/2001', freq = 'd')
sales = [np.random.randint(100) for _ in range(len(dates))]
product = [['A', 'B', 'C'][np.random.randint(3)] for _ in range(len(dates))]

df = pd.DataFrame({'Dates': dates,
                   'Sales': sales,
                   'Product': product
                   })
march = df.Dates.dt.month == 3
df = df[~march]

print(df)

monthly = pd.Grouper(key='Dates', freq='M')
sum_sales = df.groupby(['Product', monthly])['Sales'].sum()

print(sum_sales)
