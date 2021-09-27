import pandas as pd
from pandas.core.indexes.datetimes import date_range

#dr = date_range(start=pd.Timestamp("2000-01-05 05:09:15.13"), freq="1MS", periods=3)
# #print(dr)
hours = pd.Series(pd.Timestamp("2000-01-05 05:09:15.13", tz="UTC"))
# print(hours.dt.ceil('4D'))
print("#############################################################################################")
print(hours.dt.ceil('4MS'))