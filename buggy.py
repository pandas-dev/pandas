import pandas as pd
import numpy as np
# df = pd.DataFrame({"A": ["1", "", "3"]}, dtype="string")
# try:
#     whydontyouwork1 = df.where(df != "", np.nan)
#     print(whydontyouwork1)
#     arr = whydontyouwork1["A"]._values
#     print("arr",arr)
#     print("type",type(arr[1]))
# except Exception as e:
#     print(e)
# whydontyouwork2=df.where(df != "", np.nan, inplace=True)
# print(whydontyouwork2)
# arr = whydontyouwork2["A"]._values
# print("arr",arr)
# print("type:",type(arr[1]))
#################
# import pandas as pd
# import numpy as np
# df = pd.DataFrame({"HEREITIS": ["1", "", "3"]}, dtype="string")
# df2 = pd.DataFrame({"HEREITIS": ["1", "", "3"]}, dtype="string")
# whydontyouwork1 = df.where(df != "", np.nan)
# print(df2.where(df != "", np.nan, inplace=True))
# print(whydontyouwork1)
# #print(whydontyouwork2)
##################
# df = pd.DataFrame({"A": ["1", "", "3"]}, dtype="string")
# try:
#     result = df.where(df != "", np.nan)
#     arr = result["A"]._values
#     print(arr)
#     print(type(arr[1]))
# except Exception as e:
#     print("hi")
# df.where(df != "", np.nan, inplace=True)
# print(df)
# arr = df["A"]._values
# print(arr)
# print(type(arr[1]))
#################
df = pd.DataFrame({"HEREITIS": ["1", "", "3"]}, dtype="string")
df2 = pd.DataFrame({"HEREITIS": ["1", "", "3"]}, dtype="string")

df2.where(df != "", np.nan, inplace=True)
print(df2)
df = df.where(df != "", np.nan)
print(df)
#print(whydontyouwork2)