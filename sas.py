import pandas as pd
df = pd.read_sas('airline.sas7bdat')
indexNames = df['YEAR'].index
ds = df.drop(indexNames)
print(ds)
#ds.to_csv("file.csv")



