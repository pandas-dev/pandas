import pandas as pd


import pandas as pd

# Create a pandas DataFrame with a date column
df = pd.DataFrame({'date': ['4/03/2022', '5/03/2022', '6/03/2022'], 'patients': [16, 19, 11]})

# Convert the date column to a pandas `DatetimeIndex` object with the day first
#df['year'] = pd.to_datetime(df['date'], yearFirst=False)
df['year'] = pd.to_datetime(df['date'], dayfirst=True)

# Print the DataFrame
print(df)

# Create a pandas DataFrame with a date column
df = pd.DataFrame({'date': ['04/03/2004', '03/03/2005', '03/03/2007'], 'patients': [16, 19, 11]})

# Convert the date column to a pandas `DatetimeIndex` object with the day first
#df['year'] = pd.to_datetime(df['date'], yearFirst=False)
df['date'] = pd.to_datetime(df['date'], format='%d/%m/%Y').dt.strftime('%d-%m-%Y')

# Print the DataFrame
print(df)




