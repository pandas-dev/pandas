df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
df.plot(kind='bar')
pd.set_option('display.max_rows', 10)  # Display up to 10 rows
# Read CSV
df = pd.read_csv('data.csv')
# Write to Excel
df.to_excel('output.xlsx', index=False)
df['A_squared'] = df['A'].apply(lambda x: x**2)
s = pd.Series([1, 2, pd.NA], dtype="Int64")
