import pandas as pd
from markDown import markdownTable, markdownSeries

# Generic Grade book data
data = {
    'Student': [f'Student {i}' for i in range(1, 21)],
    'Math': [round(50 + i + (i % 3) * 5, 1) for i in range(1, 21)],
    'Science': [round(55 + i + (i % 4) * 3.5, 1) for i in range(1, 21)],
    'English': [round(60 + i + (i % 2) * 2.7, 1) for i in range(1, 21)],
    'History': [round(65 + i + (i % 5) * 4.2, 1) for i in range(1, 21)]
}
df = pd.DataFrame(data)

# Convert the DataFrame to a markdown table
markdown_table = markdownTable(df)
print(f"Convert the DataFrame to a markdown table\n{markdown_table}")

# Convert the 'Math' Series to a markdown table
math_series_markdown = markdownSeries(df['Math'], col1='Student', col2='Math Grade') # noqa
print(f"Convert the 'Math' Series to a markdown table\n{math_series_markdown}")

# Print value counts of 'Student' column in a markdown table format
value_counts_series = df['Student'].value_counts().sort_values(ascending=False)
markdown_table = markdownSeries(value_counts_series, col1='Student', col2='Count') # noqa
print(f"Print value counts of 'Student' column in a markdown table format\n{markdown_table}") # noqa