import pandas as pd
from markDown import markdownTable, markdownSeries

# Create a hypothetical dataset
data = {
    'Name': ['Alice', 'Bob', 'Charlie', 'David', 'Eve'],
    'Age': [25, 30, 35, 40, 45],
    'Gender': ['Female', 'Male', 'Male', 'Male', 'Female'],
    'Salary': [50000, 60000, 70000, 80000, 90000],
    'Department': ['HR', 'Finance', 'IT', 'Marketing', 'Operations']
}
df = pd.DataFrame(data)

# Create a new column 'Bonus' based on complex calculation
df['Bonus'] = df['Salary'] * 0.1 + df['Age'] * 0.05

# Perform some complex operations on the DataFrame
# For example, let's filter the DataFrame for individuals with Age > 30 and Salary > 60000
filtered_df = df[(df['Age'] > 30) & (df['Salary'] > 60000)]

# Convert the filtered DataFrame to a markdown table
filtered_markdown_table = markdownTable(filtered_df)
print("Filtered DataFrame as Markdown Table:")
print(filtered_markdown_table)

# Calculate the average Bonus for each gender
avg_bonus_by_gender = df.groupby('Gender')['Bonus'].mean()

# Convert the Series to a markdown table
avg_bonus_markdown = markdownSeries(avg_bonus_by_gender, col1='Gender', col2='Average Bonus')
print("\nAverage Bonus by Gender as Markdown Table:")
print(avg_bonus_markdown)
