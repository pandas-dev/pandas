import pandas as pd


def markdownTable(df, *column_names):
    if not column_names:
        column_names = df.columns.tolist()
    table_markdown = "| " + " | ".join(column_names) + " |\n"
    table_markdown += "| " + " | ".join(["---"] * len(column_names)) + " |\n"
    for index, row in df.iterrows():
        row_data = [str(row[col]) for col in column_names]
        table_markdown += "| " + " | ".join(row_data) + " |\n"
    return table_markdown


def markdownSeries(series, col1=None, col2=None):
    if not col1:
        col1 = series.index.name if series.index.name else "Index"
    if not col2:
        col2 = series.name if series.name else "Value"
    table_markdown = f"| {col1} | {col2} |\n|---|---|\n"
    for index, value in series.items():
        table_markdown += f"| {index} | {value} |\n"
    return table_markdown


pd.DataFrame.markdownTable = markdownTable
pd.Series.markdownTable = markdownSeries
