import pandas as pd
import numpy as np

def advanced_groupby(df, groupby_cols, agg_funcs):
    """
    Perform advanced groupby with nested grouping and custom aggregations.

    Parameters:
    df (DataFrame): The DataFrame to aggregate.
    groupby_cols (list of list): List of lists, where each inner list represents a level of grouping.
    agg_funcs (list of dict): List of dictionaries, where each dictionary corresponds to custom aggregation functions for each grouping level.

    Returns:
    DataFrame: Aggregated DataFrame.
    """
    grouped = df
    for i, cols in enumerate(groupby_cols):
        grouped = grouped.groupby(cols).agg(agg_funcs[i]).reset_index()

    return grouped

df = pd.DataFrame({
    'Region': ['North', 'North', 'South', 'South'],
    'Country': ['USA', 'Canada', 'USA', 'Canada'],
    'City': ['New York', 'Toronto', 'Los Angeles', 'Vancouver'],
    'Sales': [100, 200, 150, 250],
    'Profit': [10, 20, 15, 25]
})

groupby_cols = [['Region'], ['Region', 'Country']]

agg_funcs = [
    {'Sales': np.sum, 'Profit': np.mean},
    {'Sales': np.sum, 'Profit': np.mean}
]

result = advanced_groupby(df, groupby_cols, agg_funcs)
print(result)