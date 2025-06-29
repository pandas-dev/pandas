"""
Utility function for quick dataset diagnostics.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from pandas.util._decorators import set_module

if TYPE_CHECKING:
    from pandas import DataFrame


@set_module("pandas")
def check(df: DataFrame, round_digits: int = 2) -> DataFrame:
    """
    Provide a column-wise summary of DataFrame structure for quick diagnostics.

    This function combines several common exploratory data analysis operations
    into a single diagnostic summary, including unique values, non-null counts,
    missing value counts, and missing percentages.

    Parameters
    ----------
    df : DataFrame
        The DataFrame to analyze.
    round_digits : int, default 2
        Number of decimal places to round the missing percentage to.

    Returns
    -------
    DataFrame
        A DataFrame with columns:
        - unique: Number of unique values per column
        - non_null: Number of non-null values per column  
        - missing: Number of missing values per column
        - missing_pct: Percentage of missing values per column

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'A': [1, 2, None, 4],
    ...     'B': ['x', 'y', 'x', None],
    ...     'C': [1.0, 2.0, 3.0, 4.0]
    ... })
    >>> pd.check(df)
       unique  non_null  missing  missing_pct
    A       3         3        1        25.00
    B       2         3        1        25.00
    C       4         4        0         0.00
    """
    import pandas as pd
    
    # Calculate basic statistics for each column
    unique_counts = df.nunique()
    non_null_counts = df.count()
    missing_counts = df.isnull().sum()
    total_rows = len(df)
    missing_pct = (missing_counts / total_rows * 100).round(round_digits)
    
    # Create the result DataFrame
    result = pd.DataFrame({
        'unique': unique_counts,
        'non_null': non_null_counts,
        'missing': missing_counts,
        'missing_pct': missing_pct
    })
    
    return result