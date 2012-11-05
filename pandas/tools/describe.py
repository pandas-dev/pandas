from pandas.core.series import Series


def value_range(df):
    """
    Return the minimum and maximum of a dataframe in a series object

    Parameters
    ----------
    df : DataFrame

    Returns
    -------
    (maximum, minimum) : Series

    """
    return Series((min(df.min()), max(df.max())), ('Minimum', 'Maximum'))
