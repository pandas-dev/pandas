def read_spss(path):
    """
    Load an parquet object from the file path, returning a DataFrame.

    .. versionadded 0.24.3

    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    DataFrame
    """

    from pyreadstat import read_sav
    df, _ = read_sav(path)
    return df
