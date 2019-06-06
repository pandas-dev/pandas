from pathlib import Path
from typing import Optional, Sequence, Union

from pandas.core.api import DataFrame


def read_spss(path: Union[str, Path],
              usecols: Optional[Sequence[str]] = None,
              categorical: bool = True) -> DataFrame:
    """
    Load an SPSS file from the file path, returning a DataFrame.

    .. versionadded 0.25.0

    Parameters
    ----------
    path : string or Path
        File path
    usecols : list-like or None
        Return a subset of the columns. If None, return all columns.
    categorical : bool
        Convert categorical columns into pd.Categorical.

    Returns
    -------
    DataFrame
    """
    try:
        from pyreadstat import read_sav
    except ImportError:
        raise ImportError("pyreadstat is required to read SPSS .sav files.")
    if usecols is not None:
        if isinstance(usecols, str):
            usecols = [usecols]
    df, _ = read_sav(path, usecols=usecols,
                     apply_value_formats=categorical)
    return df
