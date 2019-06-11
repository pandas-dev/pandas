from pathlib import Path
from typing import Optional, Sequence, Union

from pandas.api.types import is_list_like
from pandas.core.api import DataFrame


def read_spss(path: Union[str, Path],
              usecols: Optional[Sequence[str]] = None,
              convert_categoricals: bool = True) -> DataFrame:
    """
    Load an SPSS file from the file path, returning a DataFrame.

    .. versionadded 0.25.0

    Parameters
    ----------
    path : string or Path
        File path
    usecols : list-like, optional
        Return a subset of the columns. If None, return all columns.
    convert_categoricals : bool, default is True
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
        if not is_list_like(usecols):
            raise TypeError("usecols must be list-like.")
        else:
            usecols = list(usecols)  # pyreadstat requires a list

    df, _ = read_sav(path, usecols=usecols,
                     apply_value_formats=convert_categoricals)
    return df
