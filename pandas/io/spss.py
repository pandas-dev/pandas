from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Sequence,
)

from pandas._config import using_nullable_dtypes

from pandas._libs import lib
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.inference import is_list_like

from pandas.io.common import stringify_path

if TYPE_CHECKING:
    from pathlib import Path

    from pandas import DataFrame


def read_spss(
    path: str | Path,
    usecols: Sequence[str] | None = None,
    convert_categoricals: bool = True,
    use_nullable_dtypes: bool | lib.NoDefault = lib.no_default,
) -> DataFrame:
    """
    Load an SPSS file from the file path, returning a DataFrame.

    Parameters
    ----------
    path : str or Path
        File path.
    usecols : list-like, optional
        Return a subset of the columns. If None, return all columns.
    convert_categoricals : bool, default is True
        Convert categorical columns into pd.Categorical.
    use_nullable_dtypes : bool = False
        Whether to use nullable dtypes as default when reading data. If
        set to True, nullable dtypes are used for all dtypes that have a nullable
        implementation, even if no nulls are present.

        .. note::

            The nullable dtype implementation can be configured by calling
            ``pd.set_option("mode.dtype_backend", "pandas")`` to use
            numpy-backed nullable dtypes or
            ``pd.set_option("mode.dtype_backend", "pyarrow")`` to use
            pyarrow-backed nullable dtypes (using ``pd.ArrowDtype``).

        .. versionadded:: 2.0

    Returns
    -------
    DataFrame
    """
    pyreadstat = import_optional_dependency("pyreadstat")

    use_nullable_dtypes = (
        use_nullable_dtypes
        if use_nullable_dtypes is not lib.no_default
        else using_nullable_dtypes()
    )

    if usecols is not None:
        if not is_list_like(usecols):
            raise TypeError("usecols must be list-like.")
        usecols = list(usecols)  # pyreadstat requires a list

    df, _ = pyreadstat.read_sav(
        stringify_path(path), usecols=usecols, apply_value_formats=convert_categoricals
    )
    if use_nullable_dtypes:
        df = df.convert_dtypes()
    return df
