from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Sequence,
)

from pandas._libs import lib
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc
from pandas.util._validators import check_dtype_backend

from pandas.core.dtypes.inference import is_list_like
from pandas.core.shared_docs import _shared_docs

from pandas.io.common import stringify_path

if TYPE_CHECKING:
    from pathlib import Path

    from pandas._typing import DtypeBackend

    from pandas import DataFrame


@doc(dtype_backend_options=_shared_docs["dtype_backend_options"])
def read_spss(
    path: str | Path,
    usecols: Sequence[str] | None = None,
    convert_categoricals: bool = True,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
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
    {dtype_backend_options}

        .. versionadded:: 2.0

    Returns
    -------
    DataFrame
    """
    pyreadstat = import_optional_dependency("pyreadstat")
    check_dtype_backend(dtype_backend)

    if usecols is not None:
        if not is_list_like(usecols):
            raise TypeError("usecols must be list-like.")
        usecols = list(usecols)  # pyreadstat requires a list

    df, _ = pyreadstat.read_sav(
        stringify_path(path), usecols=usecols, apply_value_formats=convert_categoricals
    )
    if dtype_backend is not lib.no_default:
        df = df.convert_dtypes(dtype_backend=dtype_backend)
    return df
