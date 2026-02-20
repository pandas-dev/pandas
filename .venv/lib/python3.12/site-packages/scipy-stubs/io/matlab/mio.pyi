# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

__all__ = ["loadmat", "savemat", "whosmat"]

@deprecated("will be removed in SciPy v2.0.0")
def loadmat(
    file_name: object, mdict: object = None, appendmat: object = True, *, spmatrix: bool = True, **kwargs: object
) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def savemat(
    file_name: object,
    mdict: object,
    appendmat: object = True,
    format: object = "5",
    long_field_names: object = False,
    do_compression: object = False,
    oned_as: object = "row",
) -> None: ...
@deprecated("will be removed in SciPy v2.0.0")
def whosmat(file_name: object, appendmat: object = True, **kwargs: object) -> Any: ...
