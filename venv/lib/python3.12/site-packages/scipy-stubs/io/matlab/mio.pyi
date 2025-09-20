# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

__all__ = ["loadmat", "savemat", "whosmat"]

@deprecated("will be removed in SciPy v2.0.0")
def loadmat(
    file_name: object, mdict: object = ..., appendmat: object = ..., *, spmatrix: bool = True, **kwargs: object
) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def savemat(
    file_name: object,
    mdict: object,
    appendmat: object = ...,
    format: object = ...,
    long_field_names: object = ...,
    do_compression: object = ...,
    oned_as: object = ...,
) -> None: ...
@deprecated("will be removed in SciPy v2.0.0")
def whosmat(file_name: object, appendmat: object = ..., **kwargs: object) -> Any: ...
