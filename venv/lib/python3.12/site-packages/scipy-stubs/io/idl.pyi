# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = ["readsav"]

@deprecated("will be removed in SciPy v2.0.0")
def readsav(
    file_name: object, idict: object = ..., python_dict: object = ..., uncompressed_file_name: object = ..., verbose: object = ...
) -> object: ...
