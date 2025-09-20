# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Final
from typing_extensions import deprecated

from ._arffread import MetaData as _MetaData

__all__ = ["ArffError", "MetaData", "ParseArffError", "loadarff"]

_msg: Final = "will be removed in SciPy v2.0.0"

@deprecated(_msg)
class ArffError(OSError): ...

@deprecated(_msg)
class ParseArffError(ArffError): ...  # pyright: ignore[reportDeprecated]

@deprecated(_msg)
class MetaData(_MetaData): ...

@deprecated(_msg)
def loadarff(f: object) -> object: ...
