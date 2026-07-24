# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

import optype as op

__all__ = ["hb_read", "hb_write"]

type _NoValueType = op.JustObject

@deprecated("will be removed in SciPy v2.0.0")
def hb_read(path_or_open_file: object, *, spmatrix: bool | _NoValueType = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hb_write(path_or_open_file: object, m: object, hb_info: object = None) -> object: ...
