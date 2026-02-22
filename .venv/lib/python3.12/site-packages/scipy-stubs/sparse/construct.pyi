# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

from . import _bsr, _coo, _csc, _csr, _dia

__all__ = [
    "block_diag",
    "bmat",
    "bsr_matrix",
    "coo_matrix",
    "csc_matrix",
    "csr_matrix",
    "dia_matrix",
    "diags",
    "eye",
    "get_index_dtype",
    "hstack",
    "identity",
    "issparse",
    "kron",
    "kronsum",
    "rand",
    "random",
    "spdiags",
    "vstack",
]

@deprecated("will be removed in SciPy v2.0.0")
class coo_matrix(_coo.coo_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class csc_matrix(_csc.csc_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class csr_matrix(_csr.csr_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class dia_matrix(_dia.dia_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
class bsr_matrix(_bsr.bsr_matrix): ...

@deprecated("will be removed in SciPy v2.0.0")
def issparse(x: object) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def vstack(blocks: object, format: object = None, dtype: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def spdiags(data: object, diags: object, m: object = None, n: object = None, format: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def random(
    m: object,
    n: object,
    density: object = 0.01,
    format: object = "coo",
    dtype: object = None,
    rng: object = None,
    data_rvs: object = None,
    *,
    random_state: object = None,
) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def rand(
    m: object,
    n: object,
    density: object = 0.01,
    format: object = "coo",
    dtype: object = None,
    rng: object = None,
    *,
    random_state: object = None,
) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def kron(A: object, B: object, format: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def kronsum(A: object, B: object, format: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def identity(n: object, dtype: object = "d", format: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def block_diag(mats: object, format: object = None, dtype: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def hstack(blocks: object, format: object = None, dtype: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def eye(m: object, n: object = None, k: object = 0, dtype: object = ..., format: object = None) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def diags(diagonals: object, offsets: object = 0, shape: object = None, format: object = None, dtype: object = ...) -> Any: ...
@deprecated("will be removed in SciPy v2.0.0")
def bmat(blocks: object, format: object = None, dtype: object = None) -> Any: ...

# sputils
@deprecated("will be removed in SciPy v2.0.0")
def get_index_dtype(arrays: object = (), maxval: object = None, check_contents: object = False) -> Any: ...
