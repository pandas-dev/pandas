from typing import Any, Never, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csc_matrix, sparray, spmatrix
from scipy.sparse._base import _spbase

__all__ = ["clarkson_woodruff_transform"]

_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)

_ToInt: TypeAlias = int | npc.integer
_ToSparse: TypeAlias = _spbase[_ScalarT] | sparray[_ScalarT] | spmatrix[_ScalarT]

###

def cwt_matrix(n_rows: _ToInt, n_columns: _ToInt, rng: onp.random.ToRNG | None = None) -> csc_matrix[np.int_]: ...

#
@overload
def clarkson_woodruff_transform(
    input_matrix: _ToSparse[Never],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> csc_matrix[Any]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: _ToSparse[npc.integer | np.bool_],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> csc_matrix[np.int_]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: _ToSparse[np.float32 | np.float64],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> csc_matrix[np.float64]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: _ToSparse[npc.floating80],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> csc_matrix[np.longdouble]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: _ToSparse[np.complex64 | np.complex128],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> csc_matrix[np.complex128]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: _ToSparse[npc.complexfloating160],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> csc_matrix[np.clongdouble]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: onp.CanArrayND[Never],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> onp.ArrayND[Any]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: onp.ToIntND, sketch_size: _ToInt, rng: onp.random.ToRNG | None = None, *, seed: onp.random.ToRNG | None = None
) -> onp.ArrayND[np.int_]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: onp.ToArrayND[op.JustFloat, np.float16 | np.float32 | np.float64],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: onp.ToJustLongDoubleND,
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.longdouble]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: onp.ToArrayND[op.JustComplex, np.complex64 | np.complex128],
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def clarkson_woodruff_transform(
    input_matrix: onp.ToJustCLongDoubleND,
    sketch_size: _ToInt,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.clongdouble]: ...
