# NOTE(scipy-stubs): This module only exists `if typing.TYPE_CHECKING: ...`, and has no stable API.

from typing import Literal, Protocol, SupportsIndex, TypeVar, final, type_check_only

import numpy as np
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._matrix import spmatrix

__all__ = "_CanAsAny", "_CanStack", "_CanStackAs", "_Format", "_Sparse2D", "_ToShape1D", "_ToShape2D"

###

# sparray and spmatrix must be included for them to be assignable to this type alias
type _Sparse2D[_ScalarT: npc.number | np.bool] = _spbase[_ScalarT, tuple[int, int]] | sparray[_ScalarT] | spmatrix[_ScalarT]

type _ToShape1D = tuple[SupportsIndex]  # ndim == 1
type _ToShape2D = tuple[SupportsIndex, SupportsIndex]  # ndim == 2

type _Format = Literal["bsr", "coo", "csc", "csr", "dia", "dok", "lil"]

###
# Interfaces for emulated dependent associated types

_AssocT_co = TypeVar("_AssocT_co", covariant=True)
_ScalarT_contra = TypeVar("_ScalarT_contra", bound=npc.number | np.bool, contravariant=True)

# same kind, same scalar-type, 2d shape-type
@final
@type_check_only
class _CanStack(Protocol[_AssocT_co]):
    @type_check_only
    def __assoc_stacked__(self, /) -> _AssocT_co: ...

# same kind, mapped scalar-type, 2d shape-type
@final
@type_check_only
class _CanStackAs(Protocol[_ScalarT_contra, _AssocT_co]):
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT_contra, /) -> _AssocT_co: ...

# same kind, `Any` scalar-type, same shape-type
@final
@type_check_only
class _CanAsAny(Protocol[_AssocT_co]):
    @type_check_only
    def __assoc_as_any__(self, /) -> _AssocT_co: ...
