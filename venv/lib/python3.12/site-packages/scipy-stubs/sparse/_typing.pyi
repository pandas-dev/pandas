# NOTE(scipy-stubs): This module only exists `if typing.TYPE_CHECKING: ...`, and has no stable API.

from typing import Literal, Protocol, SupportsIndex, TypeAlias, TypeVar, final, type_check_only
from typing_extensions import TypeAliasType

import numpy as np
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._matrix import spmatrix

__all__ = "_CanStack", "_CanStackAs", "_Format", "_Sparse2D", "_ToShape1D", "_ToShape2D"

###

_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)

# sparray and spmatrix must be included for them to be assignable to this type alias
_Sparse2D = TypeAliasType(
    "_Sparse2D",
    _spbase[_ScalarT, tuple[int, int]] | sparray[_ScalarT] | spmatrix[_ScalarT],
    type_params=(_ScalarT,),
)  # fmt: skip

_ToShape1D: TypeAlias = tuple[SupportsIndex]  # ndim == 1
_ToShape2D: TypeAlias = tuple[SupportsIndex, SupportsIndex]  # ndim == 2

_Format: TypeAlias = Literal["bsr", "coo", "csc", "csr", "dia", "dok", "lil"]

###
# Interfaces for emulated dependent associated types

_AssocT_co = TypeVar("_AssocT_co", covariant=True)
_ScalarT_contra = TypeVar("_ScalarT_contra", bound=npc.number | np.bool_, contravariant=True)

@final
@type_check_only
class _CanStack(Protocol[_AssocT_co]):
    @type_check_only
    def __assoc_stacked__(self, /) -> _AssocT_co: ...

@final
@type_check_only
class _CanStackAs(Protocol[_ScalarT_contra, _AssocT_co]):
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT_contra, /) -> _AssocT_co: ...
