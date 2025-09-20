from typing import Final, Generic, Literal, Self, TypeAlias, TypedDict, final, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype.numpy as onp

__all__ = [
    "MDTYPES",
    "NP_TO_MTYPES",
    "NP_TO_MXTYPES",
    "OPAQUE_DTYPE",
    "MatlabFunction",
    "MatlabObject",
    "MatlabOpaque",
    "codecs_template",
    "mat_struct",
    "mclass_dtypes_template",
    "mclass_info",
    "mdtypes_template",
    "miCOMPRESSED",
    "miDOUBLE",
    "miINT8",
    "miINT16",
    "miINT32",
    "miINT64",
    "miMATRIX",
    "miSINGLE",
    "miUINT8",
    "miUINT16",
    "miUINT32",
    "miUINT64",
    "miUTF8",
    "miUTF16",
    "miUTF32",
    "mxCELL_CLASS",
    "mxCHAR_CLASS",
    "mxDOUBLE_CLASS",
    "mxFUNCTION_CLASS",
    "mxINT8_CLASS",
    "mxINT16_CLASS",
    "mxINT32_CLASS",
    "mxINT64_CLASS",
    "mxOBJECT_CLASS",
    "mxOBJECT_CLASS_FROM_MATRIX_H",
    "mxOPAQUE_CLASS",
    "mxSINGLE_CLASS",
    "mxSPARSE_CLASS",
    "mxSTRUCT_CLASS",
    "mxUINT8_CLASS",
    "mxUINT16_CLASS",
    "mxUINT32_CLASS",
    "mxUINT64_CLASS",
]

# NOTE: explicit covariance can't be used, because shape types are invariant in `numpy<2.1`
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=tuple[int, ...])

@type_check_only
class _CodecTemplateValue(TypedDict):
    codec: _Codec
    width: Literal[1, 2, 4]

@type_check_only
class _MDTypesValue(TypedDict):
    dtypes: dict[str, np.dtype[np.generic]]
    classes: dict[str, np.dtype[np.generic]]
    codecs: dict[_MCodec, _CodecBO]

_MDTypes = TypedDict("_MDTypes", {"<": _MDTypesValue, ">": _MDTypesValue})

@type_check_only
class _NP2M(TypedDict):
    f8: Literal[9]
    c32: Literal[9]
    c24: Literal[9]
    c16: Literal[9]
    f4: Literal[7]
    c8: Literal[7]
    i8: Literal[12]
    i4: Literal[5]
    i2: Literal[3]
    i1: Literal[1]
    u8: Literal[13]
    u4: Literal[6]
    u2: Literal[4]
    u1: Literal[2]
    S1: Literal[2]
    U1: Literal[17]
    b1: Literal[2]

@type_check_only
class _NP2MX(TypedDict):
    f8: Literal[6]
    c32: Literal[6]
    c24: Literal[6]
    c16: Literal[6]
    f4: Literal[7]
    c8: Literal[7]
    i8: Literal[14]
    i4: Literal[12]
    i2: Literal[10]
    i1: Literal[8]
    u8: Literal[15]
    u4: Literal[13]
    u2: Literal[11]
    u1: Literal[9]
    S1: Literal[9]
    b1: Literal[9]

# NOTE: TypedDict doesn't support integer keys (but the python core devs are literally incapable of admitting their mistakes)
_CodecsTemplate: TypeAlias = dict[_MCodec, _CodecTemplateValue]

_ByteOrder: TypeAlias = Literal["<", ">"]
_Codec: TypeAlias = Literal["utf_8", "utf_16", "utf_32"]
_CodecBO: TypeAlias = Literal["utf_8", "utf_16_le", "utf_16_be", "utf_32_le", "utf_32_be"]

_MCodec: TypeAlias = Literal[16, 17, 18]
_MType: TypeAlias = Literal[1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 14, 15, _MCodec]
_Number: TypeAlias = Literal["i1", "i2", "i4", "i8", "u1", "u2", "u4", "u8", "f4", "f8"]
_MDTypeTemplateKey: TypeAlias = Literal[_MType, "array_flags", "file_header", "tag_full", "tag_smalldata", "U1"]
_MDTypeTemplateValueKey: TypeAlias = Literal[
    "description", "subsystem_offset", "version", "endian_test",
    "mdtype", "byte_count",
    "byte_count_mdtype", "data",
    "data_type", "flags_class", "nzmax",
]  # fmt: skip
_MDTypeTemplateValueValue: TypeAlias = Literal["S2", "u2", "s4", "u4", "i8", "S116"]
_MDTypeTemplateValueItems: TypeAlias = list[tuple[_MDTypeTemplateValueKey, _MDTypeTemplateValueValue]]
_MDTypeTemplateValue: TypeAlias = _Number | _MDTypeTemplateValueItems

_MXNumber: TypeAlias = Literal[6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
_MXType: TypeAlias = Literal[1, 2, 3, 4, 5, _MXNumber, 16, 17, 18]
_MXName: TypeAlias = Literal[
    "cell", "struct", "object", "char", "sparse",
    "single", "double", "int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64",
    "function", "opaque",
]  # fmt: skip

###

MDTYPES: Final[_MDTypes] = ...
NP_TO_MTYPES: Final[_NP2M] = ...
NP_TO_MXTYPES: Final[_NP2MX] = ...
# TODO(jorenham): use `np.dtypes.VoidDType[Literal[32]]` once we have `numpy >= 2.1.0`
OPAQUE_DTYPE: Final[np.dtype[np.void]] = ...

codecs_template: Final[_CodecsTemplate] = ...
mdtypes_template: Final[dict[_MDTypeTemplateKey, _MDTypeTemplateValue]] = ...
mclass_info: Final[dict[_MXType, _MXName]] = ...
mclass_dtypes_template: Final[dict[_MXNumber, _Number]] = ...

miINT8: Final[_MType] = 1
miUINT8: Final[_MType] = 2
miINT16: Final[_MType] = 3
miUINT16: Final[_MType] = 4
miINT32: Final[_MType] = 5
miUINT32: Final[_MType] = 6
miSINGLE: Final[_MType] = 7
miDOUBLE: Final[_MType] = 9
miINT64: Final[_MType] = 12
miUINT64: Final[_MType] = 13
miMATRIX: Final[_MType] = 14
miCOMPRESSED: Final[_MType] = 15
miUTF8: Final[_MType] = 16
miUTF16: Final[_MType] = 17
miUTF32: Final[_MType] = 18

mxCELL_CLASS: Final[_MXType] = 1
mxSTRUCT_CLASS: Final[_MXType] = 2
mxOBJECT_CLASS: Final[_MXType] = 3
mxCHAR_CLASS: Final[_MXType] = 4
mxSPARSE_CLASS: Final[_MXType] = 5
mxDOUBLE_CLASS: Final[_MXType] = 6
mxSINGLE_CLASS: Final[_MXType] = 7
mxINT8_CLASS: Final[_MXType] = 8
mxUINT8_CLASS: Final[_MXType] = 9
mxINT16_CLASS: Final[_MXType] = 10
mxUINT16_CLASS: Final[_MXType] = 11
mxINT32_CLASS: Final[_MXType] = 12
mxUINT32_CLASS: Final[_MXType] = 13
mxINT64_CLASS: Final[_MXType] = 14
mxUINT64_CLASS: Final[_MXType] = 15
mxFUNCTION_CLASS: Final[_MXType] = 16
mxOPAQUE_CLASS: Final[_MXType] = 17
mxOBJECT_CLASS_FROM_MATRIX_H: Final[_MXType] = 18

@final
class mat_struct: ...

class MatlabObject(np.ndarray[_ShapeT, np.dtype[np.void]], Generic[_ShapeT]):
    classname: Final[str | None]

    def __new__(cls, input_array: onp.AnyVoidArray, classname: str | None = None) -> Self: ...
    @override
    def __array_finalize__(self, /, obj: onp.ArrayND[np.void] | None) -> None: ...

class MatlabFunction(np.ndarray[_ShapeT, np.dtype[np.void]], Generic[_ShapeT]):
    @override
    def __new__(cls, input_array: onp.AnyVoidArray) -> Self: ...

class MatlabOpaque(np.ndarray[_ShapeT, np.dtype[np.void]], Generic[_ShapeT]):
    @override
    def __new__(cls, input_array: onp.AnyVoidArray) -> Self: ...

def _convert_codecs(template: _CodecsTemplate, byte_order: _ByteOrder) -> dict[_MCodec, _CodecBO]: ...
