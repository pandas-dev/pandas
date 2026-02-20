from typing import Literal, TypeAlias

__all__ = [
    "ByteOrder",
    "ByteOrderChar",
    "ByteOrderName",
    "Casting",
    "CastingSafe",
    "CastingUnsafe",
    "ConvolveMode",
    "Device",
    "OrderACF",
    "OrderCF",
    "OrderKACF",
    "PartitionKind",
    "SortKind",
    "SortSide",
]

Device: TypeAlias = Literal["cpu"]

ByteOrderChar: TypeAlias = Literal["<", ">", "=", "|"]
ByteOrderName: TypeAlias = Literal["little", "big", "native", "ignore", "swap"]
_ByteOrderShort: TypeAlias = Literal["L", "B", "N", "I", "S"]
ByteOrder: TypeAlias = Literal[ByteOrderChar, _ByteOrderShort, ByteOrderName]

CastingSafe: TypeAlias = Literal["no", "equiv", "safe", "same_kind"]
CastingUnsafe: TypeAlias = Literal["unsafe"]
Casting: TypeAlias = Literal[CastingSafe, CastingUnsafe]

OrderCF: TypeAlias = Literal["C", "F"]
OrderACF: TypeAlias = Literal["A", OrderCF]
OrderKACF: TypeAlias = Literal["K", OrderACF]

IndexMode: TypeAlias = Literal["raise", "wrap", "clip"]

PartitionKind: TypeAlias = Literal["introselect"]

SortKind: TypeAlias = Literal[
    "Q", "quick", "quicksort",
    "M", "merge", "mergesort",
    "H", "heap", "heapsort",
    "S", "stable", "stablesort",
]  # fmt: skip
SortSide: TypeAlias = Literal["left", "right"]

ConvolveMode: TypeAlias = Literal["full", "same", "valid"]
