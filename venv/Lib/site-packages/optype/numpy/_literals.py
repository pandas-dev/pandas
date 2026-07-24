from typing import Literal

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

type Device = Literal["cpu"]

type ByteOrderChar = Literal["<", ">", "=", "|"]
type ByteOrderName = Literal["little", "big", "native", "ignore", "swap"]
type _ByteOrderShort = Literal["L", "B", "N", "I", "S"]
type ByteOrder = Literal[ByteOrderChar, _ByteOrderShort, ByteOrderName]

type CastingSafe = Literal["no", "equiv", "safe", "same_kind"]
type CastingUnsafe = Literal["unsafe"]
type Casting = Literal[CastingSafe, CastingUnsafe]

type OrderCF = Literal["C", "F"]
type OrderACF = Literal["A", OrderCF]
type OrderKACF = Literal["K", OrderACF]

type IndexMode = Literal["raise", "wrap", "clip"]

type PartitionKind = Literal["introselect"]

type SortKind = Literal[
    "Q", "quick", "quicksort",
    "M", "merge", "mergesort",
    "H", "heap", "heapsort",
    "S", "stable", "stablesort",
]  # fmt: skip
type SortSide = Literal["left", "right"]

type ConvolveMode = Literal["full", "same", "valid"]
