import struct
from _typeshed import SupportsWrite
from collections.abc import Collection, Mapping
from typing import Any, Final

from .fragment import FragmentFD

u8: Final[struct.Struct]
u88: Final[struct.Struct]
u16: Final[struct.Struct]
u1616: Final[struct.Struct]
u32: Final[struct.Struct]
u64: Final[struct.Struct]
s88: Final[struct.Struct]
s16: Final[struct.Struct]
s1616: Final[struct.Struct]
s32: Final[struct.Struct]
unity_matrix: Final[bytes]
TRACK_ENABLED: Final = 0x1
TRACK_IN_MOVIE: Final = 0x2
TRACK_IN_PREVIEW: Final = 0x4
SELF_CONTAINED: Final = 0x1

def box(box_type: bytes, payload: bytes) -> bytes: ...
def full_box(box_type: bytes, version: int, flags: int, payload: bytes) -> bytes: ...
def write_piff_header(stream: SupportsWrite[bytes], params: Mapping[str, Any]) -> None: ...
def extract_box_data(data: bytes, box_sequence: Collection[bytes]) -> bytes | None: ...

class IsmFD(FragmentFD): ...
