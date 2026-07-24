from _typeshed import FileDescriptorOrPath
from typing import Final

from Xlib.X import (
    FamilyChaos as FamilyChaos,
    FamilyDECnet as FamilyDECnet,
    FamilyInternet as FamilyInternet,
    FamilyInternetV6 as FamilyInternetV6,
    FamilyServerInterpreted as FamilyServerInterpreted,
)

FamilyLocal: Final = 256

class Xauthority:
    entries: list[tuple[bytes, bytes, bytes, bytes, bytes]]
    def __init__(self, filename: FileDescriptorOrPath | None = None) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, i: int) -> tuple[bytes, bytes, bytes, bytes, bytes]: ...
    def get_best_auth(
        self, family: bytes, address: bytes, dispno: bytes, types: tuple[bytes, ...] = (b"MIT-MAGIC-COOKIE-1",)
    ) -> tuple[bytes, bytes]: ...
