from __future__ import annotations

__lazy_modules__ = {"struct"}

import struct

LFANEW_OFFSET = 30 * 2
FILE_HEADER_SIZE = 5 * 4
SUBSYSTEM_OFFSET = 17 * 4
IMAGE_SUBSYSTEM_WINDOWS_GUI = 2
IMAGE_SUBSYSTEM_WINDOWS_CUI = 3


def get_pe_subsystem(filename: str) -> int | None:
    with open(filename, "rb") as f:
        if f.read(2) != b"MZ":
            return None
        f.seek(LFANEW_OFFSET)
        # Explicit little-endian widths: the PE fields are fixed 32-/16-bit
        # values, whereas native "L"/"H" sizes are platform-dependent.
        lfanew = struct.unpack("<I", f.read(4))[0]
        f.seek(lfanew)
        if f.read(4) != b"PE\x00\x00":
            return None
        f.seek(FILE_HEADER_SIZE + SUBSYSTEM_OFFSET, 1)
        return int(struct.unpack("<H", f.read(2))[0])


# print(get_pe_subsystem("c:\\windows\\notepad.exe")) == 2
# print(get_pe_subsystem("c:\\python32\\python.exe")) == 3
# print(get_pe_subsystem("c:\\python32\\pythonw.exe")) == 2


__all__ = [
    "FILE_HEADER_SIZE",
    "IMAGE_SUBSYSTEM_WINDOWS_CUI",
    "IMAGE_SUBSYSTEM_WINDOWS_GUI",
    "LFANEW_OFFSET",
    "SUBSYSTEM_OFFSET",
    "get_pe_subsystem",
]


def __dir__() -> list[str]:
    return list(__all__)
