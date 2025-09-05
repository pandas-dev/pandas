from __future__ import annotations

import struct

LFANEW_OFFSET = 30 * 2
FILE_HEADER_SIZE = 5 * 4
SUBSYSTEM_OFFSET = 17 * 4
IMAGE_SUBSYSTEM_WINDOWS_GUI = 2
IMAGE_SUBSYSTEM_WINDOWS_CUI = 3


def get_pe_subsystem(filename):
    with open(filename, "rb") as f:
        if f.read(2) != b"MZ":
            return None
        f.seek(LFANEW_OFFSET)
        lfanew = struct.unpack("L", f.read(4))[0]
        f.seek(lfanew)
        if f.read(4) != b"PE\x00\x00":
            return None
        f.seek(FILE_HEADER_SIZE + SUBSYSTEM_OFFSET, 1)
        return struct.unpack("H", f.read(2))[0]


# print(get_pe_subsystem("c:\\windows\\notepad.exe")) == 2
# print(get_pe_subsystem("c:\\python32\\python.exe")) == 3
# print(get_pe_subsystem("c:\\python32\\pythonw.exe")) == 2
