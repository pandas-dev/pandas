from enum import Enum


class FileSystemType(str, Enum):
    WINDOWS = "WINDOWS"
    LUSTRE = "LUSTRE"
    ONTAP = "ONTAP"
    OPEN_ZFS = "OPENZFS"

    @classmethod
    def list_values(self) -> list[str]:
        return sorted([item.value for item in FileSystemType])
