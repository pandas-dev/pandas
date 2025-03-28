from enum import Enum
from typing import List


class FileSystemType(str, Enum):
    WINDOWS = "WINDOWS"
    LUSTRE = "LUSTRE"
    ONTAP = "ONTAP"
    OPEN_ZFS = "OPENZFS"

    @classmethod
    def list_values(self) -> List[str]:
        return sorted([item.value for item in FileSystemType])
