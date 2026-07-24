import sys
from typing import Final

VER: Final[sys._version_info]

def safe_string(value: object) -> str: ...
