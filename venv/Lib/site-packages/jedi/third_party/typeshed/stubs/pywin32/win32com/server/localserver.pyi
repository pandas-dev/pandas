from collections.abc import Iterable
from typing import Final

from _win32typing import PyIID

usage: Final[str]

def serve(clsids: Iterable[PyIID]) -> None: ...
def main() -> None: ...
