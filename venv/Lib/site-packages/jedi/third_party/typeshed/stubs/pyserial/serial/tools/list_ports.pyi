import re
import sys
from collections.abc import Generator

from serial.tools.list_ports_common import ListPortInfo

if sys.platform == "win32":
    from serial.tools.list_ports_windows import comports as comports
else:
    from serial.tools.list_ports_posix import comports as comports

def grep(regexp: str | re.Pattern[str], include_links: bool = False) -> Generator[ListPortInfo]: ...
def main() -> None: ...
