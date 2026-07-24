import sys

from serial.tools.list_ports_common import ListPortInfo

if sys.platform != "win32":
    if sys.platform == "linux":
        from serial.tools.list_ports_linux import comports as comports
    elif sys.platform == "darwin":
        from serial.tools.list_ports_osx import comports as comports
    else:
        def comports(include_links: bool = ...) -> list[ListPortInfo]: ...
