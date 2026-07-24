from serial.tools.list_ports_common import ListPortInfo

class SysFS(ListPortInfo):
    usb_device_path: str | None
    device_path: str | None
    subsystem: str | None
    usb_interface_path: str | None
    def __init__(self, device: str) -> None: ...
    def read_line(self, *args: str) -> str | None: ...

def comports(include_links: bool = False) -> list[SysFS]: ...
