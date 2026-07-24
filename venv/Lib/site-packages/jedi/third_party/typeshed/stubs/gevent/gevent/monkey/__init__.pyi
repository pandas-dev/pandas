from typing import Any

from gevent.monkey.api import get_original as get_original, patch_module as patch_module

class MonkeyPatchWarning(RuntimeWarning): ...

def is_module_patched(mod_name: str) -> bool: ...
def is_object_patched(mod_name: str, item_name: str) -> bool: ...
def patch_os() -> None: ...
def patch_queue() -> None: ...
def patch_time() -> None: ...
def patch_thread(
    threading: bool = True, _threading_local: bool = True, Event: bool = True, logging: bool = True, existing_locks: bool = True
) -> None: ...
def patch_socket(dns: bool = True, aggressive: bool = True) -> None: ...
def patch_builtins() -> None: ...
def patch_dns() -> None: ...
def patch_ssl() -> None: ...
def patch_select(aggressive: bool = True) -> None: ...
def patch_selectors(aggressive: bool = True) -> None: ...
def patch_subprocess() -> None: ...
def patch_sys(stdin: bool = True, stdout: bool = True, stderr: bool = True) -> None: ...
def patch_signal() -> None: ...
def patch_all(
    socket: bool = True,
    dns: bool = True,
    time: bool = True,
    select: bool = True,
    thread: bool = True,
    os: bool = True,
    ssl: bool = True,
    subprocess: bool = True,
    sys: bool = False,
    aggressive: bool = True,
    Event: bool = True,
    builtins: bool = True,  # does nothing on Python 3
    signal: bool = True,
    queue: bool = True,
    contextvars: bool = True,  # does nothing on Python 3.7+
    **kwargs: object,
) -> bool | None: ...
def main() -> dict[str, Any]: ...

__all__ = [
    "patch_all",
    "patch_builtins",
    "patch_dns",
    "patch_os",
    "patch_queue",
    "patch_select",
    "patch_signal",
    "patch_socket",
    "patch_ssl",
    "patch_subprocess",
    "patch_sys",
    "patch_thread",
    "patch_time",
    # query functions
    "get_original",
    "is_module_patched",
    "is_object_patched",
    # plugin API
    "patch_module",
    # module functions
    "main",
    # Errors and warnings
    "MonkeyPatchWarning",
]
