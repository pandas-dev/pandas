from _typeshed import SupportsKeysAndGetItem
from collections.abc import Iterable
from typing_extensions import TypeAlias

# Used to update PyInstaller.config.CONF
_PyIConfig: TypeAlias = (
    SupportsKeysAndGetItem[str, bool | str | list[str] | None] | Iterable[tuple[str, bool | str | list[str] | None]]
)

# https://pyinstaller.org/en/stable/usage.html#running-pyinstaller-from-python-code
def run(pyi_args: Iterable[str] | None = None, pyi_config: _PyIConfig | None = None) -> None: ...
def check_unsafe_privileges() -> None: ...
