# https://pyinstaller.org/en/stable/hooks.html#module-PyInstaller.compat
from _typeshed import FileDescriptorOrPath
from collections.abc import Iterable
from types import ModuleType
from typing import Final, Literal, overload

strict_collect_mode: bool
is_64bits: Final[bool]
is_py35: Final = True
is_py36: Final = True
is_py37: Final = True
is_py38: Final = True
is_py39: Final[bool]
is_py310: Final[bool]
is_py311: Final[bool]
is_py312: Final[bool]
is_py313: Final[bool]
is_py314: Final[bool]
is_win: Final[bool]
is_win_10: Final[bool]
is_win_11: Final[bool]
is_win_wine: Final[bool]
is_cygwin: Final[bool]
is_darwin: Final[bool]
is_linux: Final[bool]
is_solar: Final[bool]
is_aix: Final[bool]
is_freebsd: Final[bool]
is_openbsd: Final[bool]
is_hpux: Final[bool]
is_unix: Final[bool]
is_musl: Final[bool]
is_termux: Final[bool]
is_macos_11_compat: Final[bool]
is_macos_11_native: Final[bool]
is_macos_11: Final[bool]
is_nogil: Final[bool]
base_prefix: Final[str]
is_venv: Final[bool]
is_virtualenv: Final[bool]
is_conda: Final[bool]
is_pure_conda: Final[bool]
python_executable: Final[str]
is_ms_app_store: Final[bool]
BYTECODE_MAGIC: Final[bytes]
EXTENSION_SUFFIXES: Final[list[str]]
ALL_SUFFIXES: Final[list[str]]

architecture: Final[Literal["64bit", "n32bit", "32bit"]]
system: Final[Literal["Cygwin", "Linux", "Darwin", "Java", "Windows"]]
machine: Final[
    Literal["AMD64", "x86", "ARM64", "sw_64", "loongarch64", "arm", "intel", "ppc", "mips", "riscv", "s390x", "unknown"] | None
]

def is_wine_dll(filename: FileDescriptorOrPath) -> bool: ...
@overload
def getenv(name: str, default: str) -> str: ...
@overload
def getenv(name: str, default: None = None) -> str | None: ...
def setenv(name: str, value: str) -> None: ...
def unsetenv(name: str) -> None: ...
def exec_command(
    *cmdargs: str, encoding: str | None = None, raise_enoent: bool | None = None, **kwargs: int | bool | Iterable[int] | None
) -> str: ...
def exec_command_rc(*cmdargs: str, **kwargs: float | bool | Iterable[int] | None) -> int: ...
def exec_command_all(
    *cmdargs: str, encoding: str | None = None, **kwargs: int | bool | Iterable[int] | None
) -> tuple[int, str, str]: ...
def exec_python(*args: str, **kwargs: str | None) -> str: ...
def exec_python_rc(*args: str, **kwargs: str | None) -> int: ...
def getsitepackages(prefixes: Iterable[str] | None = None) -> list[str]: ...
def importlib_load_source(name: str, pathname: str) -> ModuleType: ...

PY3_BASE_MODULES: Final[set[str]]
PURE_PYTHON_MODULE_TYPES: Final[set[str]]
SPECIAL_MODULE_TYPES: Final[set[str]]
BINARY_MODULE_TYPES: Final[set[str]]
VALID_MODULE_TYPES: Final[set[str]]
BAD_MODULE_TYPES: Final[set[str]]
ALL_MODULE_TYPES: Final[set[str]]
MODULE_TYPES_TO_TOC_DICT: Final[dict[str, str]]

def check_requirements() -> None: ...
