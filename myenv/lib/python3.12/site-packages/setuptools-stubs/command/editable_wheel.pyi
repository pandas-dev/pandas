from _typeshed import Incomplete, StrPath
from collections.abc import Iterator
from enum import Enum
from pathlib import Path
from types import TracebackType
from typing import Protocol
from typing_extensions import Self, TypeAlias

from .. import Command, errors, namespaces
from ..dist import Distribution
from ..warnings import SetuptoolsWarning

# Actually from wheel.wheelfile import WheelFile
_WheelFile: TypeAlias = Incomplete

class _EditableMode(Enum):
    STRICT = "strict"
    LENIENT = "lenient"
    COMPAT = "compat"
    @classmethod
    def convert(cls, mode: str | None) -> _EditableMode: ...

class editable_wheel(Command):
    description: str
    user_options: Incomplete
    dist_dir: Incomplete
    dist_info_dir: Incomplete
    project_dir: Incomplete
    mode: Incomplete
    def initialize_options(self) -> None: ...
    package_dir: Incomplete
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...

class EditableStrategy(Protocol):
    def __call__(self, wheel: _WheelFile, files: list[str], mapping: dict[str, str]) -> None: ...
    def __enter__(self): ...
    def __exit__(
        self, _exc_type: type[BaseException] | None, _exc_value: BaseException | None, _traceback: TracebackType | None
    ) -> None: ...

class _StaticPth:
    dist: Incomplete
    name: Incomplete
    path_entries: Incomplete
    def __init__(self, dist: Distribution, name: str, path_entries: list[Path]) -> None: ...
    def __call__(self, wheel: _WheelFile, files: list[str], mapping: dict[str, str]): ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, _exc_type: type[BaseException] | None, _exc_value: BaseException | None, _traceback: TracebackType | None
    ) -> None: ...

class _LinkTree(_StaticPth):
    auxiliary_dir: Incomplete
    build_lib: Incomplete
    def __init__(self, dist: Distribution, name: str, auxiliary_dir: StrPath, build_lib: StrPath) -> None: ...
    def __call__(self, wheel: _WheelFile, files: list[str], mapping: dict[str, str]): ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, _exc_type: type[BaseException] | None, _exc_value: BaseException | None, _traceback: TracebackType | None
    ) -> None: ...

class _TopLevelFinder:
    dist: Incomplete
    name: Incomplete
    def __init__(self, dist: Distribution, name: str) -> None: ...
    def template_vars(self) -> tuple[str, str, dict[str, str], dict[str, list[str]]]: ...
    def get_implementation(self) -> Iterator[tuple[str, bytes]]: ...
    def __call__(self, wheel: _WheelFile, files: list[str], mapping: dict[str, str]): ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, _exc_type: type[BaseException] | None, _exc_value: BaseException | None, _traceback: TracebackType | None
    ) -> None: ...

class _NamespaceInstaller(namespaces.Installer):
    distribution: Incomplete
    src_root: Incomplete
    installation_dir: Incomplete
    editable_name: Incomplete
    outputs: Incomplete
    dry_run: bool
    def __init__(self, distribution, installation_dir, editable_name, src_root) -> None: ...

class InformationOnly(SetuptoolsWarning): ...
class LinksNotSupported(errors.FileError): ...
class _DebuggingTips(SetuptoolsWarning): ...
