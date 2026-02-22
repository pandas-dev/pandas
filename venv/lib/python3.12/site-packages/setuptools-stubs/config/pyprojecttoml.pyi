from _typeshed import Incomplete, StrPath
from types import TracebackType
from typing import Any
from typing_extensions import Self

from ..dist import Distribution
from . import expand

def load_file(filepath: StrPath) -> dict[Incomplete, Incomplete]: ...
def validate(config: dict[Incomplete, Incomplete], filepath: StrPath) -> bool: ...
def apply_configuration(dist: Distribution, filepath: StrPath, ignore_option_errors: bool = False) -> Distribution: ...
def read_configuration(
    filepath: StrPath, expand: bool = True, ignore_option_errors: bool = False, dist: Distribution | None = None
) -> dict[str, Any]: ...
def expand_configuration(
    config: dict[Incomplete, Incomplete],
    root_dir: StrPath | None = None,
    ignore_option_errors: bool = False,
    dist: Distribution | None = None,
) -> dict[Incomplete, Incomplete]: ...

class _ConfigExpander:
    config: dict[Incomplete, Incomplete]
    root_dir: StrPath
    project_cfg: Incomplete
    dynamic: Incomplete
    setuptools_cfg: Incomplete
    dynamic_cfg: Incomplete
    ignore_option_errors: bool
    def __init__(
        self,
        config: dict[Incomplete, Incomplete],
        root_dir: StrPath | None = None,
        ignore_option_errors: bool = False,
        dist: Distribution | None = None,
    ) -> None: ...
    def expand(self): ...

class _EnsurePackagesDiscovered(expand.EnsurePackagesDiscovered):
    def __init__(
        self, distribution: Distribution, project_cfg: dict[Incomplete, Incomplete], setuptools_cfg: dict[Incomplete, Incomplete]
    ) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_value: BaseException | None, traceback: TracebackType | None
    ) -> None: ...
