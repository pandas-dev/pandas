from _typeshed import Incomplete, StrPath
from collections.abc import Mapping
from contextlib import _GeneratorContextManager
from typing import NoReturn
from typing_extensions import TypeAlias

from . import dist

__all__ = [
    "get_requires_for_build_sdist",
    "get_requires_for_build_wheel",
    "prepare_metadata_for_build_wheel",
    "build_wheel",
    "build_sdist",
    "get_requires_for_build_editable",
    "prepare_metadata_for_build_editable",
    "build_editable",
    "__legacy__",
    "SetupRequirementsError",
]

_ConfigSettings: TypeAlias = Mapping[str, str | list[str] | None] | None

class SetupRequirementsError(BaseException):
    specifiers: Incomplete
    def __init__(self, specifiers) -> None: ...

class Distribution(dist.Distribution):
    def fetch_build_eggs(self, specifiers) -> NoReturn: ...
    @classmethod
    def patch(cls) -> _GeneratorContextManager[None]: ...

class _BuildMetaBackend:
    def run_setup(self, setup_script: str = "setup.py") -> None: ...
    def get_requires_for_build_wheel(self, config_settings: _ConfigSettings = None) -> list[str]: ...
    def get_requires_for_build_sdist(self, config_settings: _ConfigSettings = None) -> list[str]: ...
    def prepare_metadata_for_build_wheel(self, metadata_directory: StrPath, config_settings: _ConfigSettings = None) -> str: ...
    def build_wheel(
        self, wheel_directory: StrPath, config_settings: _ConfigSettings = None, metadata_directory: StrPath | None = None
    ) -> str: ...
    def build_sdist(self, sdist_directory: StrPath, config_settings: _ConfigSettings = None) -> str: ...
    def build_editable(
        self, wheel_directory: StrPath, config_settings: _ConfigSettings = None, metadata_directory: StrPath | None = None
    ) -> str: ...
    def get_requires_for_build_editable(self, config_settings: _ConfigSettings = None) -> list[str]: ...
    def prepare_metadata_for_build_editable(
        self, metadata_directory: StrPath, config_settings: _ConfigSettings = None
    ) -> str: ...

class _BuildMetaLegacyBackend(_BuildMetaBackend):
    def run_setup(self, setup_script: str = "setup.py") -> None: ...

_BACKEND: _BuildMetaBackend
get_requires_for_build_wheel = _BACKEND.get_requires_for_build_wheel
get_requires_for_build_sdist = _BACKEND.get_requires_for_build_sdist
prepare_metadata_for_build_wheel = _BACKEND.prepare_metadata_for_build_wheel
build_wheel = _BACKEND.build_wheel
build_sdist = _BACKEND.build_sdist

get_requires_for_build_editable = _BACKEND.get_requires_for_build_editable
prepare_metadata_for_build_editable = _BACKEND.prepare_metadata_for_build_editable
build_editable = _BACKEND.build_editable

__legacy__: _BuildMetaLegacyBackend
