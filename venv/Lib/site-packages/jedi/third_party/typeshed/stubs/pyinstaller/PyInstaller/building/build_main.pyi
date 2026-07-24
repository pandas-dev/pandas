from _typeshed import StrPath
from collections.abc import Iterable
from typing import Any, Literal

from PyInstaller.building import _PyiBlockCipher
from PyInstaller.building.datastruct import Target, _TOCTuple

# Referenced in: https://pyinstaller.org/en/stable/hooks.html#PyInstaller.utils.hooks.get_hook_config
# Not to be imported during runtime, but is the type reference for hooks and analysis configuration
# Also referenced in https://pyinstaller.org/en/stable/spec-files.html
# Not to be imported during runtime, but is the type reference for spec files which are executed as python code
class Analysis(Target):
    # https://pyinstaller.org/en/stable/hooks-config.html#hook-configuration-options
    hooksconfig: dict[str, dict[str, object]]
    # https://pyinstaller.org/en/stable/spec-files.html#spec-file-operation
    # https://pyinstaller.org/en/stable/feature-notes.html
    pure: list[_TOCTuple]
    zipped_data: list[_TOCTuple]
    # https://pyinstaller.org/en/stable/spec-files.html#giving-run-time-python-options
    # https://pyinstaller.org/en/stable/spec-files.html#the-splash-target
    scripts: list[_TOCTuple]
    # https://pyinstaller.org/en/stable/feature-notes.html#practical-examples
    binaries: list[_TOCTuple]
    zipfiles: list[_TOCTuple]
    datas: list[_TOCTuple]

    inputs: list[str]
    dependencies: list[_TOCTuple]
    noarchive: bool
    optimize: int
    pathex: list[StrPath]
    hiddenimports: list[str]
    hookspath: list[tuple[StrPath, int]]
    excludes: list[str]
    custom_runtime_hooks: list[StrPath]
    # https://pyinstaller.org/en/stable/hooks.html#hook-global-variables
    module_collection_mode: dict[str, str]
    def __init__(
        self,
        scripts: Iterable[StrPath],
        pathex: Iterable[StrPath] | None = None,
        binaries: Iterable[tuple[StrPath, StrPath]] | None = None,
        datas: Iterable[tuple[StrPath, StrPath]] | None = None,
        hiddenimports: Iterable[str] | None = None,
        hookspath: Iterable[StrPath] | None = None,
        hooksconfig: dict[str, dict[str, Any]] | None = None,
        excludes: Iterable[str] | None = None,
        runtime_hooks: Iterable[StrPath] | None = None,
        cipher: _PyiBlockCipher = None,
        win_no_prefer_redirects: bool = False,
        win_private_assemblies: bool = False,
        noarchive: bool = False,
        module_collection_mode: dict[str, str] | None = None,
        optimize: Literal[-1, 0, 1, 2] | None = -1,
    ) -> None: ...
