# PYZ, EXE and COLLECT referenced in https://pyinstaller.org/en/stable/spec-files.html#spec-file-operation
# MERGE is referenced in https://pyinstaller.org/en/stable/spec-files.html#example-merge-spec-file
# hide_console referenced in https://pyinstaller.org/en/stable/feature-notes.html#automatic-hiding-and-minimization-of-console-window-under-windows
# Not to be imported during runtime, but is the type reference for spec files which are executed as python code
import sys
from _typeshed import FileDescriptorOrPath, StrOrBytesPath, StrPath, Unused
from collections.abc import Iterable, Mapping, Sequence
from types import CodeType
from typing import ClassVar, Final, Literal
from typing_extensions import TypeAlias

from PyInstaller.building import _PyiBlockCipher
from PyInstaller.building.build_main import Analysis
from PyInstaller.building.datastruct import Target, _TOCTuple
from PyInstaller.building.splash import Splash
from PyInstaller.utils.win32.versioninfo import VSVersionInfo

if sys.platform == "darwin":
    _TargetArch: TypeAlias = Literal["x86_64", "arm64", "universal2"]
    _SupportedTargetArchParam: TypeAlias = _TargetArch | None
    _CodesignIdentity: TypeAlias = str | None
    _CodesignIdentityParam: TypeAlias = str | None
else:
    _TargetArch: TypeAlias = None
    _SupportedTargetArchParam: TypeAlias = Unused
    _CodesignIdentity: TypeAlias = None
    _CodesignIdentityParam: TypeAlias = Unused

if sys.platform == "win32":
    _Icon: TypeAlias = list[StrPath] | str
    _IconParam: TypeAlias = StrPath | list[StrPath] | None
elif sys.platform == "darwin":
    _Icon: TypeAlias = list[StrPath] | None
    _IconParam: TypeAlias = StrPath | list[StrPath] | None
else:
    _Icon: TypeAlias = None
    _IconParam: TypeAlias = Unused

if sys.platform == "win32":
    _VersionSrc: TypeAlias = VSVersionInfo | None
    _VersionParam: TypeAlias = VSVersionInfo | StrOrBytesPath | None
    _Manifest: TypeAlias = bytes
    _ManifestParam: TypeAlias = str | None
else:
    _VersionSrc: TypeAlias = None
    _VersionParam: TypeAlias = Unused
    _Manifest: TypeAlias = None
    _ManifestParam: TypeAlias = Unused

_HideConsole: TypeAlias = Literal["hide-early", "minimize-early", "hide-late", "minimize-late"] | None

class PYZ(Target):
    name: str
    cipher: _PyiBlockCipher
    dependencies: list[_TOCTuple]
    toc: list[_TOCTuple]
    code_dict: dict[str, CodeType]
    def __init__(self, *tocs: Iterable[_TOCTuple], name: str | None = None, cipher: _PyiBlockCipher = None) -> None: ...
    def assemble(self) -> None: ...

class PKG(Target):
    xformdict: ClassVar[dict[str, str]]
    toc: list[_TOCTuple]
    cdict: Mapping[str, bool]
    python_lib_name: str
    name: str
    exclude_binaries: bool
    strip_binaries: bool
    upx_binaries: bool
    upx_exclude: Iterable[str]
    target_arch: _TargetArch | None
    codesign_identity: _CodesignIdentity
    entitlements_file: FileDescriptorOrPath | None
    def __init__(
        self,
        toc: Iterable[_TOCTuple],
        python_lib_name: str,
        name: str | None = None,
        cdict: Mapping[str, bool] | None = None,
        exclude_binaries: bool = False,
        strip_binaries: bool = False,
        upx_binaries: bool = False,
        upx_exclude: Iterable[str] | None = None,
        target_arch: _SupportedTargetArchParam = None,
        codesign_identity: _CodesignIdentityParam = None,
        entitlements_file: FileDescriptorOrPath | None = None,
    ) -> None: ...
    def assemble(self) -> None: ...

class EXE(Target):
    exclude_binaries: bool
    bootloader_ignore_signals: bool
    console: bool
    hide_console: _HideConsole
    disable_windowed_traceback: bool
    debug: bool
    name: str
    icon: _Icon
    versrsrc: _VersionSrc
    manifest: _Manifest
    embed_manifest: bool
    resources: Sequence[str]
    strip: bool
    upx_exclude: Iterable[str]
    runtime_tmpdir: str | None
    contents_directory: str | None
    append_pkg: bool
    uac_admin: bool
    uac_uiaccess: bool
    argv_emulation: bool
    target_arch: _TargetArch
    codesign_identity: _CodesignIdentity
    entitlements_file: FileDescriptorOrPath | None
    upx: bool
    pkgname: str
    toc: list[_TOCTuple]
    pkg: PKG
    dependencies: list[_TOCTuple]
    exefiles: list[_TOCTuple]
    def __init__(
        self,
        *args: Iterable[_TOCTuple] | PYZ | Splash,
        exclude_binaries: bool = False,
        bootloader_ignore_signals: bool = False,
        console: bool = True,
        hide_console: _HideConsole = None,
        disable_windowed_traceback: bool = False,
        debug: bool = False,
        name: str | None = None,
        icon: _IconParam = None,
        version: _VersionParam = None,
        manifest: _ManifestParam = None,
        embed_manifest: Literal[True] = True,
        resources: Sequence[str] = ...,
        strip: bool = False,
        upx_exclude: Iterable[str] = ...,
        runtime_tmpdir: str | None = None,
        contents_directory: str = "_internal",
        append_pkg: bool = True,
        uac_admin: bool = False,
        uac_uiaccess: bool = False,
        argv_emulation: bool = False,
        target_arch: _SupportedTargetArchParam = None,
        codesign_identity: _CodesignIdentityParam = None,
        entitlements_file: FileDescriptorOrPath | None = None,
        upx: bool = False,
        cdict: Mapping[str, bool] | None = None,
    ) -> None: ...
    mtm: float
    def assemble(self) -> None: ...

class COLLECT(Target):
    strip_binaries: bool
    upx_exclude: Iterable[str]
    console: bool
    target_arch: _TargetArch | None
    codesign_identity: _CodesignIdentity
    entitlements_file: FileDescriptorOrPath | None
    upx_binaries: bool
    name: str
    toc: list[_TOCTuple]
    def __init__(
        self,
        *args: Iterable[_TOCTuple] | EXE,
        strip: bool = False,
        upx_exclude: Iterable[str] = ...,
        upx: bool = False,
        name: str,
    ) -> None: ...
    def assemble(self) -> None: ...

class MERGE:
    def __init__(self, *args: tuple[Analysis, Unused, str]) -> None: ...

UNCOMPRESSED: Final = False
COMPRESSED: Final = True
