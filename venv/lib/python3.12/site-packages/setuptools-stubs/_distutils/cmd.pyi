from _typeshed import BytesPath, StrOrBytesPath, StrPath, Unused
from abc import abstractmethod
from collections.abc import Callable, MutableSequence
from typing import Any, ClassVar, TypeVar, overload
from typing_extensions import TypeVarTuple, Unpack

from .dist import Distribution

_StrPathT = TypeVar("_StrPathT", bound=StrPath)
_BytesPathT = TypeVar("_BytesPathT", bound=BytesPath)
_CommandT = TypeVar("_CommandT", bound=Command)
_Ts = TypeVarTuple("_Ts")

class Command:
    dry_run: bool  # Exposed from __getattr_. Same as Distribution.dry_run
    distribution: Distribution
    # Any to work around variance issues
    sub_commands: ClassVar[list[tuple[str, Callable[[Any], bool] | None]]]
    user_options: ClassVar[
        # Specifying both because list is invariant. Avoids mypy override assignment issues
        list[tuple[str, str, str]]
        | list[tuple[str, str | None, str]]
    ]
    def __init__(self, dist: Distribution) -> None: ...
    def ensure_finalized(self) -> None: ...
    @abstractmethod
    def initialize_options(self) -> None: ...
    @abstractmethod
    def finalize_options(self) -> None: ...
    @abstractmethod
    def run(self) -> None: ...
    def announce(self, msg: str, level: int = 10) -> None: ...
    def debug_print(self, msg: str) -> None: ...
    def ensure_string(self, option: str, default: str | None = None) -> None: ...
    def ensure_string_list(self, option: str) -> None: ...
    def ensure_filename(self, option: str) -> None: ...
    def ensure_dirname(self, option: str) -> None: ...
    def get_command_name(self) -> str: ...
    def set_undefined_options(self, src_cmd: str, *option_pairs: tuple[str, str]) -> None: ...
    # NOTE: Because this is private setuptools implementation and we don't re-expose all commands here,
    # we're not overloading each and every command possibility.
    def get_finalized_command(self, command: str, create: bool = True) -> Command: ...
    @overload
    def reinitialize_command(self, command: str, reinit_subcommands: bool = False) -> Command: ...
    @overload
    def reinitialize_command(self, command: _CommandT, reinit_subcommands: bool = False) -> _CommandT: ...
    def run_command(self, command: str) -> None: ...
    def get_sub_commands(self) -> list[str]: ...
    def warn(self, msg: str) -> None: ...
    def execute(
        self, func: Callable[[Unpack[_Ts]], Unused], args: tuple[Unpack[_Ts]], msg: str | None = None, level: int = 1
    ) -> None: ...
    def mkpath(self, name: str, mode: int = 0o777) -> None: ...
    @overload
    def copy_file(
        self,
        infile: StrPath,
        outfile: _StrPathT,
        preserve_mode: bool = True,
        preserve_times: bool = True,
        link: str | None = None,
        level: Unused = 1,
    ) -> tuple[_StrPathT | str, bool]: ...
    @overload
    def copy_file(
        self,
        infile: BytesPath,
        outfile: _BytesPathT,
        preserve_mode: bool = True,
        preserve_times: bool = True,
        link: str | None = None,
        level: Unused = 1,
    ) -> tuple[_BytesPathT | bytes, bool]: ...
    def copy_tree(
        self,
        infile: StrPath,
        outfile: str,
        preserve_mode: bool = True,
        preserve_times: bool = True,
        preserve_symlinks: bool = False,
        level: Unused = 1,
    ) -> list[str]: ...
    @overload
    def move_file(self, src: StrPath, dst: _StrPathT, level: Unused = 1) -> _StrPathT | str: ...
    @overload
    def move_file(self, src: BytesPath, dst: _BytesPathT, level: Unused = 1) -> _BytesPathT | bytes: ...
    def spawn(self, cmd: MutableSequence[str], search_path: bool = True, level: Unused = 1) -> None: ...
    @overload
    def make_archive(
        self,
        base_name: str,
        format: str,
        root_dir: StrOrBytesPath | None = None,
        base_dir: str | None = None,
        owner: str | None = None,
        group: str | None = None,
    ) -> str: ...
    @overload
    def make_archive(
        self,
        base_name: StrPath,
        format: str,
        root_dir: StrOrBytesPath,
        base_dir: str | None = None,
        owner: str | None = None,
        group: str | None = None,
    ) -> str: ...
    def make_file(
        self,
        infiles: str | list[str] | tuple[str, ...],
        outfile: StrOrBytesPath,
        func: Callable[[Unpack[_Ts]], Unused],
        args: tuple[Unpack[_Ts]],
        exec_msg: str | None = None,
        skip_msg: str | None = None,
        level: Unused = 1,
    ) -> None: ...
