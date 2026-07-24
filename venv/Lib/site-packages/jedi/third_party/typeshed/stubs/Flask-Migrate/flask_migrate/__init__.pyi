# pyright: reportInvalidStubStatement=none

import sys
from _typeshed import StrPath, SupportsFlush, SupportsKeysAndGetItem, SupportsWrite
from argparse import Namespace
from collections.abc import Callable, Iterable, Sequence
from logging import Logger
from typing import Any, Protocol, TypeVar, type_check_only
from typing_extensions import ParamSpec, TypeAlias

import flask
from flask_sqlalchemy import SQLAlchemy

_T = TypeVar("_T")
_T_contra = TypeVar("_T_contra", contravariant=True)
_P = ParamSpec("_P")
_ConfigureCallback: TypeAlias = Callable[[Config], Config]
_AlembicConfigValue: TypeAlias = Any

alembic_version: tuple[int, int, int]
log: Logger

@type_check_only
class _SupportsWriteAndFlush(SupportsWrite[_T_contra], SupportsFlush, Protocol): ...

class Config:  # should inherit from alembic.config.Config which is not possible yet
    template_directory: str | None
    # Same as alembic.config.Config + template_directory kwarg
    def __init__(
        self,
        file_: StrPath | None = None,
        ini_section: str = "alembic",
        # Same as buffer argument in TextIOWrapper.__init__.buffer
        output_buffer: _SupportsWriteAndFlush[str] | None = None,
        # Same as stream argument in alembic.util.messaging
        stdout: SupportsWrite[str] = sys.stdout,
        cmd_opts: Namespace | None = None,
        config_args: SupportsKeysAndGetItem[str, _AlembicConfigValue] | Iterable[tuple[str, _AlembicConfigValue]] = ...,
        attributes: (
            SupportsKeysAndGetItem[_AlembicConfigValue, _AlembicConfigValue]
            | Iterable[tuple[_AlembicConfigValue, _AlembicConfigValue]]
            | None
        ) = None,
        *,
        template_directory: str | None = None,
    ) -> None: ...
    def get_template_directory(self) -> str: ...

class Migrate:
    configure_callbacks: list[_ConfigureCallback]
    db: SQLAlchemy | None
    directory: str
    alembic_ctx_kwargs: dict[str, _AlembicConfigValue]
    def __init__(
        self,
        app: flask.Flask | None = None,
        db: SQLAlchemy | None = None,
        directory: str = "migrations",
        command: str = "db",
        compare_type: bool = True,
        render_as_batch: bool = True,
        **kwargs: _AlembicConfigValue,
    ) -> None: ...
    def init_app(
        self,
        app: flask.Flask,
        db: SQLAlchemy | None = None,
        directory: str | None = None,
        command: str | None = None,
        compare_type: bool | None = None,
        render_as_batch: bool | None = None,
        **kwargs: _AlembicConfigValue,
    ) -> None: ...
    def configure(self, f: _ConfigureCallback) -> _ConfigureCallback: ...
    def call_configure_callbacks(self, config: Config) -> Config: ...
    def get_config(
        self, directory: str | None = None, x_arg: str | Sequence[str] | None = None, opts: Iterable[str] | None = None
    ) -> Config: ...

def catch_errors(f: Callable[_P, _T]) -> Callable[_P, _T]: ...
def list_templates() -> None: ...
def init(directory: str | None = None, multidb: bool = False, template: str | None = None, package: bool = False) -> None: ...
def revision(
    directory: str | None = None,
    message: str | None = None,
    autogenerate: bool = False,
    sql: bool = False,
    head: str = "head",
    splice: bool = False,
    branch_label: str | None = None,
    version_path: str | None = None,
    rev_id: str | None = None,
) -> None: ...
def migrate(
    directory: str | None = None,
    message: str | None = None,
    sql: bool = False,
    head: str = "head",
    splice: bool = False,
    branch_label: str | None = None,
    version_path: str | None = None,
    rev_id: str | None = None,
    x_arg: str | Sequence[str] | None = None,
) -> None: ...
def edit(directory: str | None = None, revision: str = "current") -> None: ...
def merge(
    directory: str | None = None,
    revisions: str = "",
    message: str | None = None,
    branch_label: str | None = None,
    rev_id: str | None = None,
) -> None: ...
def upgrade(
    directory: str | None = None,
    revision: str = "head",
    sql: bool = False,
    tag: str | None = None,
    x_arg: str | Sequence[str] | None = None,
) -> None: ...
def downgrade(
    directory: str | None = None,
    revision: str = "-1",
    sql: bool = False,
    tag: str | None = None,
    x_arg: str | Sequence[str] | None = None,
) -> None: ...
def show(directory: str | None = None, revision: str = "head") -> None: ...
def history(
    directory: str | None = None, rev_range: str | None = None, verbose: bool = False, indicate_current: bool = False
) -> None: ...
def heads(directory: str | None = None, verbose: bool = False, resolve_dependencies: bool = False) -> None: ...
def branches(directory: str | None = None, verbose: bool = False) -> None: ...
def current(directory: str | None = None, verbose: bool = False) -> None: ...
def stamp(
    directory: str | None = None, revision: str = "head", sql: bool = False, tag: str | None = None, purge: bool = False
) -> None: ...
def check(directory: str | None = None) -> None: ...
