from _typeshed import StrOrBytesPath
from abc import abstractmethod
from logging import Logger
from subprocess import Popen
from typing import Any, ClassVar, Literal, NewType

import setuptools.command.sdist
from fanstatic.core import Resource

logger: Logger

class CompilerError(Exception): ...

class Compiler:
    @property
    @abstractmethod
    def name(self) -> str: ...
    @property
    @abstractmethod
    def source_extension(self) -> str: ...
    def __call__(self, resource: Resource, force: bool = False) -> None: ...
    def process(self, source: StrOrBytesPath, target: StrOrBytesPath) -> Any: ...
    def should_process(self, source: StrOrBytesPath, target: StrOrBytesPath) -> bool: ...
    @property
    @abstractmethod
    def available(self) -> bool: ...
    def source_path(self, resource: Resource) -> str | None: ...
    def target_path(self, resource: Resource) -> str | None: ...

class Minifier(Compiler):
    @property
    @abstractmethod
    def name(self) -> str: ...
    @property
    @abstractmethod
    def source_extension(self) -> str: ...
    @property
    @abstractmethod
    def target_extension(self) -> str: ...
    def source_to_target(self, resource: Resource) -> str: ...

def compile_resources(argv: list[str] = ...) -> None: ...

class sdist_compile(setuptools.command.sdist.sdist): ...

class NullCompiler(Compiler):
    name: ClassVar[Literal[""]]
    source_extension = NotImplemented
    def source_path(self, resource: Resource) -> None: ...
    def target_path(self, resource: Resource) -> None: ...
    def should_process(self, source: StrOrBytesPath, target: StrOrBytesPath) -> Literal[False]: ...
    @property
    def available(self) -> Literal[False]: ...

_SourceType = NewType("_SourceType", object)
_TargetType = NewType("_TargetType", object)
SOURCE: _SourceType
TARGET: _TargetType

class CommandlineBase:
    @property
    @abstractmethod
    def command(self) -> str: ...
    arguments: ClassVar[list[str]]
    @property
    def available(self) -> bool: ...
    def process(self, source: StrOrBytesPath | _SourceType, target: StrOrBytesPath | _TargetType) -> Popen[str]: ...

class CoffeeScript(CommandlineBase, Compiler):
    name: ClassVar[Literal["coffee"]]
    command: ClassVar[Literal["coffee"]]
    source_extension = NotImplemented
    def process(  # type: ignore[override]
        self, source: StrOrBytesPath | _SourceType, target: StrOrBytesPath | _TargetType
    ) -> None: ...

COFFEE_COMPILER: CoffeeScript

class LESS(CommandlineBase, Compiler):
    name: ClassVar[Literal["less"]]
    command: ClassVar[Literal["lessc"]]
    source_extension = NotImplemented

LESS_COMPILER: LESS

class SASS(CommandlineBase, Compiler):
    name: ClassVar[Literal["sass"]]
    command: ClassVar[Literal["sass"]]
    source_extension: ClassVar[Literal[".scss"]]

SASS_COMPILER: SASS

class PythonPackageBase:
    @property
    @abstractmethod
    def package(self) -> str: ...
    @property
    def available(self) -> bool: ...

class CSSMin(PythonPackageBase, Minifier):
    name: ClassVar[Literal["cssmin"]]
    package: ClassVar[Literal["cssmin"]]
    source_extension = NotImplemented
    target_extension: ClassVar[Literal[".min.css"]]

CSSMIN_MINIFIER: CSSMin

class JSMin(PythonPackageBase, Minifier):
    name: ClassVar[Literal["jsmin"]]
    package: ClassVar[Literal["jsmin"]]
    source_extension = NotImplemented
    target_extension: ClassVar[Literal[".min.js"]]

JSMIN_MINIFIER: JSMin

class Closure(PythonPackageBase, Minifier):
    name: ClassVar[Literal["closure"]]
    package: ClassVar[Literal["closure"]]
    source_extension = NotImplemented
    target_extension: Literal[".min.js"]
    arguments: ClassVar[list[str]]
    def process(self, source: StrOrBytesPath, target: StrOrBytesPath) -> Popen[str]: ...

CLOSURE_MINIFIER: Closure
