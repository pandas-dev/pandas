from typing_extensions import TypeAlias

from yaml.composer import Composer
from yaml.constructor import BaseConstructor, Constructor, FullConstructor, SafeConstructor
from yaml.parser import Parser
from yaml.reader import Reader
from yaml.resolver import BaseResolver, Resolver
from yaml.scanner import Scanner

from .reader import _ReadStream

_Loader: TypeAlias = Loader | BaseLoader | FullLoader | SafeLoader | UnsafeLoader  # noqa: Y047  # Used in other modules

class BaseLoader(Reader, Scanner, Parser, Composer, BaseConstructor, BaseResolver):
    def __init__(self, stream: _ReadStream) -> None: ...

class FullLoader(Reader, Scanner, Parser, Composer, FullConstructor, Resolver):
    def __init__(self, stream: _ReadStream) -> None: ...

class SafeLoader(Reader, Scanner, Parser, Composer, SafeConstructor, Resolver):
    def __init__(self, stream: _ReadStream) -> None: ...

class Loader(Reader, Scanner, Parser, Composer, Constructor, Resolver):
    def __init__(self, stream: _ReadStream) -> None: ...

class UnsafeLoader(Reader, Scanner, Parser, Composer, Constructor, Resolver):
    def __init__(self, stream: _ReadStream) -> None: ...
