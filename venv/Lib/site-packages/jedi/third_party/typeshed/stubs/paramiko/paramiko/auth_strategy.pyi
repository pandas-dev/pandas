import abc
from collections.abc import Callable, Iterator
from logging import Logger
from pathlib import Path
from typing import NamedTuple

from paramiko.config import SSHConfig
from paramiko.pkey import PKey
from paramiko.ssh_exception import AuthenticationException
from paramiko.transport import Transport

class AuthSource:
    username: str
    def __init__(self, username: str) -> None: ...
    @abc.abstractmethod
    def authenticate(self, transport: Transport) -> list[str]: ...

class NoneAuth(AuthSource):
    def authenticate(self, transport: Transport) -> list[str]: ...

class Password(AuthSource):
    password_getter: Callable[[], str]
    def __init__(self, username: str, password_getter: Callable[[], str]) -> None: ...
    def authenticate(self, transport: Transport) -> list[str]: ...

class PrivateKey(AuthSource):
    def authenticate(self, transport: Transport) -> list[str]: ...

class InMemoryPrivateKey(PrivateKey):
    pkey: PKey
    def __init__(self, username: str, pkey: PKey) -> None: ...

class OnDiskPrivateKey(PrivateKey):
    source: str
    path: Path
    pkey: PKey
    def __init__(self, username: str, source: str, path: Path, pkey: PKey) -> None: ...

class SourceResult(NamedTuple):
    source: AuthSource
    result: list[str] | Exception

class AuthResult(list[SourceResult]):
    strategy: AuthStrategy
    def __init__(self, strategy: AuthStrategy, *args: SourceResult, **kwargs: object) -> None: ...

class AuthFailure(AuthenticationException):
    result: AuthResult
    def __init__(self, result: AuthResult) -> None: ...

class AuthStrategy:
    ssh_config: SSHConfig
    log: Logger
    def __init__(self, ssh_config: SSHConfig) -> None: ...
    @abc.abstractmethod
    def get_sources(self) -> Iterator[AuthSource]: ...
    def authenticate(self, transport: Transport) -> list[SourceResult]: ...
