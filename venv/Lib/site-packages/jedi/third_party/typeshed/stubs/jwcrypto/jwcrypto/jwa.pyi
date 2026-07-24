from abc import ABCMeta, abstractmethod
from collections.abc import Mapping
from typing import ClassVar

default_max_pbkdf2_iterations: int
default_enforce_hmac_key_length: bool

class JWAAlgorithm(metaclass=ABCMeta):
    @property
    @abstractmethod
    def name(self) -> str: ...
    @property
    @abstractmethod
    def description(self) -> str: ...
    @property
    @abstractmethod
    def keysize(self) -> int: ...
    @property
    @abstractmethod
    def algorithm_usage_location(self) -> str: ...
    @property
    @abstractmethod
    def algorithm_use(self) -> str: ...
    @property
    def input_keysize(self) -> int: ...

class JWA:
    algorithms_registry: ClassVar[Mapping[str, JWAAlgorithm]]
    @classmethod
    def instantiate_alg(cls, name: str, use: str | None = None) -> JWAAlgorithm: ...
    @classmethod
    def signing_alg(cls, name: str) -> JWAAlgorithm: ...
    @classmethod
    def keymgmt_alg(cls, name: str) -> JWAAlgorithm: ...
    @classmethod
    def encryption_alg(cls, name: str) -> JWAAlgorithm: ...
