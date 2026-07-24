from _typeshed import Incomplete
from collections.abc import Callable
from re import Pattern
from typing import Any, Final, Generic, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from ..rfc7517 import KeySet
from .claims import JWTClaims

_T = TypeVar("_T")

_LoadKey: TypeAlias = Callable[[Incomplete, Incomplete], Incomplete]

class JsonWebToken:
    SENSITIVE_NAMES: Final[tuple[str, ...]]
    SENSITIVE_VALUES: Final[Pattern[str]]

    def __init__(self, algorithms, private_headers=None) -> None: ...
    def check_sensitive_data(self, payload) -> None: ...
    def encode(self, header, payload, key, check: bool = True) -> bytes: ...
    @overload
    def decode(
        self,
        s: str | bytes,
        key: _LoadKey | KeySet | tuple[Incomplete, ...] | list[Incomplete] | str,
        claims_cls: None = None,
        claims_options=None,
        claims_params=None,
    ) -> JWTClaims: ...
    @overload
    def decode(
        self,
        s: str | bytes,
        key: _LoadKey | KeySet | tuple[Incomplete, ...] | list[Incomplete] | str,
        claims_cls: type[_T],
        claims_options=None,
        claims_params=None,
    ) -> _T: ...

def decode_payload(bytes_payload) -> dict[Incomplete, Incomplete]: ...

_TL = TypeVar("_TL", bound=tuple[Any, ...] | list[Any])

@type_check_only
class _Keys(TypedDict, Generic[_TL]):
    keys: _TL

@overload
def prepare_raw_key(raw: KeySet) -> KeySet: ...
@overload
def prepare_raw_key(raw: str) -> dict[str, Any] | str: ...  # dict is a JSON object
@overload
def prepare_raw_key(raw: _TL) -> _Keys[_TL]: ...
def find_encode_key(key, header): ...
def create_load_key(key: KeySet | _Keys[Incomplete] | Incomplete) -> _LoadKey: ...
