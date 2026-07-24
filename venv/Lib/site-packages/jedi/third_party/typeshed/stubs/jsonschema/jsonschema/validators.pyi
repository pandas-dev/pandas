from _typeshed import Incomplete, SupportsKeysAndGetItem
from collections.abc import Callable, Generator, Iterable, Iterator, Mapping
from contextlib import contextmanager
from typing import Any, ClassVar, overload, type_check_only
from typing_extensions import TypeAlias, deprecated

from referencing.jsonschema import Schema, SchemaRegistry
from referencing.typing import URI

from ._format import FormatChecker
from ._types import TypeChecker
from ._utils import Unset, URIDict
from .exceptions import ValidationError
from .protocols import Validator

# these type aliases do not exist at runtime, they're only defined here in the stub
_JsonObject: TypeAlias = Mapping[str, Any]
_JsonValue: TypeAlias = _JsonObject | list[Any] | str | int | float | bool | None
_ValidatorCallback: TypeAlias = Callable[[Any, Any, _JsonValue, _JsonObject], Iterator[ValidationError]]

# This class does not exist at runtime. Compatible classes are created at
# runtime by create().
@type_check_only
class _Validator(Validator):
    VALIDATORS: ClassVar[dict[Incomplete, Incomplete]]
    META_SCHEMA: ClassVar[dict[Incomplete, Incomplete]]
    TYPE_CHECKER: ClassVar[Incomplete]
    FORMAT_CHECKER: ClassVar[Incomplete]
    @staticmethod
    def ID_OF(contents: Schema) -> URI | None: ...
    schema: Schema
    format_checker: FormatChecker | None
    def __init__(
        self,
        schema: Mapping[Incomplete, Incomplete] | bool,
        resolver: Any = None,  # deprecated
        format_checker: FormatChecker | None = None,
        *,
        registry: SchemaRegistry = ...,
    ) -> None: ...
    @classmethod
    def check_schema(cls, schema: Schema, format_checker: FormatChecker | Unset = ...) -> None: ...
    @property
    @deprecated(
        "Accessing resolver() is deprecated as of v4.18.0, "
        "in favor of the https://github.com/python-jsonschema/referencing library."
    )
    def resolver(self): ...
    def evolve(self, **changes) -> _Validator: ...
    @overload
    def iter_errors(self, instance) -> Generator[Incomplete]: ...
    @overload
    @deprecated("Passing a schema to Validator.iter_errors is deprecated and will be removed in a future release.")
    def iter_errors(self, instance, _schema: Schema | None) -> Generator[Incomplete]: ...
    def descend(
        self, instance, schema: Schema, path: Incomplete | None = ..., schema_path: Incomplete | None = ..., resolver=None
    ) -> Generator[Incomplete]: ...
    def validate(self, *args, **kwargs) -> None: ...
    def is_type(self, instance, type) -> bool: ...
    @overload
    def is_valid(self, instance) -> bool: ...
    @overload
    @deprecated("Passing a schema to Validator.is_valid is deprecated and will be removed in a future release.")
    def is_valid(self, instance, _schema: Schema | None) -> bool: ...

def validates(version: str) -> Callable[[_Validator], _Validator]: ...
def create(
    meta_schema: Schema,
    validators: Mapping[str, _ValidatorCallback] | tuple[()] = (),
    version=None,
    type_checker: TypeChecker = ...,
    format_checker: FormatChecker = ...,
    id_of: Callable[[Schema], str] = ...,
    applicable_validators: Callable[[Schema], Iterable[tuple[str, _ValidatorCallback]]] = ...,
) -> type[_Validator]: ...
def extend(validator, validators=(), version=None, type_checker=None, format_checker=None): ...

# At runtime these are fields that are assigned the return values of create() calls.
class Draft3Validator(_Validator):
    __slots__ = ("_validators", "schema", "_ref_resolver", "format_checker", "_registry", "_resolver", "__weakref__")

class Draft4Validator(_Validator):
    __slots__ = ("_validators", "schema", "_ref_resolver", "format_checker", "_registry", "_resolver", "__weakref__")

class Draft6Validator(_Validator):
    __slots__ = ("_validators", "schema", "_ref_resolver", "format_checker", "_registry", "_resolver", "__weakref__")

class Draft7Validator(_Validator):
    __slots__ = ("_validators", "schema", "_ref_resolver", "format_checker", "_registry", "_resolver", "__weakref__")

class Draft201909Validator(_Validator):
    __slots__ = ("_validators", "schema", "_ref_resolver", "format_checker", "_registry", "_resolver", "__weakref__")

class Draft202012Validator(_Validator):
    __slots__ = ("_validators", "schema", "_ref_resolver", "format_checker", "_registry", "_resolver", "__weakref__")

_Handler: TypeAlias = Callable[[str], Incomplete]

@deprecated(
    "jsonschema.RefResolver is deprecated as of v4.18.0, in favor of the "
    "https://github.com/python-jsonschema/referencing library, which "
    "provides more compliant referencing behavior as well as more "
    "flexible APIs for customization. A future release will remove "
    "RefResolver. Please file a feature request (on referencing) if you "
    "are missing an API for the kind of customization you need."
)
class RefResolver:
    referrer: dict[str, Incomplete]
    cache_remote: Incomplete
    handlers: dict[str, _Handler]
    store: URIDict
    def __init__(
        self,
        base_uri: str,
        referrer: dict[str, Incomplete],
        store: Mapping[str, Mapping[str, Any]] | Iterable[tuple[str, Mapping[str, Any]]] = ...,
        cache_remote: bool = True,
        handlers: SupportsKeysAndGetItem[str, _Handler] | Iterable[tuple[str, _Handler]] = (),
        urljoin_cache=None,
        remote_cache=None,
    ) -> None: ...
    @classmethod
    def from_schema(cls, schema: Schema, id_of=..., *args, **kwargs): ...
    def push_scope(self, scope) -> None: ...
    def pop_scope(self) -> None: ...
    @property
    def resolution_scope(self): ...
    @property
    def base_uri(self) -> str: ...
    @contextmanager
    @deprecated("jsonschema.RefResolver.in_scope is deprecated and will be removed in a future release.")
    def in_scope(self, scope) -> Generator[None]: ...
    @contextmanager
    def resolving(self, ref: str) -> Generator[Incomplete]: ...
    def resolve(self, ref: str) -> tuple[str, Incomplete]: ...
    def resolve_from_url(self, url: str): ...
    def resolve_fragment(self, document, fragment: str): ...
    def resolve_remote(self, uri: str): ...

def validate(instance: object, schema: Schema, cls: type[_Validator] | None = None, *args: Any, **kwargs: Any) -> None: ...
def validator_for(schema: Schema | bool, default: type[Validator] | Unset = ...) -> type[Validator]: ...
