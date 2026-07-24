import datetime
from _typeshed import Unused
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Literal, NoReturn
from typing_extensions import Self, TypeAlias

from google.cloud.ndb import exceptions, key as key_module, query as query_module, tasklets as tasklets_module

Key = key_module.Key
Rollback = exceptions.Rollback
BlobKey: object
GeoPt: object

class KindError(exceptions.BadValueError): ...
class InvalidPropertyError(exceptions.Error): ...

BadProjectionError = InvalidPropertyError

class UnprojectedPropertyError(exceptions.Error): ...
class ReadonlyPropertyError(exceptions.Error): ...
class ComputedPropertyError(ReadonlyPropertyError): ...
class UserNotFoundError(exceptions.Error): ...

class _NotEqualMixin:
    def __ne__(self, other: object) -> bool: ...

_Direction: TypeAlias = Literal["asc", "desc"]

class IndexProperty(_NotEqualMixin):
    def __new__(cls, name: str, direction: _Direction) -> Self: ...
    @property
    def name(self) -> str: ...
    @property
    def direction(self) -> _Direction: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

class Index(_NotEqualMixin):
    def __new__(cls, kind: str, properties: list[IndexProperty], ancestor: bool) -> Self: ...
    @property
    def kind(self) -> str: ...
    @property
    def properties(self) -> list[IndexProperty]: ...
    @property
    def ancestor(self) -> bool: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...

class IndexState(_NotEqualMixin):
    def __new__(cls, definition, state, id): ...
    @property
    def definition(self): ...
    @property
    def state(self): ...
    @property
    def id(self): ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...

class ModelAdapter:
    # This actually returns NoReturn, but mypy can't handle that
    def __new__(cls, *args, **kwargs) -> Self: ...

def make_connection(*args, **kwargs) -> NoReturn: ...

class ModelAttribute: ...

class _BaseValue(_NotEqualMixin):
    b_val: object
    def __init__(self, b_val) -> None: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...

class Property(ModelAttribute):
    def __init__(
        self,
        name: str | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        required: bool | None = ...,
        default: object = None,
        choices: Iterable[object] | None = ...,
        validator: Callable[[Property, Any], object] | None = ...,
        verbose_name: str | None = ...,
        write_empty_list: bool | None = ...,
    ) -> None: ...
    def __eq__(self, value: object) -> query_module.FilterNode: ...  # type: ignore[override]
    def __ne__(self, value: object) -> query_module.FilterNode: ...  # type: ignore[override]
    def __lt__(self, value: object) -> query_module.FilterNode: ...
    def __le__(self, value: object) -> query_module.FilterNode: ...
    def __gt__(self, value: object) -> query_module.FilterNode: ...
    def __ge__(self, value: object) -> query_module.FilterNode: ...
    def IN(
        self, value: Iterable[object], server_op: bool = False
    ) -> query_module.DisjunctionNode | query_module.FilterNode | query_module.FalseNode: ...
    def NOT_IN(
        self, value: Iterable[object], server_op: bool = False
    ) -> query_module.DisjunctionNode | query_module.FilterNode | query_module.FalseNode: ...
    def __neg__(self) -> query_module.PropertyOrder: ...
    def __pos__(self) -> query_module.PropertyOrder: ...
    def __set__(self, entity: Model, value: object) -> None: ...
    def __delete__(self, entity: Model) -> None: ...

class ModelKey(Property):
    def __init__(self) -> None: ...
    def __get__(self, entity: Model, unused_cls: type[Model] | None = ...) -> key_module.Key | list[key_module.Key] | None: ...

class BooleanProperty(Property):
    def __get__(self, entity: Model, unused_cls: type[Model] | None = ...) -> bool | list[bool] | None: ...

class IntegerProperty(Property):
    def __get__(self, entity: Model, unused_cls: type[Model] | None = ...) -> int | list[int] | None: ...

class FloatProperty(Property):
    def __get__(self, entity: Model, unused_cls: type[Model] | None = ...) -> float | list[float] | None: ...

class _CompressedValue(bytes):
    z_val: bytes
    def __init__(self, z_val: bytes) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> NoReturn: ...

class BlobProperty(Property):
    def __init__(
        self,
        name: str | None = ...,
        compressed: bool | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        required: bool | None = ...,
        default: bytes | None = ...,
        choices: Iterable[bytes] | None = ...,
        validator: Callable[[Property, Any], object] | None = ...,
        verbose_name: str | None = ...,
        write_empty_list: bool | None = ...,
    ) -> None: ...
    def __get__(self, entity: Model, unused_cls: type[Model] | None = ...) -> bytes | list[bytes] | None: ...

class CompressedTextProperty(BlobProperty):
    __slots__ = ()
    def __init__(self, *args, **kwargs) -> None: ...

class TextProperty(Property):
    def __new__(cls, *args, **kwargs): ...
    def __init__(self, *args, **kwargs) -> None: ...
    def __get__(self, entity: Model, unused_cls: type[Model] | None = ...) -> str | list[str] | None: ...

class StringProperty(TextProperty):
    def __init__(self, *args, **kwargs) -> None: ...

class GeoPtProperty(Property): ...
class PickleProperty(BlobProperty): ...

class JsonProperty(BlobProperty):
    def __init__(
        self,
        name: str | None = ...,
        compressed: bool | None = ...,
        json_type: type | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        required: bool | None = ...,
        default: object = None,
        choices: Iterable[object] | None = ...,
        validator: Callable[[Property, Any], object] | None = ...,
        verbose_name: str | None = ...,
        write_empty_list: bool | None = ...,
    ) -> None: ...

class User:
    def __init__(self, email: str | None = ..., _auth_domain: str | None = ..., _user_id: str | None = ...) -> None: ...
    def nickname(self) -> str: ...
    def email(self): ...
    def user_id(self) -> str | None: ...
    def auth_domain(self) -> str: ...
    def __hash__(self) -> int: ...
    def __eq__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...

class UserProperty(Property):
    def __init__(
        self,
        name: str | None = ...,
        auto_current_user: bool | None = ...,
        auto_current_user_add: bool | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        required: bool | None = ...,
        default: bytes | None = ...,
        choices: Iterable[bytes] | None = ...,
        validator: Callable[[Property, Any], object] | None = ...,
        verbose_name: str | None = ...,
        write_empty_list: bool | None = ...,
    ) -> None: ...

class KeyProperty(Property):
    def __init__(
        self,
        name: str | None = ...,
        kind: type[Model] | str | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        required: bool | None = ...,
        default: key_module.Key | None = ...,
        choices: Iterable[key_module.Key] | None = ...,
        validator: Callable[[Property, key_module.Key], key_module.Key] | None = ...,
        verbose_name: str | None = ...,
        write_empty_list: bool | None = ...,
    ) -> None: ...

class BlobKeyProperty(Property): ...

class DateTimeProperty(Property):
    def __init__(
        self,
        name: str | None = ...,
        auto_now: bool | None = ...,
        auto_now_add: bool | None = ...,
        tzinfo: datetime.tzinfo | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        required: bool | None = ...,
        default: datetime.datetime | None = ...,
        choices: Iterable[datetime.datetime] | None = ...,
        validator: Callable[[Property, Any], object] | None = ...,
        verbose_name: str | None = ...,
        write_empty_list: bool | None = ...,
    ) -> None: ...

class DateProperty(DateTimeProperty): ...
class TimeProperty(DateTimeProperty): ...

class StructuredProperty(Property):
    def __init__(self, model_class: type, name: str | None = ..., **kwargs) -> None: ...
    def __getattr__(self, attrname: str): ...
    def IN(self, value: Iterable[object]) -> query_module.DisjunctionNode | query_module.FalseNode: ...  # type: ignore[override]

class LocalStructuredProperty(BlobProperty):
    def __init__(self, model_class: type[Model], **kwargs) -> None: ...

class GenericProperty(Property):
    def __init__(self, name: str | None = ..., compressed: bool = ..., **kwargs) -> None: ...

class ComputedProperty(GenericProperty):
    def __init__(
        self,
        func: Callable[[Model], object],
        name: str | None = ...,
        indexed: bool | None = ...,
        repeated: bool | None = ...,
        verbose_name: str | None = ...,
    ) -> None: ...

class MetaModel(type):
    def __init__(cls, name: str, bases, classdict) -> None: ...

class Model(_NotEqualMixin, metaclass=MetaModel):
    key: ModelKey
    def __init__(_self, **kwargs) -> None: ...
    def __hash__(self) -> NoReturn: ...
    def __eq__(self, other: object) -> bool: ...
    @classmethod
    def gql(cls: type[Model], query_string: str, *args, **kwargs) -> query_module.Query: ...
    def put(self, **kwargs): ...
    def put_async(self, **kwargs) -> tasklets_module.Future: ...
    @classmethod
    def query(cls: type[Model], *args, **kwargs) -> query_module.Query: ...
    @classmethod
    def allocate_ids(
        cls: type[Model],
        size: int | None = ...,
        max: int | None = ...,
        parent: key_module.Key | None = ...,
        retries: int | None = ...,
        timeout: float | None = ...,
        deadline: float | None = ...,
        use_cache: bool | None = ...,
        use_global_cache: bool | None = ...,
        global_cache_timeout: int | None = ...,
        use_datastore: bool | None = ...,
        use_memcache: bool | None = ...,
        memcache_timeout: int | None = ...,
        max_memcache_items: int | None = ...,
        force_writes: bool | None = ...,
        _options=...,
    ) -> tuple[key_module.Key, key_module.Key]: ...
    @classmethod
    def allocate_ids_async(
        cls: type[Model],
        size: int | None = ...,
        max: int | None = ...,
        parent: key_module.Key | None = ...,
        retries: int | None = ...,
        timeout: float | None = ...,
        deadline: float | None = ...,
        use_cache: bool | None = ...,
        use_global_cache: bool | None = ...,
        global_cache_timeout: int | None = ...,
        use_datastore: bool | None = ...,
        use_memcache: bool | None = ...,
        memcache_timeout: int | None = ...,
        max_memcache_items: int | None = ...,
        force_writes: bool | None = ...,
        _options=...,
    ) -> tasklets_module.Future: ...
    @classmethod
    def get_by_id(
        cls: type[Model],
        id: int | str | None,
        parent: key_module.Key | None = ...,
        namespace: str | None = ...,
        project: str | None = ...,
        app: str | None = ...,
        read_consistency: Literal["EVENTUAL"] | None = ...,
        read_policy: Literal["EVENTUAL"] | None = ...,
        transaction: bytes | None = ...,
        retries: int | None = ...,
        timeout: float | None = ...,
        deadline: float | None = ...,
        use_cache: bool | None = ...,
        use_global_cache: bool | None = ...,
        global_cache_timeout: int | None = ...,
        use_datastore: bool | None = ...,
        use_memcache: bool | None = ...,
        memcache_timeout: int | None = ...,
        max_memcache_items: int | None = ...,
        force_writes: bool | None = ...,
        _options=...,
        database: str | None = None,
    ) -> Model | None: ...
    @classmethod
    def get_by_id_async(
        cls: type[Model],
        id: int | str,
        parent: key_module.Key | None = ...,
        namespace: str | None = ...,
        project: str | None = ...,
        app: str | None = ...,
        read_consistency: Literal["EVENTUAL"] | None = ...,
        read_policy: Literal["EVENTUAL"] | None = ...,
        transaction: bytes | None = ...,
        retries: int | None = ...,
        timeout: float | None = ...,
        deadline: float | None = ...,
        use_cache: bool | None = ...,
        use_global_cache: bool | None = ...,
        global_cache_timeout: int | None = ...,
        use_datastore: bool | None = ...,
        use_memcache: bool | None = ...,
        memcache_timeout: int | None = ...,
        max_memcache_items: int | None = ...,
        force_writes: bool | None = ...,
        _options=...,
        database: str | None = None,
    ) -> tasklets_module.Future: ...
    @classmethod
    def get_or_insert(
        cls: type[Model],
        _name: str,
        parent: key_module.Key | None = ...,
        namespace: str | None = ...,
        project: str | None = ...,
        app: str | None = ...,
        read_consistency: Literal["EVENTUAL"] | None = ...,
        read_policy: Literal["EVENTUAL"] | None = ...,
        transaction: bytes | None = ...,
        retries: int | None = ...,
        timeout: float | None = ...,
        deadline: float | None = ...,
        use_cache: bool | None = ...,
        use_global_cache: bool | None = ...,
        global_cache_timeout: int | None = ...,
        use_datastore: bool | None = ...,
        use_memcache: bool | None = ...,
        memcache_timeout: int | None = ...,
        max_memcache_items: int | None = ...,
        force_writes: bool | None = ...,
        _options=...,
        **kw_model_args,
    ) -> Model: ...
    @classmethod
    def get_or_insert_async(
        cls: type[Model],
        _name: str,
        parent: key_module.Key | None = ...,
        namespace: str | None = ...,
        project: str | None = ...,
        app: str | None = ...,
        read_consistency: Literal["EVENTUAL"] | None = ...,
        read_policy: Literal["EVENTUAL"] | None = ...,
        transaction: bytes | None = ...,
        retries: int | None = ...,
        timeout: float | None = ...,
        deadline: float | None = ...,
        use_cache: bool | None = ...,
        use_global_cache: bool | None = ...,
        global_cache_timeout: int | None = ...,
        use_datastore: bool | None = ...,
        use_memcache: bool | None = ...,
        memcache_timeout: int | None = ...,
        max_memcache_items: int | None = ...,
        force_writes: bool | None = ...,
        _options=...,
        **kw_model_args,
    ) -> tasklets_module.Future: ...
    def populate(self, **kwargs) -> None: ...
    def has_complete_key(self) -> bool: ...
    def to_dict(
        self,
        include: list[object] | tuple[object, object] | set[object] | None = ...,
        exclude: list[object] | tuple[object, object] | set[object] | None = ...,
    ): ...

class Expando(Model):
    def __getattr__(self, name: str): ...
    def __setattr__(self, name: str, value) -> None: ...
    def __delattr__(self, name: str) -> None: ...

def get_multi_async(
    keys: Sequence[key_module.Key],
    read_consistency: Literal["EVENTUAL"] | None = ...,
    read_policy: Literal["EVENTUAL"] | None = ...,
    transaction: bytes | None = ...,
    retries: int | None = ...,
    timeout: float | None = ...,
    deadline: float | None = ...,
    use_cache: bool | None = ...,
    use_global_cache: bool | None = ...,
    global_cache_timeout: int | None = ...,
    use_datastore: bool | None = ...,
    use_memcache: bool | None = ...,
    memcache_timeout: int | None = ...,
    max_memcache_items: int | None = ...,
    force_writes: bool | None = ...,
    _options: object = None,
) -> list[tasklets_module.Future]: ...
def get_multi(
    keys: Sequence[key_module.Key],
    read_consistency: Literal["EVENTUAL"] | None = ...,
    read_policy: Literal["EVENTUAL"] | None = ...,
    transaction: bytes | None = ...,
    retries: int | None = ...,
    timeout: float | None = ...,
    deadline: float | None = ...,
    use_cache: bool | None = ...,
    use_global_cache: bool | None = ...,
    global_cache_timeout: int | None = ...,
    use_datastore: bool | None = ...,
    use_memcache: bool | None = ...,
    memcache_timeout: int | None = ...,
    max_memcache_items: int | None = ...,
    force_writes: bool | None = ...,
    _options: object = None,
) -> list[Model | None]: ...
def put_multi_async(
    entities: list[Model],
    retries: int | None = ...,
    timeout: float | None = ...,
    deadline: float | None = ...,
    use_cache: bool | None = ...,
    use_global_cache: bool | None = ...,
    global_cache_timeout: int | None = ...,
    use_datastore: bool | None = ...,
    use_memcache: bool | None = ...,
    memcache_timeout: int | None = ...,
    max_memcache_items: int | None = ...,
    force_writes: bool | None = ...,
    _options: object = None,
) -> list[tasklets_module.Future]: ...
def put_multi(
    entities: list[Model],
    retries: int | None = ...,
    timeout: float | None = ...,
    deadline: float | None = ...,
    use_cache: bool | None = ...,
    use_global_cache: bool | None = ...,
    global_cache_timeout: int | None = ...,
    use_datastore: bool | None = ...,
    use_memcache: bool | None = ...,
    memcache_timeout: int | None = ...,
    max_memcache_items: int | None = ...,
    force_writes: bool | None = ...,
    _options: object = None,
) -> list[key_module.Key]: ...
def delete_multi_async(
    keys: Sequence[key_module.Key],
    retries: int | None = ...,
    timeout: float | None = ...,
    deadline: float | None = ...,
    use_cache: bool | None = ...,
    use_global_cache: bool | None = ...,
    global_cache_timeout: int | None = ...,
    use_datastore: bool | None = ...,
    use_memcache: bool | None = ...,
    memcache_timeout: int | None = ...,
    max_memcache_items: int | None = ...,
    force_writes: bool | None = ...,
    _options: object = None,
) -> list[tasklets_module.Future]: ...
def delete_multi(
    keys: Sequence[key_module.Key],
    retries: int | None = ...,
    timeout: float | None = ...,
    deadline: float | None = ...,
    use_cache: bool | None = ...,
    use_global_cache: bool | None = ...,
    global_cache_timeout: int | None = ...,
    use_datastore: bool | None = ...,
    use_memcache: bool | None = ...,
    memcache_timeout: int | None = ...,
    max_memcache_items: int | None = ...,
    force_writes: bool | None = ...,
    _options: object = None,
) -> list[None]: ...
def get_indexes_async(**options: Unused) -> NoReturn: ...
def get_indexes(**options: Unused) -> NoReturn: ...
