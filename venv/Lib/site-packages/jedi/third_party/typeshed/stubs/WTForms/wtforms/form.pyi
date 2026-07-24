from _typeshed import SupportsItems
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any, ClassVar, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from wtforms.fields.core import Field, UnboundField
from wtforms.meta import DefaultMeta, _MultiDictLike

_T = TypeVar("_T")
_FormErrors: TypeAlias = dict[str, Sequence[str] | _FormErrors]

# _unbound_fields will always be a list on an instance, but on a
# class it might be None, if it never has been instantiated, or
# not instantianted after a new field had been added/removed
@type_check_only
class _UnboundFields(Protocol):
    @overload
    def __get__(self, obj: None, owner: type[object] | None = None, /) -> list[tuple[str, UnboundField[Any]]] | None: ...
    @overload
    def __get__(self, obj: object, owner: type[object] | None = None, /) -> list[tuple[str, UnboundField[Any]]]: ...

class BaseForm:
    meta: DefaultMeta
    form_errors: list[str]
    # we document this, because it's the only efficient way to introspect
    # the field names of the form, it also seems to be stable API-wise
    _fields: dict[str, Field]
    def __init__(
        self,
        fields: SupportsItems[str, UnboundField[Any]] | Iterable[tuple[str, UnboundField[Any]]],
        prefix: str = "",
        meta: DefaultMeta = ...,
    ) -> None: ...
    def __iter__(self) -> Iterator[Field]: ...
    def __contains__(self, name: str) -> bool: ...
    def __getitem__(self, name: str) -> Field: ...
    def __setitem__(self, name: str, value: UnboundField[Any]) -> None: ...
    def __delitem__(self, name: str) -> None: ...
    def populate_obj(self, obj: object) -> None: ...
    # while we would like to be more strict on extra_filters, we can't easily do that
    # without it being annoying in most situations
    def process(
        self,
        formdata: _MultiDictLike | None = None,
        obj: object | None = None,
        data: Mapping[str, Any] | None = None,
        extra_filters: Mapping[str, Sequence[Any]] | None = None,
        **kwargs: object,
    ) -> None: ...
    # same thing here with extra_validators
    def validate(self, extra_validators: Mapping[str, Sequence[Any]] | None = None) -> bool: ...
    @property
    def data(self) -> dict[str, Any]: ...
    # because of the Liskov violation in FormField.errors we need to make errors a recursive type
    @property
    def errors(self) -> _FormErrors: ...

class FormMeta(type):
    def __init__(cls, name: str, bases: Sequence[type[object]], attrs: Mapping[str, Any]) -> None: ...
    def __call__(cls: type[_T], *args: Any, **kwargs: Any) -> _T: ...
    def __setattr__(cls, name: str, value: object) -> None: ...
    def __delattr__(cls, name: str) -> None: ...

class Form(BaseForm, metaclass=FormMeta):
    # due to the metaclass this should always be a subclass of DefaultMeta
    # but if we annotate this as such, then subclasses cannot use it in the
    # intended way
    Meta: ClassVar[type[Any]]
    # this attribute is documented, so we annotate it
    _unbound_fields: _UnboundFields
    def __init__(
        self,
        formdata: _MultiDictLike | None = None,
        obj: object | None = None,
        prefix: str = "",
        data: Mapping[str, Any] | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        # same issue as with process
        extra_filters: Mapping[str, Sequence[Any]] | None = None,
        **kwargs: object,
    ) -> None: ...
    # this should emit a type_error, since it's not allowed to be called
    def __setitem__(self, name: str, value: None) -> None: ...  # type: ignore[override]
    def __delitem__(self, name: str) -> None: ...
    def __delattr__(self, name: str) -> None: ...

__all__ = ("BaseForm", "Form")
