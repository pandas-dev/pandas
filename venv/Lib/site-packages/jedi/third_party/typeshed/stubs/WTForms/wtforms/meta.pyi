from _typeshed import SupportsItems
from collections.abc import Collection, Iterator, MutableMapping
from typing import Any, Literal, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from markupsafe import Markup
from wtforms.fields.core import Field, UnboundField
from wtforms.form import BaseForm

_FieldT = TypeVar("_FieldT", bound=Field)

@type_check_only
class _SupportsGettextAndNgettext(Protocol):
    def gettext(self, string: str, /) -> str: ...
    def ngettext(self, singular: str, plural: str, n: int, /) -> str: ...

# these are the methods WTForms depends on, the dict can either provide
# a getlist or getall, if it only provides getall, it will wrapped, to
# provide getlist instead
@type_check_only
class _MultiDictLikeBase(Protocol):
    def __iter__(self) -> Iterator[str]: ...
    def __len__(self) -> int: ...
    def __contains__(self, key: Any, /) -> bool: ...

# since how file uploads are represented in formdata is implementation-specific
# we have to be generous in what we accept in the return of getlist/getall
# we can make this generic if we ever want to be more specific
@type_check_only
class _MultiDictLikeWithGetlist(_MultiDictLikeBase, Protocol):
    def getlist(self, key: str, /) -> list[Any]: ...

@type_check_only
class _MultiDictLikeWithGetall(_MultiDictLikeBase, Protocol):
    def getall(self, key: str, /) -> list[Any]: ...

_MultiDictLike: TypeAlias = _MultiDictLikeWithGetall | _MultiDictLikeWithGetlist

class DefaultMeta:
    def bind_field(self, form: BaseForm, unbound_field: UnboundField[_FieldT], options: MutableMapping[str, Any]) -> _FieldT: ...
    @overload
    def wrap_formdata(self, form: BaseForm, formdata: None) -> None: ...
    @overload
    def wrap_formdata(self, form: BaseForm, formdata: _MultiDictLike) -> _MultiDictLikeWithGetlist: ...
    def render_field(self, field: Field, render_kw: SupportsItems[str, Any]) -> Markup: ...
    csrf: bool
    csrf_field_name: str
    csrf_secret: Any | None
    csrf_context: Any | None
    csrf_class: type[Any] | None
    def build_csrf(self, form: BaseForm) -> Any: ...
    locales: Literal[False] | Collection[str]
    cache_translations: bool
    translations_cache: dict[str, _SupportsGettextAndNgettext]
    def get_translations(self, form: BaseForm) -> _SupportsGettextAndNgettext: ...
    def update_values(self, values: SupportsItems[str, Any]) -> None: ...
    # since meta can be extended with arbitrary data we add a __getattr__
    # method that returns Any
    def __getattr__(self, name: str) -> Any: ...
