from collections.abc import Iterator
from typing import Any, Generic, TypeVar

from wtforms.fields.core import Field, _Widget
from wtforms.form import BaseForm, _FormErrors
from wtforms.meta import DefaultMeta, _SupportsGettextAndNgettext

__all__ = ("FormField",)

_BoundFormT = TypeVar("_BoundFormT", bound=BaseForm)

class FormField(Field, Generic[_BoundFormT]):
    form_class: type[_BoundFormT]
    form: _BoundFormT
    separator: str
    def __init__(
        self: FormField[_BoundFormT],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        form_class: type[_BoundFormT],
        label: str | None = None,
        validators: None = None,
        separator: str = "-",
        *,
        description: str = "",
        id: str | None = None,
        default: object | None = None,
        widget: _Widget[FormField[_BoundFormT]] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...
    def __iter__(self) -> Iterator[Field]: ...
    def __getitem__(self, name: str) -> Field: ...
    def __getattr__(self, name: str) -> Field: ...
    @property
    def data(self) -> dict[str, Any]: ...
    @property
    def errors(self) -> _FormErrors: ...  # type: ignore[override]
