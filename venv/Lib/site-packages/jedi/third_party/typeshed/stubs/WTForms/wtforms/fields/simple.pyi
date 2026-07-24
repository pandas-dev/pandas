from collections.abc import Callable, Collection, Sequence
from typing import Any
from typing_extensions import Self

from wtforms.fields.core import Field, _Filter, _FormT, _Validator, _Widget
from wtforms.form import BaseForm
from wtforms.meta import DefaultMeta, _SupportsGettextAndNgettext

__all__ = (
    "BooleanField",
    "TextAreaField",
    "PasswordField",
    "FileField",
    "MultipleFileField",
    "HiddenField",
    "SearchField",
    "SubmitField",
    "StringField",
    "TelField",
    "URLField",
    "EmailField",
    "ColorField",
)

class BooleanField(Field):
    data: bool
    default: bool | Callable[[], bool] | None
    false_values: Collection[Any]
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        false_values: Collection[Any] | None = None,
        *,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: bool | Callable[[], bool] | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...

class StringField(Field):
    data: str | None
    default: str | Callable[[], str] | None
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: str | Callable[[], str] | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...

class TextAreaField(StringField): ...
class PasswordField(StringField): ...
class FileField(Field): ...

class MultipleFileField(FileField):
    data: list[Any]

class HiddenField(StringField): ...
class SubmitField(BooleanField): ...
class SearchField(StringField): ...
class TelField(StringField): ...
class URLField(StringField): ...
class EmailField(StringField): ...
class ColorField(StringField): ...
