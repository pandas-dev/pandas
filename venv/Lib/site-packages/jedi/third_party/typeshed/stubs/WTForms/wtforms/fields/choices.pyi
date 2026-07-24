from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Any
from typing_extensions import Self, TypeAlias

from wtforms.fields.core import Field, _Filter, _FormT, _Validator, _Widget
from wtforms.form import BaseForm
from wtforms.meta import DefaultMeta, _SupportsGettextAndNgettext

__all__ = ("SelectField", "SelectMultipleField", "RadioField")

# technically this allows a list, but we're more strict for type safety
_Choice: TypeAlias = tuple[Any, str] | tuple[Any, str, dict[str, Any]]
# it's too difficult to get type safety here due to to nested partially invariant collections
_GroupedChoices: TypeAlias = dict[str, Any]  # Any should be Collection[_Choice]
_FullChoice: TypeAlias = tuple[Any, str, bool, dict[str, Any]]  # value, label, selected, render_kw
_FullGroupedChoices: TypeAlias = tuple[str, Iterable[_FullChoice]]
_Option: TypeAlias = SelectFieldBase._Option

class SelectFieldBase(Field):
    option_widget: _Widget[_Option]
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        option_widget: _Widget[_Option] | None = None,
        *,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: object | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...
    def iter_choices(self) -> Iterator[_FullChoice]: ...
    def has_groups(self) -> bool: ...
    def iter_groups(self) -> Iterator[_FullGroupedChoices]: ...
    def __iter__(self) -> Iterator[_Option]: ...

    class _Option(Field):
        checked: bool

class SelectField(SelectFieldBase):
    coerce: Callable[[Any], Any]
    choices: Sequence[_Choice] | _GroupedChoices | None
    validate_choice: bool
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        coerce: Callable[[Any], Any] = ...,
        choices: Iterable[_Choice] | _GroupedChoices | Callable[[], Iterable[_Choice] | _GroupedChoices] | None = None,
        validate_choice: bool = True,
        *,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: object | None = None,
        widget: _Widget[Self] | None = None,
        option_widget: _Widget[_Option] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...
    def iter_choices(self) -> Iterator[_FullChoice]: ...
    def has_groups(self) -> bool: ...
    def iter_groups(self) -> Iterator[_FullGroupedChoices]: ...

class SelectMultipleField(SelectField):
    data: list[Any] | None

class RadioField(SelectField): ...
