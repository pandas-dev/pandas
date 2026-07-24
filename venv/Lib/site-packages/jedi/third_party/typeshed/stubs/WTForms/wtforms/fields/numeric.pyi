from collections.abc import Callable, Sequence
from decimal import Decimal
from typing import Any, Literal, overload
from typing_extensions import Self

from wtforms.fields.core import Field, _Filter, _FormT, _Validator, _Widget
from wtforms.form import BaseForm
from wtforms.meta import DefaultMeta, _SupportsGettextAndNgettext
from wtforms.utils import UnsetValue

__all__ = ("IntegerField", "DecimalField", "FloatField", "IntegerRangeField", "DecimalRangeField")

class LocaleAwareNumberField(Field):
    use_locale: bool
    number_format: Any | None
    locale: str
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        use_locale: bool = False,
        # this accepts a babel.numbers.NumberPattern, but since it
        # is an optional dependency we don't want to depend on it
        # for annotating this one argument
        number_format: str | Any | None = None,
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

class IntegerField(Field):
    data: int | None
    # technically this is not as strict and will accept anything
    # that can be passed into int(), but we might as well be
    default: int | Callable[[], int] | None
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        *,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: int | Callable[[], int] | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...

class DecimalField(LocaleAwareNumberField):
    data: Decimal | None
    # technically this is not as strict and will accept anything
    # that can be passed into Decimal(), but we might as well be
    default: Decimal | Callable[[], Decimal] | None
    places: int | None
    rounding: str | None
    @overload
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        *,
        places: UnsetValue = ...,
        rounding: None = None,
        use_locale: Literal[True],
        # this accepts a babel.numbers.NumberPattern, but since it
        # is an optional dependency we don't want to depend on it
        # for annotation this one argument
        number_format: str | Any | None = None,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: Decimal | Callable[[], Decimal] | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        places: int | UnsetValue | None = ...,
        rounding: str | None = None,
        *,
        use_locale: Literal[False] = False,
        # this accepts a babel.numbers.NumberPattern, but since it
        # is an optional dependency we don't want to depend on it
        # for annotation this one argument
        number_format: str | Any | None = None,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: Decimal | Callable[[], Decimal] | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...

class FloatField(Field):
    data: float | None
    # technically this is not as strict and will accept anything
    # that can be passed into float(), but we might as well be
    default: float | Callable[[], float] | None
    def __init__(
        self,
        label: str | None = None,
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
        *,
        filters: Sequence[_Filter] = (),
        description: str = "",
        id: str | None = None,
        default: float | Callable[[], float] | None = None,
        widget: _Widget[Self] | None = None,
        render_kw: dict[str, Any] | None = None,
        name: str | None = None,
        _form: BaseForm | None = None,
        _prefix: str = "",
        _translations: _SupportsGettextAndNgettext | None = None,
        _meta: DefaultMeta | None = None,
    ) -> None: ...

class IntegerRangeField(IntegerField): ...
class DecimalRangeField(DecimalField): ...
