from _typeshed import Unused
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any, NamedTuple
from typing_extensions import TypeAlias

from django import forms
from django.db.models import Choices
from django.forms import Widget
from django_stubs_ext import StrOrPromise

DJANGO_50: bool

# Ref: django-stubs/forms/fields.pyi
# Problem: attribute `widget` is always of type `Widget` after field instantiation.
# However, on class level it can be set to `Type[Widget]` too.
# If we annotate it as `Union[Widget, Type[Widget]]`, every code that uses field
# instances will not typecheck.
# If we annotate it as `Widget`, any widget subclasses that do e.g.
# `widget = Select` will not typecheck.
# `Any` gives too much freedom, but does not create false positives.
_ClassLevelWidget: TypeAlias = Any
# Validator parameter type depends on type of the form field used.
_ValidatorCallable: TypeAlias = Callable[[Any], None]
# Based on django-stubs utils/choices.pyi
_Choice: TypeAlias = tuple[Any, Any]
_ChoiceNamedGroup: TypeAlias = tuple[str, Iterable[_Choice]]
_Choices: TypeAlias = Iterable[_Choice | _ChoiceNamedGroup]
_ChoicesMapping: TypeAlias = Mapping[Any, Any]
_ChoicesInput: TypeAlias = _Choices | _ChoicesMapping | type[Choices] | Callable[[], _Choices | _ChoicesMapping]

class RangeField(forms.MultiValueField):
    widget: _ClassLevelWidget = ...
    def __init__(
        self,
        fields: tuple[forms.Field, forms.Field] | None = None,
        *,
        # Inherited from Django MultiValueField
        require_all_fields: bool = True,
        required: bool = ...,
        widget: Widget | type[Widget] | None = ...,
        label: StrOrPromise | None = ...,
        initial: Any | None = ...,  # Type depends on the form field used.
        help_text: StrOrPromise = ...,
        error_messages: Mapping[str, StrOrPromise] | None = ...,
        show_hidden_initial: bool = ...,
        validators: Sequence[_ValidatorCallable] = ...,
        localize: bool = ...,
        disabled: bool = ...,
        label_suffix: str | None = ...,
    ) -> None: ...  # Args/kwargs can be any field params, passes to parent
    def compress(self, data_list: list[Any] | None) -> slice | None: ...  # Data list elements can be any field value type

class DateRangeField(RangeField):
    widget: _ClassLevelWidget = ...
    def compress(self, data_list: list[Any] | None) -> slice | None: ...  # Date values in list can be any date type

class DateTimeRangeField(RangeField):
    widget: _ClassLevelWidget = ...

class IsoDateTimeRangeField(RangeField):
    widget: _ClassLevelWidget = ...

class TimeRangeField(RangeField):
    widget: _ClassLevelWidget = ...

class Lookup(NamedTuple):
    value: Any  # Lookup values can be any filterable type
    lookup_expr: str

class LookupChoiceField(forms.MultiValueField):
    def __init__(
        self,
        field: forms.Field,
        lookup_choices: Sequence[tuple[str, str]],
        *,
        empty_label: StrOrPromise = ...,
        widget: Unused = ...,
        help_text: Unused = ...,
        # Inherited from Django MultiValueField
        require_all_fields: bool = True,
        required: bool = ...,
        label: StrOrPromise | None = ...,
        initial: Any | None = ...,  # Type depends on the form field used.
        error_messages: Mapping[str, StrOrPromise] | None = ...,
        show_hidden_initial: bool = ...,
        validators: Sequence[_ValidatorCallable] = ...,
        localize: bool = ...,
        disabled: bool = ...,
        label_suffix: str | None = ...,
    ) -> None: ...  # Args/kwargs can be any field params, uses kwargs for empty_label
    def compress(self, data_list: list[Any] | None) -> Lookup | None: ...  # Data list can contain any lookup components

class IsoDateTimeField(forms.DateTimeField):
    ISO_8601: str
    input_formats: list[str]
    def strptime(self, value: str, format: str) -> Any: ...  # Returns datetime objects or parsing results

class BaseCSVField(forms.Field):
    base_widget_class: _ClassLevelWidget = ...
    def clean(self, value: Any) -> Any: ...  # Cleaned values can be any valid field type

class BaseRangeField(BaseCSVField):
    widget: _ClassLevelWidget = ...
    def clean(self, value: Any) -> Any: ...  # Input and output values can be any range type

class ChoiceIterator:
    field: ChoiceField
    choices: Sequence[tuple[Any, str]]  # Choice values can be any type (int, str, Model, etc.)
    def __init__(
        self, field: ChoiceField, choices: Sequence[tuple[Any, str]]
    ) -> None: ...  # Choice values can be any selectable type
    def __iter__(self) -> Any: ...  # Iterator yields choice tuples with any value types
    def __len__(self) -> int: ...

class ModelChoiceIterator(forms.models.ModelChoiceIterator):
    def __iter__(self) -> Any: ...  # Iterator yields choice tuples with any value types
    def __len__(self) -> int: ...

class ChoiceIteratorMixin:
    null_label: StrOrPromise | None
    null_value: Any  # Null choice values can be any type (None, empty string, etc.)
    def __init__(self, *, null_label: StrOrPromise | None, null_value: Any) -> None: ...

class ChoiceField(ChoiceIteratorMixin, forms.ChoiceField):
    iterator = ChoiceIterator
    empty_label: StrOrPromise
    def __init__(
        self,
        *,
        empty_label: StrOrPromise = ...,
        # Inherited from Django ChoiceField
        choices: _ChoicesInput = (),
        required: bool = ...,
        widget: Widget | type[Widget] | None = ...,
        label: StrOrPromise | None = ...,
        initial: Any | None = ...,  # Type depends on the form field used.
        help_text: StrOrPromise = ...,
        error_messages: Mapping[str, StrOrPromise] | None = ...,
        show_hidden_initial: bool = ...,
        validators: Sequence[_ValidatorCallable] = ...,
        localize: bool = ...,
        disabled: bool = ...,
        label_suffix: str | None = ...,
        null_label: StrOrPromise | None,
        null_value: Any,  # Type depends on the form field used.
    ) -> None: ...

class MultipleChoiceField(ChoiceIteratorMixin, forms.MultipleChoiceField):
    iterator = ChoiceIterator
    empty_label: StrOrPromise | None

class ModelChoiceField(ChoiceIteratorMixin, forms.ModelChoiceField[Any]):
    iterator = ModelChoiceIterator
    def to_python(self, value: Any) -> Any: ...  # Converts any input to Python model objects or values

class ModelMultipleChoiceField(ChoiceIteratorMixin, forms.ModelMultipleChoiceField[Any]):
    iterator = ModelChoiceIterator
