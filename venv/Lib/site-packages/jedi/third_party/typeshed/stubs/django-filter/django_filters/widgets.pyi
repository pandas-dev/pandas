from collections.abc import Mapping, Sequence
from typing import Any

from django import forms
from django.http import QueryDict
from django.utils.safestring import SafeString

class LinkWidget(forms.Widget):
    # Choice values can be any type (int, str, Model, etc.)
    choices: Sequence[tuple[Any, str]]
    # Choice values can be any selectable type
    def __init__(self, attrs: dict[str, Any] | None = None, choices: Sequence[tuple[Any, str]] = ()) -> None: ...
    data: QueryDict | dict[str, Any]
    # Return value depends on widget data type
    def value_from_datadict(self, data: Mapping[str, Any], files: Mapping[str, Any], name: str) -> Any: ...
    # Widget value and renderer can be any type, choices parameter combines with class choices
    def render(  # type: ignore[override]
        self,
        name: str,
        value: Any,
        attrs: dict[str, Any] | None = None,
        choices: Sequence[tuple[Any, str]] = (),
        renderer: Any | None = None,
    ) -> SafeString: ...
    # Choice values and selections can be any type
    def render_options(self, choices: Sequence[tuple[Any, str]], selected_choices: list[Any], name: str) -> str: ...
    # Selected choices and option values can be any type
    def render_option(self, name: str, selected_choices: list[Any], option_value: Any, option_label: str) -> str: ...
    def option_string(self) -> str: ...

class SuffixedMultiWidget(forms.MultiWidget):
    suffixes: list[str]
    def suffixed(self, name: str, suffix: str) -> str: ...
    # Widget value and context can contain any data types
    def get_context(self, name: str, value: Any, attrs: dict[str, Any] | None) -> dict[str, Any]: ...
    # Returns list of any value types from widget data
    def value_from_datadict(self, data: Mapping[str, Any], files: Mapping[str, Any], name: str) -> list[Any]: ...
    # Widget data can contain any types
    def value_omitted_from_data(self, data: Mapping[str, Any], files: Mapping[str, Any], name: str) -> bool: ...
    def replace_name(self, output: str, index: int) -> str: ...
    # Decompresses any widget value into list of components
    def decompress(self, value: Any) -> list[Any] | None: ...

class RangeWidget(SuffixedMultiWidget):
    template_name: str
    suffixes: list[str]
    # Accepts any widget attribute types
    def __init__(self, attrs: dict[str, Any] | None = None) -> None: ...
    # Decompresses any range value into list components
    def decompress(self, value: Any) -> list[Any] | None: ...

class DateRangeWidget(RangeWidget):
    suffixes: list[str]

class LookupChoiceWidget(SuffixedMultiWidget):
    suffixes: list[str]
    # Decompresses any lookup choice value into components
    def decompress(self, value: Any) -> list[Any] | None: ...

class BooleanWidget(forms.Select):
    # Accepts any widget attribute types
    def __init__(self, attrs: dict[str, Any] | None = None) -> None: ...
    # Widget value and renderer can be any type
    def render(self, name: str, value: Any, attrs: dict[str, Any] | None = None, renderer: Any | None = None) -> SafeString: ...
    # Return value type depends on widget data
    def value_from_datadict(self, data: Mapping[str, Any], files: Mapping[str, Any], name: str) -> Any: ...

class BaseCSVWidget(forms.Widget):
    # Can be widget class or instance - __init__ converts to instance via instantiation or deepcopy
    surrogate: type[Any] = ...

    # CSV widget data can contain any types
    def value_from_datadict(self, data: Mapping[str, Any], files: Mapping[str, Any], name: str) -> list[str]: ...
    # Widget value and renderer can be any type
    def render(self, name: str, value: Any, attrs: dict[str, Any] | None = None, renderer: Any | None = None) -> SafeString: ...

class CSVWidget(BaseCSVWidget, forms.TextInput): ...

class QueryArrayWidget(BaseCSVWidget, forms.TextInput):
    # Query array widget data can contain any types
    def value_from_datadict(self, data: Mapping[str, Any], files: Mapping[str, Any], name: str) -> list[str]: ...
