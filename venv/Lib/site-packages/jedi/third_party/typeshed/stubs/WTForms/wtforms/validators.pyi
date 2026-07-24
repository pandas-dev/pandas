from collections.abc import Callable, Collection, Iterable
from decimal import Decimal
from re import Match, Pattern
from typing import Any, TypeVar, overload

from wtforms.fields import Field, StringField
from wtforms.form import BaseForm

__all__ = (
    "DataRequired",
    "data_required",
    "Email",
    "email",
    "EqualTo",
    "equal_to",
    "IPAddress",
    "ip_address",
    "InputRequired",
    "input_required",
    "Length",
    "length",
    "NumberRange",
    "number_range",
    "Optional",
    "optional",
    "Regexp",
    "regexp",
    "URL",
    "url",
    "AnyOf",
    "any_of",
    "NoneOf",
    "none_of",
    "MacAddress",
    "mac_address",
    "UUID",
    "ValidationError",
    "StopValidation",
    "readonly",
    "ReadOnly",
    "disabled",
    "Disabled",
)

_ValuesT_contra = TypeVar("_ValuesT_contra", bound=Collection[Any], contravariant=True)

class ValidationError(ValueError):
    def __init__(self, message: str = "", *args: object) -> None: ...

class StopValidation(Exception):
    def __init__(self, message: str = "", *args: object) -> None: ...

class EqualTo:
    fieldname: str
    message: str | None
    def __init__(self, fieldname: str, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: Field) -> None: ...

class Length:
    min: int
    max: int
    message: str | None
    field_flags: dict[str, Any]
    def __init__(self, min: int = -1, max: int = -1, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: StringField) -> None: ...

class NumberRange:
    min: float | Decimal | None
    max: float | Decimal | None
    message: str | None
    field_flags: dict[str, Any]
    def __init__(
        self, min: float | Decimal | None = None, max: float | Decimal | None = None, message: str | None = None
    ) -> None: ...
    # any numeric field will work, for now we don't try to use a union
    # to restrict to the defined numeric fields, since user-defined fields
    # will likely not use a common base class, just like the existing
    # numeric fields
    def __call__(self, form: BaseForm, field: Field) -> None: ...

class Optional:
    string_check: Callable[[str], bool]
    field_flags: dict[str, Any]
    def __init__(self, strip_whitespace: bool = True) -> None: ...
    def __call__(self, form: BaseForm, field: Field) -> None: ...

class DataRequired:
    message: str | None
    field_flags: dict[str, Any]
    def __init__(self, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: Field) -> None: ...

class InputRequired:
    message: str | None
    field_flags: dict[str, Any]
    def __init__(self, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: Field) -> None: ...

class Regexp:
    regex: Pattern[str]
    message: str | None
    def __init__(self, regex: str | Pattern[str], flags: int = 0, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: StringField, message: str | None = None) -> Match[str]: ...

class Email:
    message: str | None
    granular_message: bool
    check_deliverability: bool
    allow_smtputf8: bool
    allow_empty_local: bool
    def __init__(
        self,
        message: str | None = None,
        granular_message: bool = False,
        check_deliverability: bool = False,
        allow_smtputf8: bool = True,
        allow_empty_local: bool = False,
    ) -> None: ...
    def __call__(self, form: BaseForm, field: StringField) -> None: ...

class IPAddress:
    ipv4: bool
    ipv6: bool
    message: str | None
    def __init__(self, ipv4: bool = True, ipv6: bool = False, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: StringField) -> None: ...
    @classmethod
    def check_ipv4(cls, value: str | None) -> bool: ...
    @classmethod
    def check_ipv6(cls, value: str | None) -> bool: ...

class MacAddress(Regexp):
    def __init__(self, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: StringField) -> None: ...  # type: ignore[override]

class URL(Regexp):
    validate_hostname: HostnameValidation
    def __init__(self, require_tld: bool = True, allow_ip: bool = True, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: StringField) -> None: ...  # type: ignore[override]

class UUID:
    message: str | None
    def __init__(self, message: str | None = None) -> None: ...
    def __call__(self, form: BaseForm, field: StringField) -> None: ...

class AnyOf:
    values: Collection[Any]
    message: str | None
    values_formatter: Callable[[Any], str]
    @overload
    def __init__(self, values: Collection[Any], message: str | None = None, values_formatter: None = None) -> None: ...
    @overload
    def __init__(
        self, values: _ValuesT_contra, message: str | None, values_formatter: Callable[[_ValuesT_contra], str]
    ) -> None: ...
    @overload
    def __init__(
        self, values: _ValuesT_contra, message: str | None = None, *, values_formatter: Callable[[_ValuesT_contra], str]
    ) -> None: ...
    def __call__(self, form: BaseForm, field: Field) -> None: ...
    @staticmethod
    def default_values_formatter(values: Iterable[object]) -> str: ...

class NoneOf:
    values: Collection[Any]
    message: str | None
    values_formatter: Callable[[Any], str]
    @overload
    def __init__(self, values: Collection[Any], message: str | None = None, values_formatter: None = None) -> None: ...
    @overload
    def __init__(
        self, values: _ValuesT_contra, message: str | None, values_formatter: Callable[[_ValuesT_contra], str]
    ) -> None: ...
    @overload
    def __init__(
        self, values: _ValuesT_contra, message: str | None = None, *, values_formatter: Callable[[_ValuesT_contra], str]
    ) -> None: ...
    def __call__(self, form: BaseForm, field: Field) -> None: ...
    @staticmethod
    def default_values_formatter(v: Iterable[object]) -> str: ...

class HostnameValidation:
    hostname_part: Pattern[str]
    tld_part: Pattern[str]
    require_tld: bool
    allow_ip: bool
    def __init__(self, require_tld: bool = True, allow_ip: bool = False) -> None: ...
    def __call__(self, hostname: str) -> bool: ...

class ReadOnly:
    def __call__(self, form: BaseForm, field: Field) -> None: ...

class Disabled:
    def __call__(self, form: BaseForm, field: Field) -> None: ...

email = Email
equal_to = EqualTo
ip_address = IPAddress
mac_address = MacAddress
length = Length
number_range = NumberRange
optional = Optional
input_required = InputRequired
data_required = DataRequired
regexp = Regexp
url = URL
any_of = AnyOf
none_of = NoneOf
readonly = ReadOnly
disabled = Disabled
