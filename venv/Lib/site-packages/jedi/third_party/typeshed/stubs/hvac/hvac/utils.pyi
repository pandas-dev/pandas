from collections.abc import Callable, Iterable, Mapping
from typing import Any, NoReturn, TypedDict, TypeVar, type_check_only
from typing_extensions import NotRequired

@type_check_only
class _DeprecateProperty(TypedDict):
    to_be_removed_in_version: str
    client_property: str
    new_property: NotRequired[str]

_T = TypeVar("_T")
_K = TypeVar("_K")
_V = TypeVar("_V")

def raise_for_error(
    method: str,
    url: str,
    status_code: int,
    message: str | None = None,
    errors: Iterable[str] | str | None = None,
    text: str | None = None,
    json: object | None = None,
) -> NoReturn: ...
def aliased_parameter(
    name: str, *aliases: str, removed_in_version: str | None, position: int | None = None, raise_on_multiple: bool = True
) -> Callable[..., Any]: ...
def generate_parameter_deprecation_message(
    to_be_removed_in_version: str, old_parameter_name: str, new_parameter_name: str | None = None, extra_notes: str | None = None
) -> str: ...
def generate_method_deprecation_message(
    to_be_removed_in_version: str, old_method_name: str, method_name: str | None = None, module_name: str | None = None
) -> str: ...
def generate_property_deprecation_message(
    to_be_removed_in_version: str, old_name: str, new_name: str, new_attribute: str, module_name: str = "Client"
) -> str: ...
def getattr_with_deprecated_properties(obj: object, item: str, deprecated_properties: dict[str, _DeprecateProperty]) -> Any: ...
def deprecated_method(to_be_removed_in_version: str, new_method: Callable[..., Any] | None = None) -> Callable[..., Any]: ...
def validate_list_of_strings_param(param_name: str, param_argument: Iterable[Any] | str) -> None: ...
def list_to_comma_delimited(list_param: Iterable[str] | None) -> str: ...
def get_token_from_env() -> str | None: ...
def comma_delimited_to_list(list_param: Iterable[_T]) -> Iterable[_T]: ...

# the docstring states that this function returns a bool, but the code does not return anything
def validate_pem_format(param_name: str, param_argument: str) -> None: ...
def remove_nones(params: Mapping[_K, _V | None]) -> Mapping[_K, _V]: ...
def format_url(
    format_str: str, *args: object, **kwargs: object
) -> str: ...  # values are passed to builtins.str, which takes an object type
