import sys
from typing import Any, ClassVar, Final

import click
from click_web.web_click_types import EmailParamType, PasswordParamType, TextAreaParamType

class FieldId:
    SEPARATOR: ClassVar[str]
    command_index: int
    param_index: int
    param_type: str
    click_type: str
    nargs: int
    form_type: str
    name: str
    key: str

    def __init__(
        self,
        command_index: int,
        param_index: int,
        param_type: str,
        click_type: str,
        nargs: int,
        form_type: str,
        name: str,
        key: str | None = None,
    ) -> None: ...
    @classmethod
    def from_string(cls, field_info_as_string: str) -> FieldId: ...

class NotSupported(ValueError): ...

class BaseInput:
    param_type_cls: type[click.types.ParamType] | None
    ctx: click.Context
    param: click.Parameter
    command_index: int
    param_index: int
    def __init__(self, ctx: click.Context, param: click.Parameter, command_index: int, param_index: int) -> None: ...
    def is_supported(self) -> bool: ...
    @property
    def fields(self) -> dict[str, Any]: ...
    @property
    def type_attrs(self) -> dict[str, Any]: ...
    def _to_cmd_line_name(self, name: str) -> str: ...
    def _build_name(self, name: str) -> str: ...

class ChoiceInput(BaseInput):
    if sys.version_info >= (3, 10):
        param_type_cls: type[click.Choice[Any]]
    else:
        param_type_cls: type[click.Choice]

class FlagInput(BaseInput):
    param_type_cls: None

class IntInput(BaseInput):
    param_type_cls: type[click.types.IntParamType]

class FloatInput(BaseInput):
    param_type_cls: type[click.types.FloatParamType]

class FolderInput(BaseInput):
    param_type_cls: None

class FileInput(BaseInput):
    param_type_cls: None

class EmailInput(BaseInput):
    param_type_cls: type[EmailParamType]

class PasswordInput(BaseInput):
    param_type_cls: type[PasswordParamType]

class TextAreaInput(BaseInput):
    param_type_cls: type[TextAreaParamType]

class DefaultInput(BaseInput):
    param_type_cls: type[click.ParamType]

INPUT_TYPES: Final[list[type[BaseInput]]]
_DEFAULT_INPUT: Final[list[type[DefaultInput]]]

def get_input_field(ctx: click.Context, param: click.Parameter, command_index: int, param_index: int) -> dict[str, Any]: ...
