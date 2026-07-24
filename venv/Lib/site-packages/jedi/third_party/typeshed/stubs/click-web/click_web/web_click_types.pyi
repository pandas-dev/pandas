import re
from typing import ClassVar, TypeVar

import click

_T = TypeVar("_T")

class EmailParamType(click.ParamType):
    EMAIL_REGEX: ClassVar[re.Pattern[str]]
    def convert(self, value: str, param: click.Parameter | None, ctx: click.Context | None) -> str: ...

class PasswordParamType(click.ParamType):
    def convert(self, value: _T, param: click.Parameter | None, ctx: click.Context | None) -> _T: ...

class TextAreaParamType(click.ParamType):
    def convert(self, value: _T, param: click.Parameter | None, ctx: click.Context | None) -> _T: ...

EMAIL_TYPE: EmailParamType
PASSWORD_TYPE: PasswordParamType
TEXTAREA_TYPE: TextAreaParamType
