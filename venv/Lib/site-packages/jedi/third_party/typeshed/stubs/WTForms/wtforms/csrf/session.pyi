from _typeshed import SupportsItemAccess
from datetime import datetime, timedelta
from typing import Any

from wtforms.csrf.core import CSRF, CSRFTokenField
from wtforms.form import BaseForm
from wtforms.meta import DefaultMeta

__all__ = ("SessionCSRF",)

class SessionCSRF(CSRF):
    TIME_FORMAT: str
    form_meta: DefaultMeta
    def generate_csrf_token(self, csrf_token_field: CSRFTokenField) -> str: ...
    def validate_csrf_token(self, form: BaseForm, field: CSRFTokenField) -> None: ...
    def now(self) -> datetime: ...
    @property
    def time_limit(self) -> timedelta: ...
    @property
    def session(self) -> SupportsItemAccess[str, Any]: ...
