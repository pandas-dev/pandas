from abc import abstractmethod
from collections.abc import Callable, Sequence
from typing import Any
from typing_extensions import Self

from wtforms.fields import HiddenField
from wtforms.fields.core import UnboundField, _Filter, _FormT, _Validator, _Widget
from wtforms.form import BaseForm
from wtforms.meta import DefaultMeta, _SupportsGettextAndNgettext

__all__ = ("CSRFTokenField", "CSRF")

class CSRFTokenField(HiddenField):
    current_token: str | None
    csrf_impl: CSRF
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
        *,
        csrf_impl: CSRF,
    ) -> None: ...

class CSRF:
    field_class: type[CSRFTokenField]
    def setup_form(self, form: BaseForm) -> list[tuple[str, UnboundField[Any]]]: ...
    @abstractmethod
    def generate_csrf_token(self, csrf_token_field: CSRFTokenField) -> str: ...
    @abstractmethod
    def validate_csrf_token(self, form: BaseForm, field: CSRFTokenField) -> None: ...
