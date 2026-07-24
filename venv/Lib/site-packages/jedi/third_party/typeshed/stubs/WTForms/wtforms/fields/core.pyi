from builtins import type as _type  # type is being shadowed in Field
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Generic, Protocol, TypeVar, overload, type_check_only
from typing_extensions import Self, TypeAlias

from markupsafe import Markup
from wtforms.form import BaseForm
from wtforms.meta import DefaultMeta, _MultiDictLikeWithGetlist, _SupportsGettextAndNgettext

_FormT = TypeVar("_FormT", bound=BaseForm)
_FieldT = TypeVar("_FieldT", bound=Field)
_FormT_contra = TypeVar("_FormT_contra", bound=BaseForm, contravariant=True)
_FieldT_contra = TypeVar("_FieldT_contra", bound=Field, contravariant=True)
# It would be nice to annotate this as invariant, i.e. input type and output type
# needs to be the same, but it will probably be too annoying to use, for now we
# trust, that people won't use it to change the type of data in a field...
_Filter: TypeAlias = Callable[[Any], Any]

@type_check_only
class _Validator(Protocol[_FormT_contra, _FieldT_contra]):
    def __call__(self, form: _FormT_contra, field: _FieldT_contra, /) -> object: ...

@type_check_only
class _Widget(Protocol[_FieldT_contra]):
    def __call__(self, field: _FieldT_contra, **kwargs: Any) -> Markup: ...

class Field:
    errors: Sequence[str]
    process_errors: Sequence[str]
    raw_data: list[Any] | None
    object_data: Any
    data: Any
    validators: Sequence[_Validator[Any, Self]]
    # even though this could be None on the base class, this should
    # never actually be None in a real field
    widget: _Widget[Self]
    do_not_call_in_templates: bool
    meta: DefaultMeta
    default: Any | None
    description: str
    render_kw: dict[str, Any]
    filters: Sequence[_Filter]
    flags: Flags
    name: str
    short_name: str
    id: str
    type: str
    label: Label
    # technically this can return UnboundField, but that is not allowed
    # by type checkers, so we use a descriptor hack to get around this
    # limitation instead
    def __new__(cls, *args: Any, **kwargs: Any) -> Self: ...
    def __init__(
        self,
        label: str | None = None,
        # for tuple we can be a bit more type safe and only accept validators
        # that would work on this or a less specific field, but in general it
        # would be too annoying to restrict to Sequence[_Validator], since mypy
        # will infer a list of mixed validators as list[object], since that is
        # the common base class between all validators
        validators: tuple[_Validator[_FormT, Self], ...] | list[Any] | None = None,
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
    def __html__(self) -> str: ...
    def __call__(self, **kwargs: object) -> Markup: ...
    @classmethod
    def check_validators(cls, validators: Iterable[_Validator[_FormT, Self]] | None) -> None: ...
    def gettext(self, string: str) -> str: ...
    def ngettext(self, singular: str, plural: str, n: int) -> str: ...
    def validate(self, form: BaseForm, extra_validators: tuple[_Validator[_FormT, Self], ...] | list[Any] = ()) -> bool: ...
    def pre_validate(self, form: BaseForm) -> None: ...
    def post_validate(self, form: BaseForm, validation_stopped: bool) -> None: ...
    def process(
        self, formdata: _MultiDictLikeWithGetlist | None, data: Any = ..., extra_filters: Sequence[_Filter] | None = None
    ) -> None: ...
    def process_data(self, value: Any) -> None: ...
    def process_formdata(self, valuelist: list[Any]) -> None: ...
    def populate_obj(self, obj: object, name: str) -> None: ...

    # this is a workaround for what is essentially illegal in static type checking
    # Field.__new__ would return an UnboundField, unless the _form parameter is
    # specified. We can't really work around it by making UnboundField a subclass
    # of Field, since all subclasses of Field still need to return an UnboundField
    # and we can't expect third parties to add a __new__ method to every field
    # they define...
    # This workaround only works for Form, not BaseForm, but we take what we can get
    # BaseForm shouldn't really be used anyways
    @overload
    def __get__(self, obj: None, owner: _type[object] | None = None) -> UnboundField[Self]: ...
    @overload
    def __get__(self, obj: object, owner: _type[object] | None = None) -> Self: ...

class UnboundField(Generic[_FieldT]):
    creation_counter: int
    field_class: type[_FieldT]
    name: str | None
    args: tuple[Any, ...]
    kwargs: dict[str, Any]
    def __init__(self, field_class: type[_FieldT], *args: object, name: str | None = None, **kwargs: object) -> None: ...
    def bind(
        self,
        form: BaseForm,
        name: str,
        prefix: str = "",
        translations: _SupportsGettextAndNgettext | None = None,
        **kwargs: object,
    ) -> _FieldT: ...

class Flags:
    # the API for this is a bit loosey goosey, the intention probably
    # was that the values should always be boolean, but __contains__
    # just returns the same thing as __getattr__ and in the widgets
    # there are fields that could accept numeric values from Flags
    def __getattr__(self, name: str) -> Any | None: ...
    def __setattr__(self, name: str, value: object) -> None: ...
    def __delattr__(self, name: str) -> None: ...
    def __contains__(self, name: str) -> Any | None: ...

class Label:
    field_id: str
    text: str
    def __init__(self, field_id: str, text: str) -> None: ...
    def __html__(self) -> str: ...
    def __call__(self, text: str | None = None, **kwargs: Any) -> Markup: ...
