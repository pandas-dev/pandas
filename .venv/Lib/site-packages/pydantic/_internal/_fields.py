"""Private logic related to fields (the `Field()` function and `FieldInfo` class), and arguments to `Annotated`."""

from __future__ import annotations as _annotations

import dataclasses
import warnings
from collections.abc import Mapping
from copy import copy
from functools import cache
from inspect import Parameter, ismethoddescriptor, signature
from re import Pattern
from typing import TYPE_CHECKING, Any, Callable, TypeVar

from pydantic_core import PydanticUndefined
from typing_extensions import TypeIs, get_origin
from typing_inspection import typing_objects
from typing_inspection.introspection import AnnotationSource

from pydantic import PydanticDeprecatedSince211
from pydantic.errors import PydanticUserError

from . import _generics, _typing_extra
from ._config import ConfigWrapper
from ._docs_extraction import extract_docstrings_from_cls
from ._import_utils import import_cached_base_model, import_cached_field_info
from ._namespace_utils import NsResolver
from ._repr import Representation
from ._utils import can_be_positional

if TYPE_CHECKING:
    from annotated_types import BaseMetadata

    from ..fields import FieldInfo
    from ..main import BaseModel
    from ._dataclasses import PydanticDataclass, StandardDataclass
    from ._decorators import DecoratorInfos


class PydanticMetadata(Representation):
    """Base class for annotation markers like `Strict`."""

    __slots__ = ()


def pydantic_general_metadata(**metadata: Any) -> BaseMetadata:
    """Create a new `_PydanticGeneralMetadata` class with the given metadata.

    Args:
        **metadata: The metadata to add.

    Returns:
        The new `_PydanticGeneralMetadata` class.
    """
    return _general_metadata_cls()(metadata)  # type: ignore


@cache
def _general_metadata_cls() -> type[BaseMetadata]:
    """Do it this way to avoid importing `annotated_types` at import time."""
    from annotated_types import BaseMetadata

    class _PydanticGeneralMetadata(PydanticMetadata, BaseMetadata):
        """Pydantic general metadata like `max_digits`."""

        def __init__(self, metadata: Any):
            self.__dict__ = metadata

    return _PydanticGeneralMetadata  # type: ignore


def _update_fields_from_docstrings(cls: type[Any], fields: dict[str, FieldInfo], use_inspect: bool = False) -> None:
    fields_docs = extract_docstrings_from_cls(cls, use_inspect=use_inspect)
    for ann_name, field_info in fields.items():
        if field_info.description is None and ann_name in fields_docs:
            field_info.description = fields_docs[ann_name]


def collect_model_fields(  # noqa: C901
    cls: type[BaseModel],
    config_wrapper: ConfigWrapper,
    ns_resolver: NsResolver | None,
    *,
    typevars_map: Mapping[TypeVar, Any] | None = None,
) -> tuple[dict[str, FieldInfo], set[str]]:
    """Collect the fields and class variables names of a nascent Pydantic model.

    The fields collection process is *lenient*, meaning it won't error if string annotations
    fail to evaluate. If this happens, the original annotation (and assigned value, if any)
    is stored on the created `FieldInfo` instance.

    The `rebuild_model_fields()` should be called at a later point (e.g. when rebuilding the model),
    and will make use of these stored attributes.

    Args:
        cls: BaseModel or dataclass.
        config_wrapper: The config wrapper instance.
        ns_resolver: Namespace resolver to use when getting model annotations.
        typevars_map: A dictionary mapping type variables to their concrete types.

    Returns:
        A two-tuple containing model fields and class variables names.

    Raises:
        NameError:
            - If there is a conflict between a field name and protected namespaces.
            - If there is a field other than `root` in `RootModel`.
            - If a field shadows an attribute in the parent model.
    """
    BaseModel = import_cached_base_model()
    FieldInfo_ = import_cached_field_info()

    bases = cls.__bases__
    parent_fields_lookup: dict[str, FieldInfo] = {}
    for base in reversed(bases):
        if model_fields := getattr(base, '__pydantic_fields__', None):
            parent_fields_lookup.update(model_fields)

    type_hints = _typing_extra.get_model_type_hints(cls, ns_resolver=ns_resolver)

    # https://docs.python.org/3/howto/annotations.html#accessing-the-annotations-dict-of-an-object-in-python-3-9-and-older
    # annotations is only used for finding fields in parent classes
    annotations = cls.__dict__.get('__annotations__', {})
    fields: dict[str, FieldInfo] = {}

    class_vars: set[str] = set()
    for ann_name, (ann_type, evaluated) in type_hints.items():
        if ann_name == 'model_config':
            # We never want to treat `model_config` as a field
            # Note: we may need to change this logic if/when we introduce a `BareModel` class with no
            # protected namespaces (where `model_config` might be allowed as a field name)
            continue

        for protected_namespace in config_wrapper.protected_namespaces:
            ns_violation: bool = False
            if isinstance(protected_namespace, Pattern):
                ns_violation = protected_namespace.match(ann_name) is not None
            elif isinstance(protected_namespace, str):
                ns_violation = ann_name.startswith(protected_namespace)

            if ns_violation:
                for b in bases:
                    if hasattr(b, ann_name):
                        if not (issubclass(b, BaseModel) and ann_name in getattr(b, '__pydantic_fields__', {})):
                            raise NameError(
                                f'Field "{ann_name}" conflicts with member {getattr(b, ann_name)}'
                                f' of protected namespace "{protected_namespace}".'
                            )
                else:
                    valid_namespaces = ()
                    for pn in config_wrapper.protected_namespaces:
                        if isinstance(pn, Pattern):
                            if not pn.match(ann_name):
                                valid_namespaces += (f're.compile({pn.pattern})',)
                        else:
                            if not ann_name.startswith(pn):
                                valid_namespaces += (pn,)

                    warnings.warn(
                        f'Field "{ann_name}" in {cls.__name__} has conflict with protected namespace "{protected_namespace}".'
                        '\n\nYou may be able to resolve this warning by setting'
                        f" `model_config['protected_namespaces'] = {valid_namespaces}`.",
                        UserWarning,
                    )
        if _typing_extra.is_classvar_annotation(ann_type):
            class_vars.add(ann_name)
            continue

        assigned_value = getattr(cls, ann_name, PydanticUndefined)

        if not is_valid_field_name(ann_name):
            continue
        if cls.__pydantic_root_model__ and ann_name != 'root':
            raise NameError(
                f"Unexpected field with name {ann_name!r}; only 'root' is allowed as a field of a `RootModel`"
            )

        # when building a generic model with `MyModel[int]`, the generic_origin check makes sure we don't get
        # "... shadows an attribute" warnings
        generic_origin = getattr(cls, '__pydantic_generic_metadata__', {}).get('origin')
        for base in bases:
            dataclass_fields = {
                field.name for field in (dataclasses.fields(base) if dataclasses.is_dataclass(base) else ())
            }
            if hasattr(base, ann_name):
                if base is generic_origin:
                    # Don't warn when "shadowing" of attributes in parametrized generics
                    continue

                if ann_name in dataclass_fields:
                    # Don't warn when inheriting stdlib dataclasses whose fields are "shadowed" by defaults being set
                    # on the class instance.
                    continue

                if ann_name not in annotations:
                    # Don't warn when a field exists in a parent class but has not been defined in the current class
                    continue

                warnings.warn(
                    f'Field name "{ann_name}" in "{cls.__qualname__}" shadows an attribute in parent '
                    f'"{base.__qualname__}"',
                    UserWarning,
                )

        if assigned_value is PydanticUndefined:  # no assignment, just a plain annotation
            if ann_name in annotations or ann_name not in parent_fields_lookup:
                # field is either:
                # - present in the current model's annotations (and *not* from parent classes)
                # - not found on any base classes; this seems to be caused by fields bot getting
                #   generated due to models not being fully defined while initializing recursive models.
                #   Nothing stops us from just creating a `FieldInfo` for this type hint, so we do this.
                field_info = FieldInfo_.from_annotation(ann_type, _source=AnnotationSource.CLASS)
                if not evaluated:
                    field_info._complete = False
                    # Store the original annotation that should be used to rebuild
                    # the field info later:
                    field_info._original_annotation = ann_type
            else:
                # The field was present on one of the (possibly multiple) base classes
                # copy the field to make sure typevar substitutions don't cause issues with the base classes
                field_info = copy(parent_fields_lookup[ann_name])

        else:  # An assigned value is present (either the default value, or a `Field()` function)
            _warn_on_nested_alias_in_annotation(ann_type, ann_name)
            if isinstance(assigned_value, FieldInfo_) and ismethoddescriptor(assigned_value.default):
                # `assigned_value` was fetched using `getattr`, which triggers a call to `__get__`
                # for descriptors, so we do the same if the `= field(default=...)` form is used.
                # Note that we only do this for method descriptors for now, we might want to
                # extend this to any descriptor in the future (by simply checking for
                # `hasattr(assigned_value.default, '__get__')`).
                assigned_value.default = assigned_value.default.__get__(None, cls)

            # The `from_annotated_attribute()` call below mutates the assigned `Field()`, so make a copy:
            original_assignment = (
                assigned_value._copy() if not evaluated and isinstance(assigned_value, FieldInfo_) else assigned_value
            )

            field_info = FieldInfo_.from_annotated_attribute(ann_type, assigned_value, _source=AnnotationSource.CLASS)
            # Store the original annotation and assignment value that should be used to rebuild the field info later.
            # Note that the assignment is always stored as the annotation might contain a type var that is later
            #  parameterized with an unknown forward reference (and we'll need it to rebuild the field info):
            field_info._original_assignment = original_assignment
            if not evaluated:
                field_info._complete = False
                field_info._original_annotation = ann_type
            elif 'final' in field_info._qualifiers and not field_info.is_required():
                warnings.warn(
                    f'Annotation {ann_name!r} is marked as final and has a default value. Pydantic treats {ann_name!r} as a '
                    'class variable, but it will be considered as a normal field in V3 to be aligned with dataclasses. If you '
                    f'still want {ann_name!r} to be considered as a class variable, annotate it as: `ClassVar[<type>] = <default>.`',
                    category=PydanticDeprecatedSince211,
                    # Incorrect when `create_model` is used, but the chance that final with a default is used is low in that case:
                    stacklevel=4,
                )
                class_vars.add(ann_name)
                continue

            # attributes which are fields are removed from the class namespace:
            # 1. To match the behaviour of annotation-only fields
            # 2. To avoid false positives in the NameError check above
            try:
                delattr(cls, ann_name)
            except AttributeError:
                pass  # indicates the attribute was on a parent class

        # Use cls.__dict__['__pydantic_decorators__'] instead of cls.__pydantic_decorators__
        # to make sure the decorators have already been built for this exact class
        decorators: DecoratorInfos = cls.__dict__['__pydantic_decorators__']
        if ann_name in decorators.computed_fields:
            raise TypeError(
                f'Field {ann_name!r} of class {cls.__name__!r} overrides symbol of same name in a parent class. '
                'This override with a computed_field is incompatible.'
            )
        fields[ann_name] = field_info

    if typevars_map:
        for field in fields.values():
            if field._complete:
                field.apply_typevars_map(typevars_map)

    if config_wrapper.use_attribute_docstrings:
        _update_fields_from_docstrings(cls, fields)
    return fields, class_vars


def _warn_on_nested_alias_in_annotation(ann_type: type[Any], ann_name: str) -> None:
    FieldInfo = import_cached_field_info()

    args = getattr(ann_type, '__args__', None)
    if args:
        for anno_arg in args:
            if typing_objects.is_annotated(get_origin(anno_arg)):
                for anno_type_arg in _typing_extra.get_args(anno_arg):
                    if isinstance(anno_type_arg, FieldInfo) and anno_type_arg.alias is not None:
                        warnings.warn(
                            f'`alias` specification on field "{ann_name}" must be set on outermost annotation to take effect.',
                            UserWarning,
                        )
                        return


def rebuild_model_fields(
    cls: type[BaseModel],
    *,
    ns_resolver: NsResolver,
    typevars_map: Mapping[TypeVar, Any],
) -> dict[str, FieldInfo]:
    """Rebuild the (already present) model fields by trying to reevaluate annotations.

    This function should be called whenever a model with incomplete fields is encountered.

    Raises:
        NameError: If one of the annotations failed to evaluate.

    Note:
        This function *doesn't* mutate the model fields in place, as it can be called during
        schema generation, where you don't want to mutate other model's fields.
    """
    FieldInfo_ = import_cached_field_info()

    rebuilt_fields: dict[str, FieldInfo] = {}
    with ns_resolver.push(cls):
        for f_name, field_info in cls.__pydantic_fields__.items():
            if field_info._complete:
                rebuilt_fields[f_name] = field_info
            else:
                existing_desc = field_info.description
                ann = _typing_extra.eval_type(
                    field_info._original_annotation,
                    *ns_resolver.types_namespace,
                )
                ann = _generics.replace_types(ann, typevars_map)

                if (assign := field_info._original_assignment) is PydanticUndefined:
                    new_field = FieldInfo_.from_annotation(ann, _source=AnnotationSource.CLASS)
                else:
                    new_field = FieldInfo_.from_annotated_attribute(ann, assign, _source=AnnotationSource.CLASS)
                # The description might come from the docstring if `use_attribute_docstrings` was `True`:
                new_field.description = new_field.description if new_field.description is not None else existing_desc
                rebuilt_fields[f_name] = new_field

    return rebuilt_fields


def collect_dataclass_fields(
    cls: type[StandardDataclass],
    *,
    ns_resolver: NsResolver | None = None,
    typevars_map: dict[Any, Any] | None = None,
    config_wrapper: ConfigWrapper | None = None,
) -> dict[str, FieldInfo]:
    """Collect the fields of a dataclass.

    Args:
        cls: dataclass.
        ns_resolver: Namespace resolver to use when getting dataclass annotations.
            Defaults to an empty instance.
        typevars_map: A dictionary mapping type variables to their concrete types.
        config_wrapper: The config wrapper instance.

    Returns:
        The dataclass fields.
    """
    FieldInfo_ = import_cached_field_info()

    fields: dict[str, FieldInfo] = {}
    ns_resolver = ns_resolver or NsResolver()
    dataclass_fields = cls.__dataclass_fields__

    # The logic here is similar to `_typing_extra.get_cls_type_hints`,
    # although we do it manually as stdlib dataclasses already have annotations
    # collected in each class:
    for base in reversed(cls.__mro__):
        if not dataclasses.is_dataclass(base):
            continue

        with ns_resolver.push(base):
            for ann_name, dataclass_field in dataclass_fields.items():
                if ann_name not in base.__dict__.get('__annotations__', {}):
                    # `__dataclass_fields__`contains every field, even the ones from base classes.
                    # Only collect the ones defined on `base`.
                    continue

                globalns, localns = ns_resolver.types_namespace
                ann_type, evaluated = _typing_extra.try_eval_type(dataclass_field.type, globalns, localns)

                if _typing_extra.is_classvar_annotation(ann_type):
                    continue

                if (
                    not dataclass_field.init
                    and dataclass_field.default is dataclasses.MISSING
                    and dataclass_field.default_factory is dataclasses.MISSING
                ):
                    # TODO: We should probably do something with this so that validate_assignment behaves properly
                    #   Issue: https://github.com/pydantic/pydantic/issues/5470
                    continue

                if isinstance(dataclass_field.default, FieldInfo_):
                    if dataclass_field.default.init_var:
                        if dataclass_field.default.init is False:
                            raise PydanticUserError(
                                f'Dataclass field {ann_name} has init=False and init_var=True, but these are mutually exclusive.',
                                code='clashing-init-and-init-var',
                            )

                        # TODO: same note as above re validate_assignment
                        continue
                    field_info = FieldInfo_.from_annotated_attribute(
                        ann_type, dataclass_field.default, _source=AnnotationSource.DATACLASS
                    )
                    field_info._original_assignment = dataclass_field.default
                else:
                    field_info = FieldInfo_.from_annotated_attribute(
                        ann_type, dataclass_field, _source=AnnotationSource.DATACLASS
                    )
                    field_info._original_assignment = dataclass_field

                if not evaluated:
                    field_info._complete = False
                    field_info._original_annotation = ann_type

                fields[ann_name] = field_info

                if field_info.default is not PydanticUndefined and isinstance(
                    getattr(cls, ann_name, field_info), FieldInfo_
                ):
                    # We need this to fix the default when the "default" from __dataclass_fields__ is a pydantic.FieldInfo
                    setattr(cls, ann_name, field_info.default)

    if typevars_map:
        for field in fields.values():
            # We don't pass any ns, as `field.annotation`
            # was already evaluated. TODO: is this method relevant?
            # Can't we juste use `_generics.replace_types`?
            field.apply_typevars_map(typevars_map)

    if config_wrapper is not None and config_wrapper.use_attribute_docstrings:
        _update_fields_from_docstrings(
            cls,
            fields,
            # We can't rely on the (more reliable) frame inspection method
            # for stdlib dataclasses:
            use_inspect=not hasattr(cls, '__is_pydantic_dataclass__'),
        )

    return fields


def rebuild_dataclass_fields(
    cls: type[PydanticDataclass],
    *,
    config_wrapper: ConfigWrapper,
    ns_resolver: NsResolver,
    typevars_map: Mapping[TypeVar, Any],
) -> dict[str, FieldInfo]:
    """Rebuild the (already present) dataclass fields by trying to reevaluate annotations.

    This function should be called whenever a dataclass with incomplete fields is encountered.

    Raises:
        NameError: If one of the annotations failed to evaluate.

    Note:
        This function *doesn't* mutate the dataclass fields in place, as it can be called during
        schema generation, where you don't want to mutate other dataclass's fields.
    """
    FieldInfo_ = import_cached_field_info()

    rebuilt_fields: dict[str, FieldInfo] = {}
    with ns_resolver.push(cls):
        for f_name, field_info in cls.__pydantic_fields__.items():
            if field_info._complete:
                rebuilt_fields[f_name] = field_info
            else:
                existing_desc = field_info.description
                ann = _typing_extra.eval_type(
                    field_info._original_annotation,
                    *ns_resolver.types_namespace,
                )
                ann = _generics.replace_types(ann, typevars_map)
                new_field = FieldInfo_.from_annotated_attribute(
                    ann,
                    field_info._original_assignment,
                    _source=AnnotationSource.DATACLASS,
                )

                # The description might come from the docstring if `use_attribute_docstrings` was `True`:
                new_field.description = new_field.description if new_field.description is not None else existing_desc
                rebuilt_fields[f_name] = new_field

    return rebuilt_fields


def is_valid_field_name(name: str) -> bool:
    return not name.startswith('_')


def is_valid_privateattr_name(name: str) -> bool:
    return name.startswith('_') and not name.startswith('__')


def takes_validated_data_argument(
    default_factory: Callable[[], Any] | Callable[[dict[str, Any]], Any],
) -> TypeIs[Callable[[dict[str, Any]], Any]]:
    """Whether the provided default factory callable has a validated data parameter."""
    try:
        sig = signature(default_factory)
    except (ValueError, TypeError):
        # `inspect.signature` might not be able to infer a signature, e.g. with C objects.
        # In this case, we assume no data argument is present:
        return False

    parameters = list(sig.parameters.values())

    return len(parameters) == 1 and can_be_positional(parameters[0]) and parameters[0].default is Parameter.empty
