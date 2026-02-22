"""Object loader for autodoc"""

from __future__ import annotations

import re
from inspect import Parameter
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, NewType, TypeVar

from sphinx.ext.autodoc._dynamic._docstrings import (
    _docstring_lines_for_props,
    _get_docstring_lines,
)
from sphinx.ext.autodoc._dynamic._importer import _import_object
from sphinx.ext.autodoc._dynamic._mock import ismock
from sphinx.ext.autodoc._dynamic._signatures import _format_signatures
from sphinx.ext.autodoc._dynamic._type_comments import (
    _ensure_annotations_from_type_comments,
    _update_annotations_using_type_comments,
)
from sphinx.ext.autodoc._names import _parse_name
from sphinx.ext.autodoc._property_types import (
    _AssignStatementProperties,
    _ClassDefProperties,
    _FunctionDefProperties,
    _ItemProperties,
    _ModuleProperties,
    _TypeStatementProperties,
)
from sphinx.ext.autodoc._sentinels import (
    RUNTIME_INSTANCE_ATTRIBUTE,
    SLOTS_ATTR,
    UNINITIALIZED_ATTR,
)
from sphinx.ext.autodoc._shared import LOGGER, _get_render_mode
from sphinx.locale import __
from sphinx.util import inspect
from sphinx.util.inspect import safe_getattr
from sphinx.util.typing import get_type_hints, restify, stringify_annotation

if TYPE_CHECKING:
    from collections.abc import Mapping, MutableSet, Sequence
    from typing import Any

    from sphinx.environment import _CurrentDocument
    from sphinx.events import EventManager
    from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
    from sphinx.ext.autodoc._dynamic._importer import _ImportedObject
    from sphinx.ext.autodoc._property_types import _AutodocFuncProperty, _AutodocObjType
    from sphinx.ext.autodoc._shared import _AttrGetter, _AutodocConfig


_hide_value_re = re.compile(r'^:meta \s*hide-value:( +|$)')


def _load_object_by_name(
    *,
    name: str,
    objtype: _AutodocObjType,
    current_document: _CurrentDocument,
    config: _AutodocConfig,
    events: EventManager,
    get_attr: _AttrGetter,
    options: _AutoDocumenterOptions,
    parent_modname: str | None = None,
    ref_context: Mapping[str, str | None],
    reread_always: MutableSet[str],
) -> _ItemProperties | None:
    """Import and load the object given by *name*."""
    parsed = _parse_name(
        name=name,
        objtype=objtype,
        current_document=current_document,
        ref_context=ref_context,
    )
    if parsed is None:
        return None
    module_name, parts, args, retann = parsed

    # Import the module and get the object to document
    im = _import_object(
        module_name=module_name,
        obj_path=parts,
        mock_imports=config.autodoc_mock_imports,
        get_attr=get_attr,
        obj_type=objtype,
        type_aliases=config.autodoc_type_aliases,
    )
    if im is None:
        # See BuildEnvironment.note_reread()
        reread_always.add(current_document.docname)
        return None

    # Assemble object properties from the imported object.
    parent = im.parent
    props = _make_props_from_imported_object(
        im,
        config=config,
        events=events,
        get_attr=get_attr,
        module_name=module_name,
        objtype=objtype,
        parts=parts,
    )
    if props is None:
        return None

    if options.class_doc_from is not None:
        class_doc_from = options.class_doc_from
    else:
        class_doc_from = config.autoclass_content

    docstrings = _get_docstring_lines(
        props,
        class_doc_from=class_doc_from,
        get_attr=get_attr,
        inherit_docstrings=config.autodoc_inherit_docstrings,
        parent=parent,
        tab_width=options._tab_width,
    )
    if docstrings:
        for docstring_lines in docstrings:
            for line in docstring_lines:
                if _hide_value_re.match(line):
                    props._docstrings_has_hide_value = True
                    break

    # format the object's signature, if any
    try:
        signatures = _format_signatures(
            args=args,
            retann=retann,
            autodoc_annotations=current_document.autodoc_annotations,
            config=config,
            docstrings=docstrings,
            events=events,
            get_attr=get_attr,
            parent=parent,
            options=options,
            props=props,
        )
    except Exception as exc:
        msg = __('error while formatting signature for %s: %s')
        LOGGER.warning(msg, props.full_name, exc, type='autodoc')
        return None
    props.signatures = tuple(
        f'{args} -> {retann}' if retann else str(args) for args, retann in signatures
    )

    props.docstring_lines = _docstring_lines_for_props(
        docstrings,
        props=props,
        parent_modname=parent_modname,
        events=events,
        options=options,
    )

    return props


def _make_props_from_imported_object(
    im: _ImportedObject,
    *,
    config: _AutodocConfig,
    events: EventManager,
    get_attr: _AttrGetter,
    module_name: str,
    objtype: _AutodocObjType,
    parts: tuple[str, ...],
) -> _ItemProperties | None:
    parent = im.parent
    object_name = im.object_name
    obj = im.obj
    obj_properties: set[_AutodocFuncProperty] = set()
    render_mode = _get_render_mode(config.autodoc_typehints_format)

    if objtype == 'module':
        try:
            mod_origin = im.module.__spec__.origin  # type: ignore[union-attr]
        except AttributeError:
            file_path = None
        else:
            file_path = Path(mod_origin) if mod_origin is not None else None

        mod_all = safe_getattr(obj, '__all__', None)
        if isinstance(mod_all, (list, tuple)) and all(
            isinstance(e, str) for e in mod_all
        ):
            mod_all = tuple(mod_all)
        elif mod_all is not None:
            # Invalid __all__ found.
            msg = __('Ignoring invalid __all__ in module %s: %r')
            LOGGER.warning(msg, module_name, mod_all, type='autodoc')
            mod_all = None

        return _ModuleProperties(
            obj_type=objtype,
            module_name=module_name,
            docstring_lines=(),
            file_path=file_path,
            all=mod_all,
            _obj=obj,
            _obj___module__=obj.__name__,
        )

    if objtype in {'class', 'exception'}:
        if isinstance(obj, (NewType, TypeVar)):
            obj_module_name = getattr(obj, '__module__', module_name)
            if obj_module_name != module_name and module_name.startswith(
                obj_module_name
            ):
                bases = module_name[len(obj_module_name) :].strip('.').split('.')
                parts = tuple(bases) + parts
                module_name = obj_module_name

        if orig_bases := inspect.getorigbases(obj):
            # A subclass of generic types
            # refs: PEP-560 <https://peps.python.org/pep-0560/>
            obj_bases = list(orig_bases)
        elif hasattr(obj, '__bases__') and obj.__bases__:
            # A normal class
            obj_bases = list(obj.__bases__)
        else:
            obj_bases = []
        full_name = '.'.join((module_name, *parts))
        events.emit(
            'autodoc-process-bases',
            full_name,
            obj,
            SimpleNamespace(),
            obj_bases,
        )
        base_classes = tuple(restify(cls, mode=render_mode) for cls in obj_bases)

        return _ClassDefProperties(
            obj_type=objtype,  # type: ignore[arg-type]
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            bases=getattr(obj, '__bases__', None),
            _obj=obj,
            _obj___module__=get_attr(obj, '__module__', None),
            _obj___name__=getattr(obj, '__name__', None),
            _obj___qualname__=getattr(obj, '__qualname__', None),
            _obj_bases=base_classes,
            _obj_is_new_type=isinstance(obj, NewType),
            _obj_is_typevar=isinstance(obj, TypeVar),
        )

    if objtype in {'function', 'decorator'}:
        if inspect.isstaticmethod(obj, cls=parent, name=object_name):
            obj_properties.add('staticmethod')
        if inspect.isclassmethod(obj):
            obj_properties.add('classmethod')
        if inspect.iscoroutinefunction(obj) or inspect.isasyncgenfunction(obj):
            obj_properties.add('async')

        return _FunctionDefProperties(
            obj_type=objtype,  # type: ignore[arg-type]
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            properties=frozenset(obj_properties),
            _obj=obj,
            _obj___module__=get_attr(obj, '__module__', None),
            _obj___name__=getattr(obj, '__name__', None),
            _obj___qualname__=getattr(obj, '__qualname__', None),
        )

    if objtype == 'method':
        # to distinguish classmethod/staticmethod
        obj_ = parent.__dict__.get(object_name, obj)
        if inspect.isstaticmethod(obj_, cls=parent, name=object_name):
            obj_properties.add('staticmethod')
        elif (
            inspect.is_classmethod_like(obj_)
            or inspect.is_singledispatch_method(obj_)
            and inspect.is_classmethod_like(obj_.func)
        ):
            obj_properties.add('classmethod')
        if inspect.isabstractmethod(obj_):
            obj_properties.add('abstractmethod')
        if inspect.iscoroutinefunction(obj_) or inspect.isasyncgenfunction(obj_):
            obj_properties.add('async')

        return _FunctionDefProperties(
            obj_type=objtype,
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            properties=frozenset(obj_properties),
            _obj=obj,
            _obj___module__=get_attr(obj, '__module__', None),
            _obj___name__=getattr(obj, '__name__', None),
            _obj___qualname__=getattr(obj, '__qualname__', None),
        )

    if objtype == 'property':
        if not inspect.isproperty(obj):
            # Support for class properties. Note: these only work on Python 3.9.
            __dict__ = safe_getattr(parent, '__dict__', {})
            obj = __dict__.get(parts[-1])
            if isinstance(obj, classmethod) and inspect.isproperty(obj.__func__):
                obj = obj.__func__
                obj_properties.add('classmethod')
            else:
                return None
        if inspect.isabstractmethod(obj):
            obj_properties.add('abstractmethod')

        # get property return type annotation
        obj_property_type_annotation = None
        if safe_getattr(obj, 'fget', None):  # property
            func = obj.fget  # type: ignore[union-attr]
        elif safe_getattr(obj, 'func', None):  # cached_property
            func = obj.func  # type: ignore[union-attr]
        else:
            func = None
        if func is not None:
            # update the annotations of the property getter
            if config.autodoc_use_type_comments:
                _update_annotations_using_type_comments(func, False)

            try:
                signature = inspect.signature(
                    func, type_aliases=config.autodoc_type_aliases
                )
            except TypeError as exc:
                full_name = '.'.join((module_name, *parts))
                LOGGER.warning(
                    __('Failed to get a function signature for %s: %s'),
                    full_name,
                    exc,
                )
                pass
            except ValueError:
                pass
            else:
                if signature.return_annotation is not Parameter.empty:
                    short_literals = config.python_display_short_literal_types
                    obj_property_type_annotation = stringify_annotation(
                        signature.return_annotation,
                        render_mode,
                        short_literals=short_literals,
                    )

        return _FunctionDefProperties(
            obj_type=objtype,
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            properties=frozenset(obj_properties),
            _obj=obj,
            _obj___module__=get_attr(parent or obj, '__module__', None) or module_name,
            _obj___name__=getattr(parent or obj, '__name__', None),
            _obj___qualname__=getattr(parent or obj, '__qualname__', None),
            _obj_property_type_annotation=obj_property_type_annotation,
        )

    if objtype == 'data':
        # Update __annotations__ to support type_comment and so on
        if config.autodoc_use_type_comments:
            _ensure_annotations_from_type_comments(parent)

        # obtain annotation
        annotations = get_type_hints(
            parent,
            None,
            config.autodoc_type_aliases,
            include_extras=True,
        )
        if parts[-1] in annotations:
            short_literals = config.python_display_short_literal_types
            type_annotation = stringify_annotation(
                annotations[parts[-1]], render_mode, short_literals=short_literals
            )
        else:
            type_annotation = None

        if (
            obj is RUNTIME_INSTANCE_ATTRIBUTE
            or obj is SLOTS_ATTR
            or obj is UNINITIALIZED_ATTR
        ):
            obj_sentinel = obj
        else:
            obj_sentinel = None

        return _AssignStatementProperties(
            obj_type=objtype,
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            value=...,
            annotation='',
            class_var=False,
            instance_var=False,
            _obj=obj,
            _obj___module__=get_attr(parent or obj, '__module__', None) or module_name,
            _obj_is_generic_alias=inspect.isgenericalias(obj),
            _obj_is_attribute_descriptor=inspect.isattributedescriptor(obj),
            _obj_is_mock=ismock(obj),
            _obj_is_sentinel=obj_sentinel,
            _obj_repr_rst=inspect.object_description(obj),
            _obj_type_annotation=type_annotation,
        )

    if objtype == 'attribute':
        if _is_slots_attribute(parent=parent, obj_path=parts):
            obj = SLOTS_ATTR
        elif inspect.isenumattribute(obj):
            obj = obj.value
        if parent and config.autodoc_use_type_comments:
            # Update __annotations__ to support type_comment and so on
            _ensure_annotations_from_type_comments(parent)

        # obtain annotation
        annotations = get_type_hints(
            parent,
            None,
            config.autodoc_type_aliases,
            include_extras=True,
        )
        if parts[-1] in annotations:
            short_literals = config.python_display_short_literal_types
            type_annotation = stringify_annotation(
                annotations[parts[-1]], render_mode, short_literals=short_literals
            )
        else:
            type_annotation = None

        if (
            obj is RUNTIME_INSTANCE_ATTRIBUTE
            or obj is SLOTS_ATTR
            or obj is UNINITIALIZED_ATTR
        ):
            obj_sentinel = obj
        else:
            obj_sentinel = None

        return _AssignStatementProperties(
            obj_type=objtype,
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            value=...,
            annotation='',
            class_var=False,
            instance_var=False,
            _obj=obj,
            _obj___module__=get_attr(obj, '__module__', None),
            _obj_is_generic_alias=inspect.isgenericalias(obj),
            _obj_is_attribute_descriptor=inspect.isattributedescriptor(obj),
            _obj_is_mock=ismock(obj),
            _obj_is_sentinel=obj_sentinel,
            _obj_repr_rst=inspect.object_description(obj),
            _obj_type_annotation=type_annotation,
        )

    if objtype == 'type':
        obj_module_name = getattr(obj, '__module__', module_name)
        if obj_module_name != module_name and module_name.startswith(obj_module_name):
            bases = module_name[len(obj_module_name) :].strip('.').split('.')
            parts = tuple(bases) + parts
            module_name = obj_module_name

        short_literals = config.python_display_short_literal_types
        ann = stringify_annotation(
            obj.__value__, render_mode, short_literals=short_literals
        )
        return _TypeStatementProperties(
            obj_type=objtype,
            module_name=module_name,
            parts=parts,
            docstring_lines=(),
            _obj=obj,
            _obj___module__=get_attr(obj, '__module__', None),
            _obj___name__=getattr(obj, '__name__', None),
            _obj___qualname__=getattr(obj, '__qualname__', None),
            _obj___value__=ann,
        )

    return _ItemProperties(
        obj_type=objtype,
        module_name=module_name,
        parts=parts,
        docstring_lines=(),
        _obj=obj,
        _obj___module__=get_attr(obj, '__module__', None),
    )


def _is_slots_attribute(*, parent: Any, obj_path: Sequence[str]) -> bool:
    """Check the subject is an attribute in __slots__."""
    try:
        if parent___slots__ := inspect.getslots(parent):
            return obj_path[-1] in parent___slots__
        else:
            return False
    except (ValueError, TypeError):
        return False
