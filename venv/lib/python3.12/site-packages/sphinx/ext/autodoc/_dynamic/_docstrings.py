from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from sphinx.errors import PycodeError
from sphinx.ext.autodoc._dynamic._importer import (
    _get_attribute_comment,
    _is_runtime_instance_attribute_not_commented,
)
from sphinx.ext.autodoc._property_types import _ClassDefProperties
from sphinx.ext.autodoc._sentinels import (
    RUNTIME_INSTANCE_ATTRIBUTE,
    SLOTS_ATTR,
    UNINITIALIZED_ATTR,
)
from sphinx.ext.autodoc._shared import LOGGER
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import inspect
from sphinx.util.docstrings import prepare_docstring
from sphinx.util.inspect import getdoc

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping
    from typing import Any, Literal

    from sphinx.events import EventManager
    from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
    from sphinx.ext.autodoc._property_types import _ItemProperties
    from sphinx.ext.autodoc._shared import _AttrGetter


def _docstring_lines_for_props(
    docstrings: list[list[str]] | None,
    /,
    *,
    props: _ItemProperties,
    parent_modname: str | None,
    events: EventManager,
    options: _AutoDocumenterOptions,
) -> tuple[str, ...]:
    attr_docs = _attr_docs_for_props(props, parent_modname=parent_modname)
    prepared_docstrings = _prepare_docstrings(
        docstrings, props=props, attr_docs=attr_docs
    )
    docstring_lines = _process_docstrings(
        prepared_docstrings,
        events=events,
        props=props,
        options=options,
    )
    return tuple(docstring_lines)


def _attr_docs_for_props(
    props: _ItemProperties, *, parent_modname: str | None
) -> Mapping[tuple[str, str], list[str]]:
    if props.obj_type in {'class', 'exception'}:
        # If a class gets imported into the module ``parent_modname``
        # the analyzer won't find the source of the class,
        # if it looks in ``parent_modname``.
        real_modname = props.module_name
    elif parent_modname is None:
        real_modname = props.canonical_module_name
    else:
        real_modname = parent_modname

    try:
        analyzer = ModuleAnalyzer.for_module(real_modname)
        # parse right now, to get PycodeErrors on parsing (results will
        # be cached anyway)
        analyzer.analyze()
    except PycodeError as exc:
        LOGGER.debug('[autodoc] module analyzer failed: %s', exc)
        # no source file -- e.g. for builtin and C modules
        attr_docs = {}
    else:
        attr_docs = analyzer.attr_docs
    return attr_docs


def _prepare_docstrings(
    docstrings: list[list[str]] | None,
    *,
    props: _ItemProperties,
    attr_docs: Mapping[tuple[str, str], list[str]],
) -> list[list[str]] | None:
    """Add content from docstrings, attribute documentation and user."""
    # add content from attribute documentation
    if props.obj_type not in {'data', 'attribute', 'type'} and props.parts:
        key = ('.'.join(props.parent_names), props.name)
        try:
            # make a copy of docstring for attributes to avoid cache
            # the change of autodoc-process-docstring event.
            return [list(attr_docs[key])]
        except KeyError:
            pass

    if docstrings is None:
        return None
    if not docstrings:
        # append at least a dummy docstring, so that the event
        # autodoc-process-docstring is fired and can add some
        # content if desired
        docstrings.append([])
    return docstrings


def _process_docstrings(
    docstrings: list[list[str]] | None,
    *,
    events: EventManager,
    props: _ItemProperties,
    options: _AutoDocumenterOptions,
) -> Iterator[str]:
    """Let the user process the docstrings before adding them."""
    if docstrings is None:
        return
    for docstring_lines in docstrings:
        # let extensions pre-process docstrings
        events.emit(
            'autodoc-process-docstring',
            props.obj_type,
            props.full_name,
            props._obj,
            options,
            docstring_lines,
        )

        yield from docstring_lines
        if docstring_lines and docstring_lines[-1]:
            # ensure the docstring ends with a blank line
            yield ''


def _get_docstring_lines(
    props: _ItemProperties,
    *,
    class_doc_from: Literal['both', 'class', 'init'],
    get_attr: _AttrGetter,
    inherit_docstrings: bool,
    parent: Any,
    tab_width: int,
) -> list[list[str]] | None:
    """Decode and return lines of the docstring(s) for the object.

    When it returns None, autodoc-process-docstring will not be called for this
    object.
    """
    obj = props._obj

    if props.obj_type in {'class', 'exception'}:
        assert isinstance(props, _ClassDefProperties)

        if isinstance(obj, TypeVar):
            if obj.__doc__ == TypeVar.__doc__:
                return []
        if props.doc_as_attr:
            # Don't show the docstring of the class when it is an alias.
            if _class_variable_comment(props):
                return []
            return None

        docstrings = []
        if attr_docstring := getdoc(obj):
            docstrings.append(attr_docstring)

        # for classes, what the "docstring" is can be controlled via a
        # config value; the default is only the class docstring
        if class_doc_from in {'both', 'init'}:
            __init__ = get_attr(obj, '__init__', None)
            init_docstring = getdoc(
                __init__,
                allow_inherited=inherit_docstrings,
                cls=obj,  # TODO: object or obj?
                name='__init__',
            )
            # no __init__ means default __init__
            if init_docstring == object.__init__.__doc__:
                init_docstring = None
            if not init_docstring:
                # try __new__
                __new__ = get_attr(obj, '__new__', None)
                init_docstring = getdoc(
                    __new__,
                    allow_inherited=inherit_docstrings,
                    cls=object,  # TODO: object or obj?
                    name='__new__',
                )
                # no __new__ means default __new__
                if init_docstring == object.__new__.__doc__:
                    init_docstring = None
            if init_docstring:
                if class_doc_from == 'init':
                    docstrings = [init_docstring]
                else:
                    docstrings.append(init_docstring)

        return [prepare_docstring(docstring, tab_width) for docstring in docstrings]

    if props.obj_type == 'method':
        docstring = getdoc(
            obj,
            allow_inherited=inherit_docstrings,
            cls=parent,
            name=props.object_name,
        )
        if (
            not docstring
            or (props.name == '__init__' and docstring == object.__init__.__doc__)
            or (props.name == '__new__' and docstring == object.__new__.__doc__)
        ):
            return []
        return [prepare_docstring(docstring, tab_width)]

    if props.obj_type == 'data':
        # Check the variable has a docstring-comment

        # get_module_comment()
        comment = None
        try:
            analyzer = ModuleAnalyzer.for_module(props.module_name)
            analyzer.analyze()
            key = ('', props.name)
            if key in analyzer.attr_docs:
                comment = list(analyzer.attr_docs[key])
        except PycodeError:
            pass

        if comment:
            return [comment]

        if obj is UNINITIALIZED_ATTR:
            return []

        docstring = getdoc(
            obj,
            allow_inherited=inherit_docstrings,
            cls=parent,
            name=props.object_name,
        )
        if not docstring:
            return []
        return [prepare_docstring(docstring, tab_width)]

    if props.obj_type == 'attribute':
        # Check the attribute has a docstring-comment
        comment = _get_attribute_comment(
            parent=parent, obj_path=props.parts, attrname=props.parts[-1]
        )
        if comment:
            return [comment]

        # Disable `autodoc_inherit_docstring` to avoid to obtain
        # a docstring from the value which descriptor returns unexpectedly.
        # See: https://github.com/sphinx-doc/sphinx/issues/7805
        inherit_docstrings = False

        if obj is SLOTS_ATTR:
            # support for __slots__
            try:
                parent___slots__ = inspect.getslots(parent)
                if parent___slots__ and (docstring := parent___slots__.get(props.name)):
                    return [prepare_docstring(docstring)]
                return []
            except ValueError as exc:
                LOGGER.warning(
                    __('Invalid __slots__ found on %s. Ignored.'),
                    (parent.__qualname__, exc),
                    type='autodoc',
                )
                return []

        if (
            obj is RUNTIME_INSTANCE_ATTRIBUTE
            and _is_runtime_instance_attribute_not_commented(
                parent=parent, obj_path=props.parts
            )
        ):
            return None

        if obj is UNINITIALIZED_ATTR:
            return None

        if not inspect.isattributedescriptor(obj):
            # the docstring of non-data descriptor is very probably
            # the wrong thing to display
            return None

        docstring = getdoc(
            obj,
            allow_inherited=inherit_docstrings,
            cls=parent,
            name=props.object_name,
        )
        if not docstring:
            return []
        return [prepare_docstring(docstring, tab_width)]

    if props.obj_type == 'type':
        try:
            analyzer = ModuleAnalyzer.for_module(props.module_name)
            analyzer.analyze()
        except PycodeError:
            return None

        key = ('', props.name)
        if key in analyzer.attr_docs:
            if comment := list(analyzer.attr_docs[key]):
                return [comment]

        return None

    docstring = getdoc(
        obj,
        allow_inherited=inherit_docstrings,
        cls=parent,
        name=props.object_name,
    )
    if not docstring:
        return []
    return [prepare_docstring(docstring, tab_width)]


def _class_variable_comment(props: _ItemProperties) -> bool:
    try:
        analyzer = ModuleAnalyzer.for_module(props.module_name)
        analyzer.analyze()
        key = ('', props.dotted_parts)
        return bool(analyzer.attr_docs.get(key, False))
    except PycodeError:
        return False
