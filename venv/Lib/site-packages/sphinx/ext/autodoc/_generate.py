from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from docutils.statemachine import StringList

from sphinx.errors import PycodeError
from sphinx.ext.autodoc._dynamic._loader import _load_object_by_name
from sphinx.ext.autodoc._dynamic._member_finder import _gather_members
from sphinx.ext.autodoc._dynamic._mock import ismock
from sphinx.ext.autodoc._renderer import _add_content, _directive_header_lines
from sphinx.ext.autodoc._sentinels import ALL
from sphinx.ext.autodoc._shared import LOGGER, _get_render_mode
from sphinx.locale import _, __
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import inspect
from sphinx.util.typing import restify, stringify_annotation

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping, MutableSet

    from sphinx.environment import _CurrentDocument
    from sphinx.events import EventManager
    from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
    from sphinx.ext.autodoc._property_types import _AutodocObjType, _ItemProperties
    from sphinx.ext.autodoc._shared import _AttrGetter, _AutodocConfig
    from sphinx.util.typing import _RestifyMode


def _auto_document_object(
    *,
    config: _AutodocConfig,
    current_document: _CurrentDocument,
    events: EventManager,
    get_attr: _AttrGetter,
    more_content: StringList | None,
    name: str,
    obj_type: _AutodocObjType,
    options: _AutoDocumenterOptions,
    record_dependencies: MutableSet[str],
    ref_context: Mapping[str, str | None],
    reread_always: MutableSet[str],
) -> StringList | None:
    props = _load_object_by_name(
        name=name,
        objtype=obj_type,
        current_document=current_document,
        config=config,
        events=events,
        get_attr=get_attr,
        options=options,
        ref_context=ref_context,
        reread_always=reread_always,
    )
    if props is None:
        return None

    result = StringList()
    _generate_directives(
        more_content=more_content,
        config=config,
        current_document=current_document,
        events=events,
        get_attr=get_attr,
        indent='',
        options=options,
        props=props,
        record_dependencies=record_dependencies,
        ref_context=ref_context,
        reread_always=reread_always,
        result=result,
    )
    return result


def _generate_directives(
    more_content: StringList | None = None,
    parent_modname: str | None = None,
    check_module: bool = False,
    all_members: bool = False,
    *,
    config: _AutodocConfig,
    current_document: _CurrentDocument,
    events: EventManager,
    get_attr: _AttrGetter,
    indent: str,
    options: _AutoDocumenterOptions,
    props: _ItemProperties,
    record_dependencies: MutableSet[str],
    ref_context: Mapping[str, str | None],
    reread_always: MutableSet[str],
    result: StringList,
) -> None:
    """Generate reST for the object given by *props*, and possibly for its members.

    If *more_content* is given, include that content. If *parent_modname* is
    given, use that module name to find attribute docs. If *check_module* is
    True, only generate if the object is defined in the module name it is
    imported from. If *all_members* is True, document all members.
    """
    # If there is no parent module specified, figure out which to use.
    # The real module is used in the module analyzer to look up the module
    # where the attribute documentation would actually be found in.
    # This is used for situations where you have a module that collects the
    # functions and classes of internal submodules.
    if parent_modname is None or props.obj_type in {'class', 'exception'}:
        # If a class gets imported into the module ``parent_modname``
        # the analyzer won't find the source of the class,
        # if it looks in ``parent_modname``.
        real_modname = props.canonical_module_name
    else:
        real_modname = parent_modname

    # try to also get a source code analyzer for attribute docs
    try:
        analyzer = ModuleAnalyzer.for_module(real_modname)
        # parse right now, to get PycodeErrors on parsing (results will
        # be cached anyway)
        analyzer.analyze()
        record_dependencies.add(analyzer.srcname)
    except PycodeError as exc:
        LOGGER.debug('[autodoc] module analyzer failed: %s', exc)
        # no source file -- e.g. for builtin and C modules
        analyzer = None
        # at least add the module source file as a dependency
        if props.module_name:
            try:
                module_spec = sys.modules[props.module_name].__spec__
            except (AttributeError, KeyError):
                pass
            else:
                if (
                    module_spec is not None
                    and module_spec.has_location
                    and module_spec.origin
                ):
                    record_dependencies.add(module_spec.origin)

    if real_modname != props.canonical_module_name:
        # Add module to dependency list if target object is defined in other module.
        try:
            srcname, _ = ModuleAnalyzer.get_module_source(props.canonical_module_name)
            record_dependencies.add(str(srcname))
        except PycodeError:
            pass

    has_docstring = bool(props.docstring_lines)
    if ismock(props._obj) and not has_docstring:
        LOGGER.warning(
            __('A mocked object is detected: %r'),
            props.full_name,
            type='autodoc',
            subtype='mocked_object',
        )

    # check __module__ of object (for members not given explicitly)
    if check_module and not options.imported_members:
        subject = inspect.unpartial(props._obj)
        modname = get_attr(subject, '__module__', None)
        if modname and modname != props.module_name:
            return

    # add all content (from docstrings, attribute docs etc.)
    analyzer_source = '' if analyzer is None else analyzer.srcname
    _add_directive_lines(
        more_content=more_content,
        is_final=analyzer is not None and props.dotted_parts in analyzer.finals,
        config=config,
        indent=indent,
        options=options,
        props=props,
        result=result,
        source_name=_docstring_source_name(props=props, source=analyzer_source),
    )

    # document members, if possible
    _document_members(
        all_members=all_members,
        analyzer_order=analyzer.tagorder if analyzer is not None else {},
        attr_docs=analyzer.attr_docs if analyzer is not None else {},
        config=config,
        current_document=current_document,
        events=events,
        get_attr=get_attr,
        indent=indent,
        options=options,
        props=props,
        real_modname=real_modname,
        record_dependencies=record_dependencies,
        ref_context=ref_context,
        reread_always=reread_always,
        result=result,
    )


def _add_directive_lines(
    *,
    more_content: StringList | None,
    is_final: bool,
    config: _AutodocConfig,
    indent: str,
    options: _AutoDocumenterOptions,
    props: _ItemProperties,
    result: StringList,
    source_name: str,
) -> None:
    # generate the directive header and options, if applicable
    lines = _directive_header_lines(
        autodoc_typehints=config.autodoc_typehints,
        directive_name=(
            'py:attribute'
            if props.obj_type in {'class', 'exception'} and props.doc_as_attr  # type: ignore[attr-defined]
            else f'py:{props.obj_type}'
        ),
        is_final=is_final,
        options=options,
        props=props,
    )
    header_lines = StringList(list(lines), source='')

    # add content from docstrings or attribute documentation
    docstring_lines = StringList(props.docstring_lines, source=source_name)

    # add alias information, if applicable
    lines = _body_alias_lines(
        render_mode=_get_render_mode(config.autodoc_typehints_format),
        short_literals=config.python_display_short_literal_types,
        props=props,
    )
    alias_lines = StringList(list(lines), source='')

    # make sure that the result starts with an empty line.  This is
    # necessary for some situations where another directive preprocesses
    # reST and no starting newline is present
    result.append('', '')
    _add_content(header_lines, result=result, indent=indent)
    result.append('', '')
    _add_content(docstring_lines, result=result, indent=indent + '   ')
    if more_content is not None:
        # add additional content from the directive, if present
        _add_content(more_content, result=result, indent=indent + '   ')
    _add_content(alias_lines, result=result, indent=indent + '   ')


def _document_members(
    *,
    all_members: bool,
    analyzer_order: dict[str, int],
    attr_docs: dict[tuple[str, str], list[str]],
    config: _AutodocConfig,
    current_document: _CurrentDocument,
    events: EventManager,
    get_attr: _AttrGetter,
    indent: str,
    options: _AutoDocumenterOptions,
    props: _ItemProperties,
    real_modname: str,
    record_dependencies: MutableSet[str],
    ref_context: Mapping[str, str | None],
    reread_always: MutableSet[str],
    result: StringList,
) -> None:
    """Generate reST for member documentation.

    If *all_members* is True, document all members, else those given by
    *self.options.members*.
    """
    has_members = props.obj_type == 'module' or (
        props.obj_type in {'class', 'exception'} and not props.doc_as_attr  # type: ignore[attr-defined]
    )
    if not has_members:
        return

    want_all = bool(all_members or options.inherited_members or options.members is ALL)
    member_documenters = _gather_members(
        want_all=want_all,
        indent=indent,
        analyzer_order=analyzer_order,
        attr_docs=attr_docs,
        config=config,
        current_document=current_document,
        events=events,
        get_attr=get_attr,
        options=options,
        parent_modname=real_modname,
        props=props,
        ref_context=ref_context,
        reread_always=reread_always,
    )

    # for implicit module members, check __module__ to avoid
    # documenting imported objects
    members_check_module = bool(
        props.obj_type == 'module'
        and want_all
        and (options.ignore_module_all or props.all is None)  # type: ignore[attr-defined]
    )
    for member_props, is_attr, member_indent in member_documenters:
        assert member_props.module_name
        # Note that those two methods above do not emit events, so
        # whatever objects we deduced should not have changed.
        _generate_directives(
            more_content=None,
            parent_modname=real_modname,
            check_module=members_check_module and not is_attr,
            all_members=True,
            config=config,
            current_document=current_document,
            events=events,
            get_attr=get_attr,
            indent=member_indent,
            options=options,
            props=member_props,
            record_dependencies=record_dependencies,
            ref_context=ref_context,
            reread_always=reread_always,
            result=result,
        )


def _body_alias_lines(
    *, props: _ItemProperties, render_mode: _RestifyMode, short_literals: bool
) -> Iterator[str]:
    """Add content from docstrings, attribute documentation and user."""
    if props.obj_type in {'data', 'attribute'}:
        from sphinx.ext.autodoc._property_types import _AssignStatementProperties

        assert isinstance(props, _AssignStatementProperties)

        # Support for documenting GenericAliases
        if props._obj_is_generic_alias:
            alias = restify(props._obj, mode=render_mode)
            yield _('alias of %s') % alias
            yield ''
            return
        return

    if props.obj_type in {'class', 'exception'}:
        from sphinx.ext.autodoc._property_types import _ClassDefProperties

        assert isinstance(props, _ClassDefProperties)

        obj = props._obj

        if props._obj_is_new_type:
            supertype = restify(obj.__supertype__, mode=render_mode)
            yield _('alias of %s') % supertype
            yield ''
            return

        if props._obj_is_typevar:
            attrs = [
                repr(obj.__name__),
                *(
                    stringify_annotation(
                        constraint, render_mode, short_literals=short_literals
                    )
                    for constraint in obj.__constraints__
                ),
            ]
            if obj.__bound__:
                attrs.append(rf'bound=\ {restify(obj.__bound__, mode=render_mode)}')
            if obj.__covariant__:
                attrs.append('covariant=True')
            if obj.__contravariant__:
                attrs.append('contravariant=True')

            alias = f'TypeVar({", ".join(attrs)})'
            yield _('alias of %s') % alias
            yield ''
            return

        if props.doc_as_attr:
            try:
                analyzer = ModuleAnalyzer.for_module(props.module_name)
                analyzer.analyze()
                key = ('', props.dotted_parts)
                class_var_doc_comment = key in analyzer.attr_docs
            except PycodeError:
                class_var_doc_comment = False

            if class_var_doc_comment:
                return
            alias = restify(obj, mode=render_mode)
            yield _('alias of %s') % alias
            return

        return

    return


def _docstring_source_name(*, props: _ItemProperties, source: str) -> str:
    obj_module = inspect.safe_getattr(props._obj, '__module__', None)
    obj_qualname = inspect.safe_getattr(props._obj, '__qualname__', None)
    if obj_module and obj_qualname:
        # Get the correct location of docstring from props._obj
        # to support inherited methods
        fullname = f'{obj_module}.{obj_qualname}'
    else:
        fullname = props.full_name

    if source:
        return f'{source}:docstring of {fullname}'
    return f'docstring of {fullname}'
