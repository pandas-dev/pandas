from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.ext.autodoc._property_types import (
    _AssignStatementProperties,
    _ClassDefProperties,
    _FunctionDefProperties,
    _TypeStatementProperties,
)
from sphinx.ext.autodoc._sentinels import SUPPRESS
from sphinx.locale import _

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Literal

    from docutils.statemachine import StringList

    from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
    from sphinx.ext.autodoc._property_types import _ItemProperties


def _directive_header_lines(
    *,
    autodoc_typehints: Literal['signature', 'description', 'none', 'both'],
    directive_name: str,
    is_final: bool,
    options: _AutoDocumenterOptions,
    props: _ItemProperties,
) -> Iterator[str]:
    """Emit the directive header and option lines."""
    # normally the name doesn't contain the module
    # (except for module directives of course)
    name = props.dotted_parts or props.module_name

    # emit one signature per line
    # the first line contains the directive prefix
    sig_line, *sig_lines = props.signatures or ('',)
    prefix = f'.. {directive_name}:: '
    yield f'{prefix}{name}{sig_line}'
    # emit remaining lines, indented to the same column
    prefix = ' ' * len(prefix)
    for sig_line in sig_lines:
        yield f'{prefix}{name}{sig_line}'

    if options.no_index or options.noindex:
        yield '   :no-index:'
    if options.no_index_entry:
        yield '   :no-index-entry:'
    if props.parts:
        # Be explicit about the module, this is necessary since .. class::
        # etc. don't support a prepended module name
        yield f'   :module: {props.module_name}'

    if props.obj_type == 'module':
        # add some module-specific options
        if options.synopsis:
            yield f'   :synopsis: {options.synopsis}'
        if options.platform:
            yield f'   :platform: {options.platform}'
        if options.deprecated:
            yield '   :deprecated:'

    if props.obj_type in {'class', 'exception'}:
        assert isinstance(props, _ClassDefProperties)

        if props._obj_is_new_type or props._obj_is_typevar:
            return

        if is_final:
            yield '   :final:'

        canonical_fullname = props.canonical_full_name
        if (
            not props.doc_as_attr
            and not props._obj_is_new_type
            and canonical_fullname
            and props.full_name != canonical_fullname
        ):
            yield f'   :canonical: {canonical_fullname}'

        # add inheritance info, if wanted
        if not props.doc_as_attr and options.show_inheritance:
            yield ''
            yield '   ' + _('Bases: %s') % ', '.join(props._obj_bases)

    if props.obj_type in {'function', 'decorator'}:
        assert isinstance(props, _FunctionDefProperties)

        if props.is_async:
            yield '   :async:'

    if props.obj_type == 'method':
        assert isinstance(props, _FunctionDefProperties)

        if props.is_abstractmethod:
            yield '   :abstractmethod:'
        if props.is_async:
            yield '   :async:'
        if props.is_classmethod:
            yield '   :classmethod:'
        if props.is_staticmethod:
            yield '   :staticmethod:'
        if props.is_final or is_final:
            yield '   :final:'

    if props.obj_type == 'property':
        assert isinstance(props, _FunctionDefProperties)

        if props.is_abstractmethod:
            yield '   :abstractmethod:'
        if props.is_classmethod:
            yield '   :classmethod:'

        objrepr = props._obj_property_type_annotation
        if autodoc_typehints != 'none' and objrepr is not None:
            yield f'   :type: {objrepr}'

    if props.obj_type == 'data':
        assert isinstance(props, _AssignStatementProperties)

        if options.annotation is SUPPRESS or props._obj_is_generic_alias:
            pass
        elif options.annotation:
            yield f'   :annotation: {options.annotation}'
        else:
            type_annotation = props._obj_type_annotation
            if autodoc_typehints != 'none' and type_annotation is not None:
                yield f'   :type: {type_annotation}'

            if (
                not options.no_value
                and props._obj_is_sentinel is None  # not any sentinel
                and not props._docstrings_has_hide_value
                and not props._obj_is_mock
            ):
                yield f'   :value: {props._obj_repr_rst}'

    if props.obj_type == 'attribute':
        assert isinstance(props, _AssignStatementProperties)

        if (
            options.annotation
            and options.annotation is not SUPPRESS
            and not props._obj_is_generic_alias
        ):
            yield f'   :annotation: {options.annotation}'
        else:
            type_annotation = props._obj_type_annotation
            if autodoc_typehints != 'none' and type_annotation is not None:
                yield f'   :type: {type_annotation}'

            if (
                not options.no_value
                and props._obj_is_sentinel is None  # not any sentinel
                and not props._obj_is_attribute_descriptor
                and not props._obj_is_generic_alias
                and not props._docstrings_has_hide_value
                and not props._obj_is_mock
            ):
                yield f'   :value: {props._obj_repr_rst}'

    if props.obj_type == 'type':
        assert isinstance(props, _TypeStatementProperties)

        if not options.no_value and not props._docstrings_has_hide_value:
            yield f'   :canonical: {props._obj___value__}'


def _add_content(content: StringList, *, result: StringList, indent: str) -> None:
    for line, src in zip(content.data, content.items, strict=True):
        if line.strip():  # not a blank line
            result.append(indent + line, src[0], src[1])
        else:
            result.append('', src[0], src[1])
