"""Generating content for autodoc using typehints"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, cast

from docutils import nodes

from sphinx import addnodes
from sphinx.ext.autodoc._dynamic._type_annotations import _record_typehints

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any

    from docutils.nodes import Element

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc._legacy_class_based._directive_options import Options
    from sphinx.ext.autodoc._property_types import _AutodocObjType


# Retained: legacy class-based
def record_typehints(
    app: Sphinx,
    objtype: str,
    name: str,
    obj: Any,
    options: Options,
    args: str,
    retann: str,
) -> None:
    """Record type hints to env object."""
    _record_typehints(
        autodoc_annotations=app.env.current_document.autodoc_annotations,
        name=name,
        obj=obj,
        short_literals=app.config.python_display_short_literal_types,
        type_aliases=app.config.autodoc_type_aliases,
        unqualified_typehints=app.config.autodoc_typehints_format == 'short',
    )


def _merge_typehints(
    app: Sphinx, domain: str, obj_type: _AutodocObjType, contentnode: Element
) -> None:
    if domain != 'py':
        return
    if app.config.autodoc_typehints not in {'both', 'description'}:
        return

    try:
        signature = cast('addnodes.desc_signature', contentnode.parent[0])
        if signature['module']:
            fullname = f'{signature["module"]}.{signature["fullname"]}'
        else:
            fullname = signature['fullname']
    except KeyError:
        # signature node does not have valid context info for the target object
        return

    annotations = app.env.current_document.autodoc_annotations
    if annotations.get(fullname, {}):
        field_lists = [n for n in contentnode if isinstance(n, nodes.field_list)]
        if field_lists == []:
            field_list = insert_field_list(contentnode)
            field_lists.append(field_list)

        for field_list in field_lists:
            if app.config.autodoc_typehints_description_target == 'all':
                if obj_type == 'class':
                    modify_field_list(
                        field_list, annotations[fullname], suppress_rtype=True
                    )
                else:
                    modify_field_list(field_list, annotations[fullname])
            elif app.config.autodoc_typehints_description_target == 'documented_params':
                augment_descriptions_with_types(
                    field_list, annotations[fullname], force_rtype=True
                )
            else:
                augment_descriptions_with_types(
                    field_list, annotations[fullname], force_rtype=False
                )


def insert_field_list(node: Element) -> nodes.field_list:
    field_list = nodes.field_list()
    desc = [n for n in node if isinstance(n, addnodes.desc)]
    if desc:
        # insert just before sub object descriptions (ex. methods, nested classes, etc.)
        index = node.index(desc[0])
        node.insert(index - 1, [field_list])
    else:
        node += field_list

    return field_list


def modify_field_list(
    node: nodes.field_list, annotations: dict[str, str], suppress_rtype: bool = False
) -> None:
    arguments: dict[str, dict[str, bool]] = {}
    fields = cast('Iterable[nodes.field]', node)
    for field in fields:
        field_name = field[0].astext()
        parts = re.split(' +', field_name)
        if parts[0] == 'param':
            if len(parts) == 2:
                # :param xxx:
                arg = arguments.setdefault(parts[1], {})
                arg['param'] = True
            elif len(parts) > 2:
                # :param xxx yyy:
                name = ' '.join(parts[2:])
                arg = arguments.setdefault(name, {})
                arg['param'] = True
                arg['type'] = True
        elif parts[0] == 'type':
            name = ' '.join(parts[1:])
            arg = arguments.setdefault(name, {})
            arg['type'] = True
        elif parts[0] == 'rtype':
            arguments['return'] = {'type': True}

    for name, annotation in annotations.items():
        if name == 'return':
            continue

        if '*' + name in arguments:
            name = '*' + name
            arguments.get(name)
        elif '**' + name in arguments:
            name = '**' + name
            arguments.get(name)
        else:
            arg = arguments.get(name, {})

        if not arg.get('type'):
            field = nodes.field()
            field += nodes.field_name('', 'type ' + name)
            field += nodes.field_body('', nodes.paragraph('', annotation))
            node += field
        if not arg.get('param'):
            field = nodes.field()
            field += nodes.field_name('', 'param ' + name)
            field += nodes.field_body('', nodes.paragraph('', ''))
            node += field

    if 'return' in annotations and 'return' not in arguments:
        annotation = annotations['return']
        if annotation == 'None' and suppress_rtype:
            return

        field = nodes.field()
        field += nodes.field_name('', 'rtype')
        field += nodes.field_body('', nodes.paragraph('', annotation))
        node += field


def augment_descriptions_with_types(
    node: nodes.field_list,
    annotations: dict[str, str],
    force_rtype: bool,
) -> None:
    fields = cast('Iterable[nodes.field]', node)
    has_description: set[str] = set()
    has_type: set[str] = set()
    for field in fields:
        field_name = field[0].astext()
        parts = re.split(' +', field_name)
        if parts[0] == 'param':
            if len(parts) == 2:
                # :param xxx:
                has_description.add(parts[1])
            elif len(parts) > 2:
                # :param xxx yyy:
                name = ' '.join(parts[2:])
                has_description.add(name)
                has_type.add(name)
        elif parts[0] == 'type':
            name = ' '.join(parts[1:])
            has_type.add(name)
        elif parts[0] in {'return', 'returns'}:
            has_description.add('return')
        elif parts[0] == 'rtype':
            has_type.add('return')

    # Add 'type' for parameters with a description but no declared type.
    for name, annotation in annotations.items():
        if name in {'return', 'returns'}:
            continue

        if '*' + name in has_description:
            name = '*' + name
        elif '**' + name in has_description:
            name = '**' + name

        if name in has_description and name not in has_type:
            field = nodes.field()
            field += nodes.field_name('', 'type ' + name)
            field += nodes.field_body('', nodes.paragraph('', annotation))
            node += field

    # Add 'rtype' if 'return' is present and 'rtype' isn't.
    if 'return' in annotations:
        rtype = annotations['return']
        if 'return' not in has_type and (
            'return' in has_description or (force_rtype and rtype != 'None')
        ):
            field = nodes.field()
            field += nodes.field_name('', 'rtype')
            field += nodes.field_body('', nodes.paragraph('', rtype))
            node += field
