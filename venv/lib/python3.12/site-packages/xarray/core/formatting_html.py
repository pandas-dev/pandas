from __future__ import annotations

import uuid
from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass
from functools import lru_cache, partial
from html import escape
from importlib.resources import files
from math import ceil
from typing import TYPE_CHECKING

from xarray.core.formatting import (
    _coord_sort_key,
    filter_nondefault_indexes,
    inherited_vars,
    inline_index_repr,
    inline_variable_array_repr,
    short_data_repr,
)
from xarray.core.options import OPTIONS, _get_boolean_with_default

STATIC_FILES = (
    ("xarray.static.html", "icons-svg-inline.html"),
    ("xarray.static.css", "style.css"),
)

if TYPE_CHECKING:
    from xarray.core.datatree import DataTree


@lru_cache(None)
def _load_static_files():
    """Lazily load the resource files into memory the first time they are needed"""
    return [
        files(package).joinpath(resource).read_text(encoding="utf-8")
        for package, resource in STATIC_FILES
    ]


def short_data_repr_html(array) -> str:
    """Format "data" for DataArray and Variable."""
    internal_data = getattr(array, "variable", array)._data
    if hasattr(internal_data, "_repr_html_"):
        return internal_data._repr_html_()
    text = escape(short_data_repr(array))
    return f"<pre>{text}</pre>"


def format_dims(dim_sizes, dims_with_index) -> str:
    if not dim_sizes:
        return ""

    dim_css_map = {
        dim: " class='xr-has-index'" if dim in dims_with_index else ""
        for dim in dim_sizes
    }

    dims_li = "".join(
        f"<li><span{dim_css_map[dim]}>{escape(str(dim))}</span>: {size}</li>"
        for dim, size in dim_sizes.items()
    )

    return f"<ul class='xr-dim-list'>{dims_li}</ul>"


def summarize_attrs(attrs) -> str:
    attrs_dl = "".join(
        f"<dt><span>{escape(str(k))} :</span></dt><dd>{escape(str(v))}</dd>"
        for k, v in attrs.items()
    )

    return f"<dl class='xr-attrs'>{attrs_dl}</dl>"


def _icon(icon_name) -> str:
    # icon_name should be defined in xarray/static/html/icon-svg-inline.html
    return (
        f"<svg class='icon xr-{icon_name}'><use xlink:href='#{icon_name}'></use></svg>"
    )


def summarize_variable(name, var, is_index=False, dtype=None) -> str:
    variable = var.variable if hasattr(var, "variable") else var

    cssclass_idx = " class='xr-has-index'" if is_index else ""
    dims_str = f"({', '.join(escape(dim) for dim in var.dims)})"
    name = escape(str(name))
    dtype = dtype or escape(str(var.dtype))

    # "unique" ids required to expand/collapse subsections
    attrs_id = "attrs-" + str(uuid.uuid4())
    data_id = "data-" + str(uuid.uuid4())
    disabled = "" if len(var.attrs) else "disabled"

    preview = escape(inline_variable_array_repr(variable, 35))
    attrs_ul = summarize_attrs(var.attrs)
    data_repr = short_data_repr_html(variable)

    attrs_icon = _icon("icon-file-text2")
    data_icon = _icon("icon-database")

    return (
        f"<div class='xr-var-name'><span{cssclass_idx}>{name}</span></div>"
        f"<div class='xr-var-dims'>{dims_str}</div>"
        f"<div class='xr-var-dtype'>{dtype}</div>"
        f"<div class='xr-var-preview xr-preview'>{preview}</div>"
        f"<input id='{attrs_id}' class='xr-var-attrs-in' "
        f"type='checkbox' {disabled}>"
        f"<label for='{attrs_id}' title='Show/Hide attributes'>"
        f"{attrs_icon}</label>"
        f"<input id='{data_id}' class='xr-var-data-in' type='checkbox'>"
        f"<label for='{data_id}' title='Show/Hide data repr'>"
        f"{data_icon}</label>"
        f"<div class='xr-var-attrs'>{attrs_ul}</div>"
        f"<div class='xr-var-data'>{data_repr}</div>"
    )


def summarize_coords(variables) -> str:
    li_items = []
    dims = tuple(variables._data.dims)
    dim_ordered_coords = sorted(
        variables.items(), key=partial(_coord_sort_key, dims=dims)
    )
    for k, v in dim_ordered_coords:
        li_content = summarize_variable(k, v, is_index=k in variables.xindexes)
        li_items.append(f"<li class='xr-var-item'>{li_content}</li>")

    vars_li = "".join(li_items)

    return f"<ul class='xr-var-list'>{vars_li}</ul>"


def summarize_vars(variables) -> str:
    vars_li = "".join(
        f"<li class='xr-var-item'>{summarize_variable(k, v)}</li>"
        for k, v in variables.items()
    )

    return f"<ul class='xr-var-list'>{vars_li}</ul>"


def short_index_repr_html(index) -> str:
    if hasattr(index, "_repr_html_"):
        return index._repr_html_()

    return f"<pre>{escape(repr(index))}</pre>"


def summarize_index(coord_names, index) -> str:
    name = "<br>".join([escape(str(n)) for n in coord_names])

    index_id = f"index-{uuid.uuid4()}"
    preview = escape(inline_index_repr(index, max_width=70))
    details = short_index_repr_html(index)

    data_icon = _icon("icon-database")

    return (
        f"<div class='xr-index-name'><div>{name}</div></div>"
        f"<div class='xr-index-preview'>{preview}</div>"
        # need empty input + label here to conform to the fixed CSS grid layout
        f"<input type='checkbox' disabled/>"
        f"<label></label>"
        f"<input id='{index_id}' class='xr-index-data-in' type='checkbox'/>"
        f"<label for='{index_id}' title='Show/Hide index repr'>{data_icon}</label>"
        f"<div class='xr-index-data'>{details}</div>"
    )


def summarize_indexes(indexes) -> str:
    indexes_li = "".join(
        f"<li class='xr-var-item'>{summarize_index(v, i)}</li>"
        for v, i in indexes.items()
    )
    return f"<ul class='xr-var-list'>{indexes_li}</ul>"


def collapsible_section(
    header: str,
    inline_details="",
    details="",
    n_items=None,
    enabled=True,
    collapsed=False,
    span_grid=False,
) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())

    has_items = n_items is not None and n_items
    n_items_span = "" if n_items is None else f" <span>({n_items})</span>"
    enabled_attr = "" if enabled and has_items else " disabled"
    collapsed_attr = "" if collapsed or not has_items else " checked"
    span_grid_attr = " xr-span-grid" if span_grid else ""
    tip = " title='Expand/collapse section'" if enabled_attr == "" else ""
    html = (
        f"<input id='{data_id}' class='xr-section-summary-in' type='checkbox'{enabled_attr}{collapsed_attr} />"
        f"<label for='{data_id}' class='xr-section-summary{span_grid_attr}'{tip}>{header}{n_items_span}</label>"
        f"<div class='xr-section-inline-details'>{inline_details}</div>"
    )
    if details:
        html += f"<div class='xr-section-details'>{details}</div>"
    return html


def _mapping_section(
    mapping,
    name,
    details_func,
    max_items_collapse,
    expand_option_name,
    enabled=True,
    **kwargs,
) -> str:
    n_items = len(mapping)
    expanded = max_items_collapse is None or _get_boolean_with_default(
        expand_option_name, n_items < max_items_collapse
    )
    collapsed = not expanded

    return collapsible_section(
        f"{name}:",
        details=details_func(mapping, **kwargs),
        n_items=n_items,
        enabled=enabled,
        collapsed=collapsed,
    )


def dim_section(obj) -> str:
    dim_list = format_dims(obj.sizes, obj.xindexes.dims)

    return collapsible_section(
        "Dimensions:", inline_details=dim_list, enabled=False, collapsed=True
    )


def array_section(obj) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())
    collapsed = (
        "checked"
        if _get_boolean_with_default("display_expand_data", default=True)
        else ""
    )
    variable = getattr(obj, "variable", obj)
    preview = escape(inline_variable_array_repr(variable, max_width=70))
    data_repr = short_data_repr_html(obj)
    data_icon = _icon("icon-database")

    return (
        "<div class='xr-array-wrap'>"
        f"<input id='{data_id}' class='xr-array-in' type='checkbox' {collapsed}>"
        f"<label for='{data_id}' title='Show/hide data repr'>{data_icon}</label>"
        f"<div class='xr-array-preview xr-preview'><span>{preview}</span></div>"
        f"<div class='xr-array-data'>{data_repr}</div>"
        "</div>"
    )


coord_section = partial(
    _mapping_section,
    name="Coordinates",
    details_func=summarize_coords,
    max_items_collapse=25,
    expand_option_name="display_expand_coords",
)

datavar_section = partial(
    _mapping_section,
    name="Data variables",
    details_func=summarize_vars,
    max_items_collapse=15,
    expand_option_name="display_expand_data_vars",
)

index_section = partial(
    _mapping_section,
    name="Indexes",
    details_func=summarize_indexes,
    max_items_collapse=0,
    expand_option_name="display_expand_indexes",
)

attr_section = partial(
    _mapping_section,
    name="Attributes",
    details_func=summarize_attrs,
    max_items_collapse=10,
    expand_option_name="display_expand_attrs",
)


def _get_indexes_dict(indexes):
    return {
        tuple(index_vars.keys()): idx for idx, index_vars in indexes.group_by_index()
    }


def _sections_repr(sections: list[str]) -> str:
    section_items = "".join(f"<li class='xr-section-item'>{s}</li>" for s in sections)
    return f"<ul class='xr-sections'>{section_items}</ul>"


def _obj_repr(obj, header_components, sections):
    """Return HTML repr of an xarray object.

    If CSS is not injected (untrusted notebook), fallback to the plain text repr.

    """
    header = f"<div class='xr-header'>{''.join(h for h in header_components)}</div>"

    icons_svg, css_style = _load_static_files()
    return (
        "<div>"
        f"{icons_svg}<style>{css_style}</style>"
        f"<pre class='xr-text-repr-fallback'>{escape(repr(obj))}</pre>"
        "<div class='xr-wrap' style='display:none'>"
        f"{header}"
        f"{_sections_repr(sections)}"
        "</div>"
        "</div>"
    )


def array_repr(arr) -> str:
    dims = OrderedDict((k, v) for k, v in zip(arr.dims, arr.shape, strict=True))
    if hasattr(arr, "xindexes"):
        indexed_dims = arr.xindexes.dims
    else:
        indexed_dims = {}

    obj_type = f"xarray.{type(arr).__name__}"
    arr_name = escape(repr(arr.name)) if getattr(arr, "name", None) else ""

    header_components = [
        f"<div class='xr-obj-type'>{obj_type}</div>",
        f"<div class='xr-obj-name'>{arr_name}</div>",
        format_dims(dims, indexed_dims),
    ]

    sections = [array_section(arr)]

    if hasattr(arr, "coords"):
        if arr.coords:
            sections.append(coord_section(arr.coords))

    if hasattr(arr, "xindexes"):
        display_default_indexes = _get_boolean_with_default(
            "display_default_indexes", False
        )
        xindexes = filter_nondefault_indexes(
            _get_indexes_dict(arr.xindexes), not display_default_indexes
        )
        if xindexes:
            indexes = _get_indexes_dict(arr.xindexes)
            sections.append(index_section(indexes))

    if arr.attrs:
        sections.append(attr_section(arr.attrs))

    return _obj_repr(arr, header_components, sections)


def dataset_repr(ds) -> str:
    obj_type = f"xarray.{type(ds).__name__}"

    header_components = [f"<div class='xr-obj-type'>{escape(obj_type)}</div>"]

    sections = []

    sections.append(dim_section(ds))

    if ds.coords:
        sections.append(coord_section(ds.coords))

    sections.append(datavar_section(ds.data_vars))

    display_default_indexes = _get_boolean_with_default(
        "display_default_indexes", False
    )
    xindexes = filter_nondefault_indexes(
        _get_indexes_dict(ds.xindexes), not display_default_indexes
    )
    if xindexes:
        sections.append(index_section(xindexes))

    if ds.attrs:
        sections.append(attr_section(ds.attrs))

    return _obj_repr(ds, header_components, sections)


inherited_coord_section = partial(
    _mapping_section,
    name="Inherited coordinates",
    details_func=summarize_coords,
    max_items_collapse=25,
    expand_option_name="display_expand_coords",
)


def _datatree_node_sections(node: DataTree, root: bool) -> tuple[list[str], int]:
    from xarray.core.coordinates import Coordinates

    ds = node._to_dataset_view(rebuild_dims=False, inherit=True)
    node_coords = node.to_dataset(inherit=False).coords

    # use this class to get access to .xindexes property
    inherited_coords = Coordinates(
        coords=inherited_vars(node._coord_variables),
        indexes=inherited_vars(node._indexes),
    )

    # Only show dimensions if also showing a variable or coordinates section.
    show_dims = node_coords or (root and inherited_coords) or ds.data_vars

    display_default_indexes = _get_boolean_with_default(
        "display_default_indexes", False
    )
    xindexes = filter_nondefault_indexes(
        _get_indexes_dict(ds.xindexes), not display_default_indexes
    )

    sections = []
    if show_dims:
        sections.append(dim_section(ds))
    if node_coords:
        sections.append(coord_section(node_coords))
    if root and inherited_coords:
        sections.append(inherited_coord_section(inherited_coords))
    if ds.data_vars:
        sections.append(datavar_section(ds.data_vars))
    if xindexes:
        sections.append(index_section(xindexes))
    if ds.attrs:
        sections.append(attr_section(ds.attrs))

    displayed_line_count = (
        len(node.children)
        + int(bool(show_dims))
        + int(bool(node_coords))
        + len(node_coords)
        + int(root) * (int(bool(inherited_coords)) + len(inherited_coords))
        + int(bool(ds.data_vars))
        + len(ds.data_vars)
        + int(bool(xindexes))
        + len(xindexes)
        + int(bool(ds.attrs))
        + len(ds.attrs)
    )

    return sections, displayed_line_count


def _tree_item_count(node: DataTree, cache: dict[int, int]) -> int:
    if id(node) in cache:
        return cache[id(node)]

    node_ds = node.to_dataset(inherit=False)
    node_count = len(node_ds.variables) + len(node_ds.attrs)
    child_count = sum(
        _tree_item_count(child, cache) for child in node.children.values()
    )
    total = node_count + child_count
    cache[id(node)] = total
    return total


@dataclass
class _DataTreeDisplay:
    node: DataTree
    sections: list[str]
    item_count: int
    collapsed: bool
    disabled: bool


def _build_datatree_displays(tree: DataTree) -> dict[str, _DataTreeDisplay]:
    displayed_line_count = 0
    html_line_count = 0
    displays: dict[str, _DataTreeDisplay] = {}
    item_count_cache: dict[int, int] = {}
    root = True
    collapsed = False
    disabled = False

    html_limit = OPTIONS["display_max_html_elements"]
    uncollapsed_limit = OPTIONS["display_max_items"]

    too_many_items_section = collapsible_section(
        "<em>Too many items to display (display_max_html_elements exceeded)</em>",
        enabled=False,
        collapsed=True,
        span_grid=True,
    )

    for node in tree.subtree:  # breadth-first
        parent = node.parent
        if parent is not None:
            parent_display = displays.get(parent.path, None)
            if parent_display is not None and parent_display.disabled:
                break  # no need to build display

        item_count = _tree_item_count(node, item_count_cache)

        sections, node_line_count = _datatree_node_sections(node, root)
        new_displayed_count = displayed_line_count + node_line_count
        new_html_count = html_line_count + node_line_count

        disabled = not root and (disabled or new_html_count > html_limit)
        if disabled:
            sections = [too_many_items_section]
            collapsed = True
        else:
            html_line_count = new_html_count

        collapsed = not root and (collapsed or new_displayed_count > uncollapsed_limit)
        if not collapsed:
            displayed_line_count = new_displayed_count

        displays[node.path] = _DataTreeDisplay(
            node, sections, item_count, collapsed, disabled
        )
        root = False

    # If any node is collapsed, ensure its immediate siblings are also collapsed
    for display in displays.values():
        if not display.disabled:
            if any(
                displays[child.path].collapsed
                for child in display.node.children.values()
            ):
                for child in display.node.children.values():
                    displays[child.path].collapsed = True

    return displays


def _ellipsis_element() -> str:
    """Create an ellipsis element for truncated children."""
    return (
        "<div class='xr-group-box'>"
        "<div class='xr-group-box-vline' style='height: 100%'></div>"
        "<div class='xr-group-box-contents'>"
        "<div class='xr-group-box-ellipsis'>⋮</div>"
        "</div>"
        "</div>"
    )


def children_section(
    children: Mapping[str, DataTree], displays: dict[str, _DataTreeDisplay]
) -> str:
    child_elements = []
    children_list = list(children.values())
    nchildren = len(children_list)
    max_children = int(OPTIONS["display_max_children"])

    if nchildren <= max_children:
        # Render all children
        for i, child in enumerate(children_list):
            is_last = i == nchildren - 1
            child_elements.append(datatree_child_repr(child, displays, end=is_last))
    else:
        # Truncate: show first ceil(max/2), ellipsis, last floor(max/2)
        first_n = ceil(max_children / 2)
        last_n = max_children - first_n

        child_elements.extend(
            datatree_child_repr(children_list[i], displays, end=False)
            for i in range(first_n)
        )

        child_elements.append(_ellipsis_element())

        child_elements.extend(
            datatree_child_repr(children_list[i], displays, end=(i == nchildren - 1))
            for i in range(nchildren - last_n, nchildren)
        )

    children_html = "".join(child_elements)
    return f"<div class='xr-children'>{children_html}</div>"


def datatree_sections(
    node: DataTree, displays: dict[str, _DataTreeDisplay]
) -> list[str]:
    display = displays[node.path]
    sections = []
    if node.children and not display.disabled:
        sections.append(children_section(node.children, displays))
    sections.extend(display.sections)
    return sections


def datatree_child_repr(
    node: DataTree,
    displays: dict[str, _DataTreeDisplay],
    end: bool,
) -> str:
    # Wrap DataTree HTML representation with a tee to the left of it.
    #
    # Enclosing HTML tag is a <div> with :code:`display: inline-grid` style.
    #
    # Turns:
    # [    title    ]
    # |   details   |
    # |_____________|
    #
    # into (A):
    # |─ [    title    ]
    # |  |   details   |
    # |  |_____________|
    #
    # or (B):
    # └─ [    title    ]
    #    |   details   |
    #    |_____________|

    vline_height = "1.2em" if end else "100%"

    path = escape(node.path)
    display = displays[node.path]

    group_id = "group-" + str(uuid.uuid4())
    collapsed = " checked" if display.collapsed else ""
    tip = " title='Expand/collapse group'" if not display.disabled else ""

    sections = datatree_sections(node, displays)
    sections_html = _sections_repr(sections) if sections else ""

    html = f"""
        <div class='xr-group-box'>
            <div class='xr-group-box-vline' style='height: {vline_height}'></div>
            <div class='xr-group-box-hline'></div>
            <div class='xr-group-box-contents'>
                <input id='{group_id}' type='checkbox'{collapsed} />
                <label for='{group_id}'{tip}>
                    {path}
                    <span>({display.item_count})</span>
                </label>
                {sections_html}
            </div>
        </div>
    """
    return "".join(t.strip() for t in html.split("\n"))


def datatree_repr(node: DataTree) -> str:
    displays = _build_datatree_displays(node)
    header_components = [
        f"<div class='xr-obj-type'>xarray.{type(node).__name__}</div>",
    ]
    if node.name is not None:
        name = escape(repr(node.name))
        header_components.append(f"<div class='xr-obj-name'>{name}</div>")
    sections = datatree_sections(node, displays)
    return _obj_repr(node, header_components, sections)
