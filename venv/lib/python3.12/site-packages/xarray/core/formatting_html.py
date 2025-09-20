from __future__ import annotations

import uuid
from collections import OrderedDict
from collections.abc import Mapping
from functools import lru_cache, partial
from html import escape
from importlib.resources import files
from math import ceil
from typing import TYPE_CHECKING, Literal

from xarray.core.formatting import (
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
    for k, v in variables.items():
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
    name, inline_details="", details="", n_items=None, enabled=True, collapsed=False
) -> str:
    # "unique" id to expand/collapse the section
    data_id = "section-" + str(uuid.uuid4())

    has_items = n_items is not None and n_items
    n_items_span = "" if n_items is None else f" <span>({n_items})</span>"
    enabled = "" if enabled and has_items else "disabled"
    collapsed = "" if collapsed or not has_items else "checked"
    tip = " title='Expand/collapse section'" if enabled else ""

    return (
        f"<input id='{data_id}' class='xr-section-summary-in' "
        f"type='checkbox' {enabled} {collapsed}>"
        f"<label for='{data_id}' class='xr-section-summary' {tip}>"
        f"{name}:{n_items_span}</label>"
        f"<div class='xr-section-inline-details'>{inline_details}</div>"
        f"<div class='xr-section-details'>{details}</div>"
    )


def _mapping_section(
    mapping,
    name,
    details_func,
    max_items_collapse,
    expand_option_name,
    enabled=True,
    max_option_name: Literal["display_max_children"] | None = None,
) -> str:
    n_items = len(mapping)
    expanded = _get_boolean_with_default(
        expand_option_name, n_items < max_items_collapse
    )
    collapsed = not expanded

    inline_details = ""
    if max_option_name and max_option_name in OPTIONS:
        max_items = int(OPTIONS[max_option_name])
        if n_items > max_items:
            inline_details = f"({max_items}/{n_items})"

    return collapsible_section(
        name,
        inline_details=inline_details,
        details=details_func(mapping),
        n_items=n_items,
        enabled=enabled,
        collapsed=collapsed,
    )


def dim_section(obj) -> str:
    dim_list = format_dims(obj.sizes, obj.xindexes.dims)

    return collapsible_section(
        "Dimensions", inline_details=dim_list, enabled=False, collapsed=True
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


def _obj_repr(obj, header_components, sections):
    """Return HTML repr of an xarray object.

    If CSS is not injected (untrusted notebook), fallback to the plain text repr.

    """
    header = f"<div class='xr-header'>{''.join(h for h in header_components)}</div>"
    sections = "".join(f"<li class='xr-section-item'>{s}</li>" for s in sections)

    icons_svg, css_style = _load_static_files()
    return (
        "<div>"
        f"{icons_svg}<style>{css_style}</style>"
        f"<pre class='xr-text-repr-fallback'>{escape(repr(obj))}</pre>"
        "<div class='xr-wrap' style='display:none'>"
        f"{header}"
        f"<ul class='xr-sections'>{sections}</ul>"
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
    arr_name = f"'{arr.name}'" if getattr(arr, "name", None) else ""

    header_components = [
        f"<div class='xr-obj-type'>{obj_type}</div>",
        f"<div class='xr-array-name'>{arr_name}</div>",
        format_dims(dims, indexed_dims),
    ]

    sections = [array_section(arr)]

    if hasattr(arr, "coords"):
        sections.append(coord_section(arr.coords))

    if hasattr(arr, "xindexes"):
        indexes = _get_indexes_dict(arr.xindexes)
        sections.append(index_section(indexes))

    sections.append(attr_section(arr.attrs))

    return _obj_repr(arr, header_components, sections)


def dataset_repr(ds) -> str:
    obj_type = f"xarray.{type(ds).__name__}"

    header_components = [f"<div class='xr-obj-type'>{escape(obj_type)}</div>"]

    sections = [
        dim_section(ds),
        coord_section(ds.coords),
        datavar_section(ds.data_vars),
        index_section(_get_indexes_dict(ds.xindexes)),
        attr_section(ds.attrs),
    ]

    return _obj_repr(ds, header_components, sections)


def summarize_datatree_children(children: Mapping[str, DataTree]) -> str:
    MAX_CHILDREN = OPTIONS["display_max_children"]
    n_children = len(children)

    children_html = []
    for i, (n, c) in enumerate(children.items()):
        if i < ceil(MAX_CHILDREN / 2) or i >= ceil(n_children - MAX_CHILDREN / 2):
            is_last = i == (n_children - 1)
            children_html.append(
                _wrap_datatree_repr(datatree_node_repr(n, c), end=is_last)
            )
        elif n_children > MAX_CHILDREN and i == ceil(MAX_CHILDREN / 2):
            children_html.append("<div>...</div>")

    return "".join(
        [
            "<div style='display: inline-grid; grid-template-columns: 100%; grid-column: 1 / -1'>",
            "".join(children_html),
            "</div>",
        ]
    )


children_section = partial(
    _mapping_section,
    name="Groups",
    details_func=summarize_datatree_children,
    max_items_collapse=1,
    max_option_name="display_max_children",
    expand_option_name="display_expand_groups",
)

inherited_coord_section = partial(
    _mapping_section,
    name="Inherited coordinates",
    details_func=summarize_coords,
    max_items_collapse=25,
    expand_option_name="display_expand_coords",
)


def datatree_node_repr(group_title: str, node: DataTree, show_inherited=False) -> str:
    from xarray.core.coordinates import Coordinates

    header_components = [f"<div class='xr-obj-type'>{escape(group_title)}</div>"]

    ds = node._to_dataset_view(rebuild_dims=False, inherit=True)
    node_coords = node.to_dataset(inherit=False).coords

    # use this class to get access to .xindexes property
    inherited_coords = Coordinates(
        coords=inherited_vars(node._coord_variables),
        indexes=inherited_vars(node._indexes),
    )

    sections = [
        children_section(node.children),
        dim_section(ds),
        coord_section(node_coords),
    ]

    # only show inherited coordinates on the root
    if show_inherited:
        sections.append(inherited_coord_section(inherited_coords))

    sections += [
        datavar_section(ds.data_vars),
        attr_section(ds.attrs),
    ]

    return _obj_repr(ds, header_components, sections)


def _wrap_datatree_repr(r: str, end: bool = False) -> str:
    """
    Wrap HTML representation with a tee to the left of it.

    Enclosing HTML tag is a <div> with :code:`display: inline-grid` style.

    Turns:
    [    title    ]
    |   details   |
    |_____________|

    into (A):
    |─ [    title    ]
    |  |   details   |
    |  |_____________|

    or (B):
    └─ [    title    ]
       |   details   |
       |_____________|

    Parameters
    ----------
    r: str
        HTML representation to wrap.
    end: bool
        Specify if the line on the left should continue or end.

        Default is True.

    Returns
    -------
    str
        Wrapped HTML representation.

        Tee color is set to the variable :code:`--xr-border-color`.
    """
    # height of line
    end = bool(end)
    height = "100%" if end is False else "1.2em"
    return "".join(
        [
            "<div style='display: inline-grid; grid-template-columns: 0px 20px auto; width: 100%;'>",
            "<div style='",
            "grid-column-start: 1;",
            "border-right: 0.2em solid;",
            "border-color: var(--xr-border-color);",
            f"height: {height};",
            "width: 0px;",
            "'>",
            "</div>",
            "<div style='",
            "grid-column-start: 2;",
            "grid-row-start: 1;",
            "height: 1em;",
            "width: 20px;",
            "border-bottom: 0.2em solid;",
            "border-color: var(--xr-border-color);",
            "'>",
            "</div>",
            "<div style='",
            "grid-column-start: 3;",
            "'>",
            r,
            "</div>",
            "</div>",
        ]
    )


def datatree_repr(dt: DataTree) -> str:
    obj_type = f"xarray.{type(dt).__name__}"
    return datatree_node_repr(obj_type, dt, show_inherited=True)
