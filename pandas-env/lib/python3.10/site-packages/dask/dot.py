from __future__ import annotations

import os
import re
from functools import partial

from dask.core import get_dependencies, ishashable, istask
from dask.utils import apply, funcname, import_required, key_split


def task_label(task):
    """Label for a task on a dot graph.

    Examples
    --------
    >>> from operator import add
    >>> task_label((add, 1, 2))
    'add'
    >>> task_label((add, (add, 1, 2), 3))
    'add(...)'
    """
    func = task[0]
    if func is apply:
        func = task[1]
    if hasattr(func, "funcs"):
        if len(func.funcs) > 1:
            return f"{funcname(func.funcs[0])}(...)"
        else:
            head = funcname(func.funcs[0])
    else:
        head = funcname(func)
    if any(has_sub_tasks(i) for i in task[1:]):
        return f"{head}(...)"
    else:
        return head


def has_sub_tasks(task):
    """Returns True if the task has sub tasks"""
    if istask(task):
        return True
    elif isinstance(task, list):
        return any(has_sub_tasks(i) for i in task)
    else:
        return False


def name(x):
    try:
        return str(hash(x))
    except TypeError:
        return str(hash(str(x)))


_HASHPAT = re.compile("([0-9a-z]{32})")
_UUIDPAT = re.compile("([0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12})")


def label(x, cache=None):
    """

    >>> label('x')
    'x'

    >>> label(('x', 1))
    "('x', 1)"

    >>> from hashlib import md5
    >>> x = 'x-%s-hello' % md5(b'1234').hexdigest()
    >>> x
    'x-81dc9bdb52d04dc20036dbd8313ed055-hello'

    >>> label(x)
    'x-#-hello'

    >>> from uuid import uuid1
    >>> x = 'x-%s-hello' % uuid1()
    >>> x  # doctest: +SKIP
    'x-4c1a3d7e-0b45-11e6-8334-54ee75105593-hello'

    >>> label(x)
    'x-#-hello'
    """
    s = str(x)
    for pattern in (_HASHPAT, _UUIDPAT):
        m = re.search(pattern, s)
        if m is not None:
            for h in m.groups():
                if cache is not None:
                    n = cache.get(h, len(cache))
                    label = f"#{n}"
                    # cache will be overwritten destructively
                    cache[h] = n
                else:
                    label = "#"
                s = s.replace(h, label)
    return s


def box_label(key, verbose=False):
    """Label boxes in graph by chunk index

    >>> box_label(('x', 1, 2, 3))
    '(1, 2, 3)'
    >>> box_label(('x', 123))
    '123'
    >>> box_label('x')
    ''
    """
    if isinstance(key, tuple):
        key = key[1:]
        if len(key) == 1:
            [key] = key
        return str(key)
    elif verbose:
        return str(key)
    else:
        return ""


def to_graphviz(
    dsk,
    data_attributes=None,
    function_attributes=None,
    rankdir="BT",
    graph_attr=None,
    node_attr=None,
    edge_attr=None,
    collapse_outputs=False,
    verbose=False,
    **kwargs,
):
    graphviz = import_required(
        "graphviz",
        "Drawing dask graphs with the graphviz engine requires the `graphviz` "
        "python library and the `graphviz` system library.\n\n"
        "Please either conda or pip install as follows:\n\n"
        "  conda install python-graphviz     # either conda install\n"
        "  python -m pip install graphviz    # or pip install and follow installation instructions",
    )

    data_attributes = data_attributes or {}
    function_attributes = function_attributes or {}
    graph_attr = graph_attr or {}
    node_attr = node_attr or {}
    edge_attr = edge_attr or {}

    graph_attr["rankdir"] = rankdir
    node_attr["fontname"] = "helvetica"

    graph_attr.update(kwargs)
    g = graphviz.Digraph(
        graph_attr=graph_attr, node_attr=node_attr, edge_attr=edge_attr
    )

    seen = set()
    connected = set()

    for k, v in dsk.items():
        k_name = name(k)
        if istask(v):
            func_name = name((k, "function")) if not collapse_outputs else k_name
            if collapse_outputs or func_name not in seen:
                seen.add(func_name)
                attrs = function_attributes.get(k, {}).copy()
                attrs.setdefault("label", key_split(k))
                attrs.setdefault("shape", "circle")
                g.node(func_name, **attrs)
            if not collapse_outputs:
                g.edge(func_name, k_name)
                connected.add(func_name)
                connected.add(k_name)

            for dep in get_dependencies(dsk, k):
                dep_name = name(dep)
                if dep_name not in seen:
                    seen.add(dep_name)
                    attrs = data_attributes.get(dep, {}).copy()
                    attrs.setdefault("label", box_label(dep, verbose))
                    attrs.setdefault("shape", "box")
                    g.node(dep_name, **attrs)
                g.edge(dep_name, func_name)
                connected.add(dep_name)
                connected.add(func_name)

        elif ishashable(v) and v in dsk:
            v_name = name(v)
            g.edge(v_name, k_name)
            connected.add(v_name)
            connected.add(k_name)

        if (not collapse_outputs or k_name in connected) and k_name not in seen:
            seen.add(k_name)
            attrs = data_attributes.get(k, {}).copy()
            attrs.setdefault("label", box_label(k, verbose))
            attrs.setdefault("shape", "box")
            g.node(k_name, **attrs)
    return g


IPYTHON_IMAGE_FORMATS = frozenset(["jpeg", "png"])
IPYTHON_NO_DISPLAY_FORMATS = frozenset(["dot", "pdf"])


def _get_display_cls(format):
    """
    Get the appropriate IPython display class for `format`.

    Returns `IPython.display.SVG` if format=='svg', otherwise
    `IPython.display.Image`.

    If IPython is not importable, return dummy function that swallows its
    arguments and returns None.
    """
    dummy = lambda *args, **kwargs: None
    try:
        import IPython.display as display
    except ImportError:
        # Can't return a display object if no IPython.
        return dummy

    if format in IPYTHON_NO_DISPLAY_FORMATS:
        # IPython can't display this format natively, so just return None.
        return dummy
    elif format in IPYTHON_IMAGE_FORMATS:
        # Partially apply `format` so that `Image` and `SVG` supply a uniform
        # interface to the caller.
        return partial(display.Image, format=format)
    elif format == "svg":
        return display.SVG
    else:
        raise ValueError("Unknown format '%s' passed to `dot_graph`" % format)


def dot_graph(dsk, filename="mydask", format=None, **kwargs):
    """
    Render a task graph using dot.

    If `filename` is not None, write a file to disk with the specified name and extension.
    If no extension is specified, '.png' will be used by default.

    Parameters
    ----------
    dsk : dict
        The graph to display.
    filename : str or None, optional
        The name of the file to write to disk. If the provided `filename`
        doesn't include an extension, '.png' will be used by default.
        If `filename` is None, no file will be written, and we communicate
        with dot using only pipes.  Default is 'mydask'.
    format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
        Format in which to write output file.  Default is 'png'.
    **kwargs
        Additional keyword arguments to forward to `to_graphviz`.

    Returns
    -------
    result : None or IPython.display.Image or IPython.display.SVG  (See below.)

    Notes
    -----
    If IPython is installed, we return an IPython.display object in the
    requested format.  If IPython is not installed, we just return None.

    We always return None if format is 'pdf' or 'dot', because IPython can't
    display these formats natively. Passing these formats with filename=None
    will not produce any useful output.

    See Also
    --------
    dask.dot.to_graphviz
    """
    g = to_graphviz(dsk, **kwargs)
    return graphviz_to_file(g, filename, format)


def graphviz_to_file(g, filename, format):
    fmts = [".png", ".pdf", ".dot", ".svg", ".jpeg", ".jpg"]

    if (
        format is None
        and filename is not None
        and any(filename.lower().endswith(fmt) for fmt in fmts)
    ):
        filename, format = os.path.splitext(filename)
        format = format[1:].lower()

    if format is None:
        format = "png"

    data = g.pipe(format=format)
    if not data:
        raise RuntimeError(
            "Graphviz failed to properly produce an image. "
            "This probably means your installation of graphviz "
            "is missing png support. See: "
            "https://github.com/ContinuumIO/anaconda-issues/"
            "issues/485 for more information."
        )

    display_cls = _get_display_cls(format)

    if filename is None:
        return display_cls(data=data)

    full_filename = ".".join([filename, format])
    with open(full_filename, "wb") as f:
        f.write(data)

    return display_cls(filename=full_filename)


def _to_cytoscape_json(
    dsk,
    data_attributes=None,
    function_attributes=None,
    collapse_outputs=False,
    verbose=False,
    **kwargs,
):
    """
    Convert a dask graph to Cytoscape JSON:
    https://js.cytoscape.org/#notation/elements-json
    """
    nodes = []
    edges = []
    data = {"nodes": nodes, "edges": edges}

    data_attributes = data_attributes or {}
    function_attributes = function_attributes or {}

    seen = set()
    connected = set()

    for k, v in dsk.items():
        k_name = name(k)
        if istask(v):
            func_name = name((k, "function")) if not collapse_outputs else k_name
            if collapse_outputs or func_name not in seen:
                seen.add(func_name)
                attrs = function_attributes.get(k, {}).copy()
                nodes.append(
                    {
                        "data": {
                            "id": func_name,
                            "label": key_split(k),
                            "shape": "ellipse",
                            "color": "gray",
                            **attrs,
                        }
                    }
                )
            if not collapse_outputs:
                edges.append({"data": {"source": func_name, "target": k_name}})

                connected.add(func_name)
                connected.add(k_name)

            for dep in get_dependencies(dsk, k):
                dep_name = name(dep)
                if dep_name not in seen:
                    seen.add(dep_name)
                    attrs = data_attributes.get(dep, {}).copy()
                    nodes.append(
                        {
                            "data": {
                                "id": dep_name,
                                "label": box_label(dep, verbose),
                                "shape": "rectangle",
                                "color": "gray",
                                **attrs,
                            }
                        }
                    )
                edges.append(
                    {
                        "data": {
                            "source": dep_name,
                            "target": func_name,
                        }
                    }
                )
                connected.add(dep_name)
                connected.add(func_name)

        elif ishashable(v) and v in dsk:
            v_name = name(v)
            edges.append(
                {
                    "data": {
                        "source": v_name,
                        "target": k_name,
                    }
                }
            )
            connected.add(v_name)
            connected.add(k_name)

        if (not collapse_outputs or k_name in connected) and k_name not in seen:
            seen.add(k_name)
            attrs = data_attributes.get(k, {}).copy()
            nodes.append(
                {
                    "data": {
                        "id": k_name,
                        "label": box_label(k, verbose),
                        "shape": "rectangle",
                        "color": "gray",
                        **attrs,
                    }
                }
            )
    return data


def cytoscape_graph(
    dsk,
    filename: str | None = "mydask",
    format: str | None = None,
    *,
    rankdir: str = "BT",
    node_sep: float = 10,
    edge_sep: float = 10,
    spacing_factor: float = 1,
    node_style: dict[str, str] | None = None,
    edge_style: dict[str, str] | None = None,
    **kwargs,
):
    """
    Create an ipycytoscape widget for a dask graph.

    If `filename` is not None, write an HTML file to disk with the specified name.

    This uses the Cytoscape dagre layout algorithm. Options for that are documented here:
    https://github.com/cytoscape/cytoscape.js-dagre#api

    Parameters
    ----------
    dsk : dict
        The graph to display.
    filename : str or None, optional
        The name of the HTML file to write to disk.
    format : str, optional
        Not used in this engine.
    rankdir: str
        The direction in which to orient the visualization.
    node_sep: float
        The separation (in px) between nodes in the DAG layout.
    edge_sep: float
        The separation (in px) between edges in the DAG layout.
    spacing_factor: float
        An overall scaling factor to increase (>1) or decrease (<1) the spacing
        of the layout.
    node_style: dict[str, str], optional
        A dictionary of style attributes for nodes (refer to Cytoscape JSON docs
        for available options: https://js.cytoscape.org/#notation/elements-json)
    edge_style: dict[str, str], optional
        A dictionary of style attributes for edges (refer to Cytoscape JSON docs
        for available options: https://js.cytoscape.org/#notation/elements-json)
    **kwargs
        Additional keyword arguments to forward to `_to_cytoscape_json`.

    Returns
    -------
    result : ipycytoscape.CytoscapeWidget
    """
    ipycytoscape = import_required(
        "ipycytoscape",
        "Drawing dask graphs with the cytoscape engine requires the `ipycytoscape` "
        "python library.\n\n"
        "Please either conda or pip install as follows:\n\n"
        "  conda install ipycytoscape            # either conda install\n"
        "  python -m pip install ipycytoscape    # or pip install",
    )

    node_style = node_style or {}
    edge_style = edge_style or {}

    data = _to_cytoscape_json(dsk, **kwargs)
    # TODO: it's not easy to programmatically increase the height of the widget.
    # Ideally we would make it a bit bigger, but that will probably require upstreaming
    # some fixes.
    g = ipycytoscape.CytoscapeWidget(
        layout={"height": "400px"},
    )
    g.set_layout(
        name="dagre",
        rankDir=rankdir,
        nodeSep=node_sep,
        edgeSep=edge_sep,
        spacingFactor=spacing_factor,
        nodeDimensionsIncludeLabels=True,
    )
    g.graph.add_graph_from_json(
        data,
        directed=True,
    )
    g.set_style(
        [
            {
                "selector": "node",
                "style": {
                    "font-family": "helvetica",
                    "font-size": "24px",
                    "font-weight": "bold",
                    "color": "black",
                    "background-color": "#eee",
                    "border-color": "data(color)",
                    "border-width": 4,
                    "opacity": "1.0",
                    "text-valign": "center",
                    "text-halign": "center",
                    "label": "data(label)",
                    "shape": "data(shape)",
                    "width": 64,
                    "height": 64,
                    **node_style,
                },
            },
            {
                "selector": "edge",
                "style": {
                    "width": 8,
                    "line-color": "gray",
                    "target-arrow-shape": "triangle",
                    "target-arrow-color": "gray",
                    "curve-style": "bezier",
                    **edge_style,
                },
            },
        ],
    )
    # Tweak the zoom sensitivity
    z = g.zoom
    g.max_zoom = z * 2.0
    g.min_zoom = z / 10.0
    g.wheel_sensitivity = 0.1

    if filename is not None:
        from ipywidgets.embed import embed_minimal_html

        filename = filename if filename.endswith(".html") else filename + ".html"
        embed_minimal_html(filename, views=[g], title="Dask task graph")
    return g
