# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
#
#
# Parts of this code is from IPyVolume (24.05.2017), used here under
# this copyright and license with permission from the author
# (see https://github.com/jupyter-widgets/ipywidgets/pull/1387)

"""
Functions for generating embeddable HTML/javascript of a widget.
"""

import json
import re
from .widgets import Widget, DOMWidget, widget as widget_module
from .widgets.widget_link import Link
from .widgets.docutils import doc_subst
from ._version import __html_manager_version__

snippet_template = """
{load}
<script type="application/vnd.jupyter.widget-state+json">
{json_data}
</script>
{widget_views}
"""

load_template = """<script src="{embed_url}"{use_cors}></script>"""

load_requirejs_template = """
<!-- Load require.js. Delete this if your page already loads require.js -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js" integrity="sha256-Ae2Vz/4ePdIu6ZyI/5ZGsYnb+m0JlOmKPjt6XZ9JJkA=" crossorigin="anonymous"></script>
<script src="{embed_url}"{use_cors}></script>
"""

requirejs_snippet_template = """
<script type="application/vnd.jupyter.widget-state+json">
{json_data}
</script>
{widget_views}
"""



html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
</head>
<body>
{snippet}
</body>
</html>
"""

widget_view_template = """<script type="application/vnd.jupyter.widget-view+json">
{view_spec}
</script>"""

DEFAULT_EMBED_SCRIPT_URL = 'https://cdn.jsdelivr.net/npm/@jupyter-widgets/html-manager@%s/dist/embed.js'%__html_manager_version__
DEFAULT_EMBED_REQUIREJS_URL = 'https://cdn.jsdelivr.net/npm/@jupyter-widgets/html-manager@%s/dist/embed-amd.js'%__html_manager_version__

_doc_snippets = {}
_doc_snippets['views_attribute'] = """
    views: widget or collection of widgets or None
        The widgets to include views for. If None, all DOMWidgets are
        included (not just the displayed ones).
"""
_doc_snippets['embed_kwargs'] = """
    drop_defaults: boolean
        Whether to drop default values from the widget states.
    state: dict or None (default)
        The state to include. When set to None, the state of all widgets
        know to the widget manager is included. Otherwise it uses the
        passed state directly. This allows for end users to include a
        smaller state, under the responsibility that this state is
        sufficient to reconstruct the embedded views.
    indent: integer, string or None
        The indent to use for the JSON state dump. See `json.dumps` for
        full description.
    embed_url: string or None
        Allows for overriding the URL used to fetch the widget manager
        for the embedded code. This defaults (None) to a `jsDelivr` CDN url.
    requirejs: boolean (True)
        Enables the requirejs-based embedding, which allows for custom widgets.
        If True, the embed_url should point to an AMD module.
    cors: boolean (True)
        If True avoids sending user credentials while requesting the scripts.
        When opening an HTML file from disk, some browsers may refuse to load
        the scripts.
"""


def _find_widget_refs_by_state(widget, state):
    """Find references to other widgets in a widget's state"""
    # Copy keys to allow changes to state during iteration:
    keys = tuple(state.keys())
    for key in keys:
        value = getattr(widget, key)
        # Trivial case: Direct references to other widgets:
        if isinstance(value, Widget):
            yield value
        # Also check for buried references in known, JSON-able structures
        # Note: This might miss references buried in more esoteric structures
        elif isinstance(value, (list, tuple)):
            for item in value:
                if isinstance(item, Widget):
                    yield item
        elif isinstance(value, dict):
            for item in value.values():
                if isinstance(item, Widget):
                    yield item


def _get_recursive_state(widget, store=None, drop_defaults=False):
    """Gets the embed state of a widget, and all other widgets it refers to as well"""
    if store is None:
        store = dict()
    state = widget._get_embed_state(drop_defaults=drop_defaults)
    store[widget.model_id] = state

    # Loop over all values included in state (i.e. don't consider excluded values):
    for ref in _find_widget_refs_by_state(widget, state['state']):
        if ref.model_id not in store:
            _get_recursive_state(ref, store, drop_defaults=drop_defaults)
    return store


def add_resolved_links(store, drop_defaults):
    """Adds the state of any link models between two models in store"""
    for widget_id, widget in widget_module._instances.items(): # go over all widgets
        if isinstance(widget, Link) and widget_id not in store:
            if widget.source[0].model_id in store and widget.target[0].model_id in store:
                store[widget.model_id] = widget._get_embed_state(drop_defaults=drop_defaults)


def dependency_state(widgets, drop_defaults=True):
    """Get the state of all widgets specified, and their dependencies.

    This uses a simple dependency finder, including:
     - any widget directly referenced in the state of an included widget
     - any widget in a list/tuple attribute in the state of an included widget
     - any widget in a dict attribute in the state of an included widget
     - any jslink/jsdlink between two included widgets
    What this alogorithm does not do:
     - Find widget references in nested list/dict structures
     - Find widget references in other types of attributes

    Note that this searches the state of the widgets for references, so if
    a widget reference is not included in the serialized state, it won't
    be considered as a dependency.

    Parameters
    ----------
    widgets: single widget or list of widgets.
       This function will return the state of every widget mentioned
       and of all their dependencies.
    drop_defaults: boolean
        Whether to drop default values from the widget states.

    Returns
    -------
    A dictionary with the state of the widgets and any widget they
    depend on.
    """
    # collect the state of all relevant widgets
    if widgets is None:
        # Get state of all widgets, no smart resolution needed.
        state = Widget.get_manager_state(drop_defaults=drop_defaults, widgets=None)['state']
    else:
        try:
            widgets[0]
        except (IndexError, TypeError):
            widgets = [widgets]
        state = {}
        for widget in widgets:
            _get_recursive_state(widget, state, drop_defaults)
        # Add any links between included widgets:
        add_resolved_links(state, drop_defaults)
    return state


@doc_subst(_doc_snippets)
def embed_data(views, drop_defaults=True, state=None):
    """Gets data for embedding.

    Use this to get the raw data for embedding if you have special
    formatting needs.

    Parameters
    ----------
    {views_attribute}
    drop_defaults: boolean
        Whether to drop default values from the widget states.
    state: dict or None (default)
        The state to include. When set to None, the state of all widgets
        know to the widget manager is included. Otherwise it uses the
        passed state directly. This allows for end users to include a
        smaller state, under the responsibility that this state is
        sufficient to reconstruct the embedded views.

    Returns
    -------
    A dictionary with the following entries:
        manager_state: dict of the widget manager state data
        view_specs: a list of widget view specs
    """
    if views is None:
        views = [w for w in widget_module._instances.values() if isinstance(w, DOMWidget)]
    else:
        try:
            views[0]
        except (IndexError, TypeError):
            views = [views]

    if state is None:
        # Get state of all known widgets
        state = Widget.get_manager_state(drop_defaults=drop_defaults, widgets=None)['state']

    # Rely on ipywidget to get the default values
    json_data = Widget.get_manager_state(widgets=[])
    # but plug in our own state
    json_data['state'] = state

    view_specs = [w.get_view_spec() for w in views]

    return dict(manager_state=json_data, view_specs=view_specs)

script_escape_re = re.compile(r'<(script|/script|!--)', re.IGNORECASE)
def escape_script(s):
    """Escape a string that will be the content of an HTML script tag.

    We replace the opening bracket of <script, </script, and <!-- with the unicode
    equivalent. This is inspired by the documentation for the script tag at
    https://html.spec.whatwg.org/multipage/scripting.html#restrictions-for-contents-of-script-elements

    We only replace these three cases so that most html or other content
    involving `<` is readable.
    """
    return script_escape_re.sub(r'\\u003c\1', s)

@doc_subst(_doc_snippets)
def embed_snippet(views,
                  drop_defaults=True,
                  state=None,
                  indent=2,
                  embed_url=None,
                  requirejs=True,
                  cors=True
                 ):
    """Return a snippet that can be embedded in an HTML file.

    Parameters
    ----------
    {views_attribute}
    {embed_kwargs}

    Returns
    -------
    A unicode string with an HTML snippet containing several `<script>` tags.
    """

    data = embed_data(views, drop_defaults=drop_defaults, state=state)

    widget_views = '\n'.join(
        widget_view_template.format(view_spec=escape_script(json.dumps(view_spec)))
        for view_spec in data['view_specs']
    )

    if embed_url is None:
        embed_url = DEFAULT_EMBED_REQUIREJS_URL if requirejs else DEFAULT_EMBED_SCRIPT_URL

    load = load_requirejs_template if requirejs else load_template

    use_cors = ' crossorigin="anonymous"' if cors else ''
    values = {
        'load': load.format(embed_url=embed_url, use_cors=use_cors),
        'json_data': escape_script(json.dumps(data['manager_state'], indent=indent)),
        'widget_views': widget_views,
    }

    return snippet_template.format(**values)


@doc_subst(_doc_snippets)
def embed_minimal_html(fp, views, title='IPyWidget export', template=None, **kwargs):
    """Write a minimal HTML file with widget views embedded.

    Parameters
    ----------
    fp: filename or file-like object
        The file to write the HTML output to.
    {views_attribute}
    title: title of the html page.
    template: Template in which to embed the widget state.
        This should be a Python string with placeholders
        `{{title}}` and `{{snippet}}`. The `{{snippet}}` placeholder
        will be replaced by all the widgets.
    {embed_kwargs}
    """
    snippet = embed_snippet(views, **kwargs)

    values = {
        'title': title,
        'snippet': snippet,
    }
    if template is None:
        template = html_template

    html_code = template.format(**values)

    # Check if fp is writable:
    if hasattr(fp, 'write'):
        fp.write(html_code)
    else:
        # Assume fp is a filename:
        with open(fp, "w") as f:
            f.write(html_code)
