"""Allow `MathJax`_ to be used to display math in Sphinx's HTML writer.

This requires the MathJax JavaScript library on your webserver/computer.

.. _MathJax: https://www.mathjax.org/
"""

from __future__ import annotations

import json
from types import NoneType
from typing import TYPE_CHECKING, cast

from docutils import nodes

import sphinx
from sphinx.errors import ExtensionError
from sphinx.locale import _
from sphinx.util.math import get_node_equation_number

if TYPE_CHECKING:
    from typing import Any

    from sphinx.application import Sphinx
    from sphinx.builders.html import StandaloneHTMLBuilder
    from sphinx.util.typing import ExtensionMetadata
    from sphinx.writers.html5 import HTML5Translator

# more information for mathjax secure url is here:
# https://docs.mathjax.org/en/latest/web/start.html#using-mathjax-from-a-content-delivery-network-cdn
MATHJAX_URL = 'https://cdn.jsdelivr.net/npm/mathjax@4/tex-mml-chtml.js'

logger = sphinx.util.logging.getLogger(__name__)


def html_visit_math(self: HTML5Translator, node: nodes.math) -> None:
    self.body.append(
        self.starttag(node, 'span', '', CLASS='math notranslate nohighlight')
    )
    self.body.append(
        self.builder.config.mathjax_inline[0]
        + self.encode(node.astext())
        + self.builder.config.mathjax_inline[1]
        + '</span>'
    )
    raise nodes.SkipNode


def html_visit_displaymath(self: HTML5Translator, node: nodes.math_block) -> None:
    self.body.append(self.starttag(node, 'div', CLASS='math notranslate nohighlight'))
    if node.get('no-wrap', node.get('nowrap', False)):
        self.body.append(self.encode(node.astext()))
        self.body.append('</div>')
        raise nodes.SkipNode

    # necessary to e.g. set the id property correctly
    if node['number']:
        number = get_node_equation_number(self, node)
        self.body.append('<span class="eqno">(%s)' % number)
        self.add_permalink_ref(node, _('Link to this equation'))
        self.body.append('</span>')
    self.body.append(self.builder.config.mathjax_display[0])
    parts = [prt for prt in node.astext().split('\n\n') if prt.strip()]
    if len(parts) > 1:  # Add alignment if there are more than 1 equation
        self.body.append(r' \begin{align}\begin{aligned}')
    for i, part in enumerate(parts):
        part = self.encode(part)
        if r'\\' in part:
            self.body.append(r'\begin{split}' + part + r'\end{split}')
        else:
            self.body.append(part)
        if i < len(parts) - 1:  # append new line if not the last equation
            self.body.append(r'\\')
    if len(parts) > 1:  # Add alignment if there are more than 1 equation
        self.body.append(r'\end{aligned}\end{align} ')
    self.body.append(self.builder.config.mathjax_display[1])
    self.body.append('</div>\n')
    raise nodes.SkipNode


def install_mathjax(
    app: Sphinx,
    pagename: str,
    templatename: str,
    context: dict[str, Any],
    event_arg: Any,
) -> None:
    if app.builder.format != 'html':
        return
    if app.builder.math_renderer_name != 'mathjax':  # type: ignore[attr-defined]
        return
    if not app.config.mathjax_path:
        msg = 'mathjax_path config value must be set for the mathjax extension to work'
        raise ExtensionError(msg)

    builder = cast('StandaloneHTMLBuilder', app.builder)
    page_has_equations = context.get('has_maths_elements', False)

    # Enable mathjax only if equations exists
    if app.registry.html_assets_policy == 'always' or page_has_equations:
        if app.config.mathjax2_config:
            if app.config.mathjax_path == MATHJAX_URL:
                logger.warning(
                    'mathjax_config/mathjax2_config does not work '
                    'for the current MathJax version, use mathjax4_config instead'
                )
            body = f'MathJax.Hub.Config({json.dumps(app.config.mathjax2_config)})'
            builder.add_js_file('', type='text/x-mathjax-config', body=body)

        if app.config.mathjax3_config:
            body = f'window.MathJax = {json.dumps(app.config.mathjax3_config)}'
            builder.add_js_file('', body=body)

        if app.config.mathjax4_config:
            body = f'window.MathJax = {json.dumps(app.config.mathjax4_config)}'
            builder.add_js_file('', body=body)

        if app.config.mathjax_config_path:
            config_path = app.confdir / app.config.mathjax_config_path
            if not config_path.exists():
                msg = f'mathjax_config_path file not found: {config_path}'
                raise ExtensionError(msg)
            if not config_path.is_file() or config_path.suffix != '.js':
                msg = f'mathjax_config_path: expected a .js file, but got {config_path}'
                raise ExtensionError(msg)
            body = config_path.read_text(encoding='utf-8')
            builder.add_js_file('', body=body)

        options = {}
        if app.config.mathjax_options:
            options.update(app.config.mathjax_options)
        if 'async' not in options and 'defer' not in options:
            if app.config.mathjax2_config or app.config.mathjax_config:
                # Load old MathJax versions via the 'async' method
                options['async'] = 'async'
            else:
                # Load MathJax v3+ via the 'defer' method
                options['defer'] = 'defer'
        builder.add_js_file(app.config.mathjax_path, **options)


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_html_math_renderer(
        'mathjax',
        inline_renderers=(html_visit_math, None),
        block_renderers=(html_visit_displaymath, None),
    )

    app.add_config_value('mathjax_path', MATHJAX_URL, 'html', types=frozenset({str}))
    app.add_config_value('mathjax_options', {}, 'html', types=frozenset({dict}))
    app.add_config_value(
        'mathjax_inline', [r'\(', r'\)'], 'html', types=frozenset({list, tuple})
    )
    app.add_config_value(
        'mathjax_display', [r'\[', r'\]'], 'html', types=frozenset({list, tuple})
    )
    app.add_config_value(
        'mathjax_config', None, 'html', types=frozenset({dict, NoneType})
    )
    app.add_config_value(
        'mathjax2_config',
        lambda c: c.mathjax_config,
        'html',
        types=frozenset({dict, NoneType}),
    )
    app.add_config_value(
        'mathjax3_config', None, 'html', types=frozenset({dict, NoneType})
    )
    app.add_config_value(
        'mathjax4_config', None, 'html', types=frozenset({dict, NoneType})
    )
    app.add_config_value('mathjax_config_path', '', 'html', types=frozenset({str}))
    app.connect('html-page-context', install_mathjax)

    return {
        'version': sphinx.__display_version__,
        'parallel_read_safe': True,
    }
