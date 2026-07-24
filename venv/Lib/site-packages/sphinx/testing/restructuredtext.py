from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.parsers import RSTParser
from sphinx.util.docutils import _parse_str_to_doctree

if TYPE_CHECKING:
    from docutils import nodes

    from sphinx.application import Sphinx


def parse(app: Sphinx, text: str, docname: str = 'index') -> nodes.document:
    """Parse a string as reStructuredText with Sphinx."""
    config = app.config
    env = app.env
    registry = app.registry
    srcdir = app.srcdir

    # Get settings
    settings_overrides = {
        'env': env,
        'gettext_compact': True,
        'input_encoding': 'utf-8',
        'output_encoding': 'unicode',
        'traceback': True,
    }

    # Create parser
    parser = RSTParser()
    parser._config = config
    parser._env = env

    env.current_document.docname = docname
    try:
        return _parse_str_to_doctree(
            text,
            filename=srcdir / f'{docname}.rst',
            default_settings=settings_overrides,
            env=env,
            parser=parser,
            transforms=registry.get_transforms(),
        )
    finally:
        env.current_document.docname = ''
