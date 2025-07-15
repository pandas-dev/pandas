"""
reST directive for syntax-highlighting ipython interactive sessions.

"""

from sphinx import highlighting
from ipython_pygments_lexers import IPyLexer


def setup(app):
    """Setup as a sphinx extension."""

    # This is only a lexer, so adding it below to pygments appears sufficient.
    # But if somebody knows what the right API usage should be to do that via
    # sphinx, by all means fix it here.  At least having this setup.py
    # suppresses the sphinx warning we'd get without it.
    metadata = {"parallel_read_safe": True, "parallel_write_safe": True}
    return metadata


# Register the extension as a valid pygments lexer.
# Alternatively, we could register the lexer with pygments instead. This would
# require using setuptools entrypoints: http://pygments.org/docs/plugins

ipy3 = IPyLexer()

highlighting.lexers["ipython"] = ipy3
highlighting.lexers["ipython3"] = ipy3
