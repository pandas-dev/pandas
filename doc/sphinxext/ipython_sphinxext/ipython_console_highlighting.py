"""
reST directive for syntax-highlighting ipython interactive sessions.

"""

from IPython.lib.lexers import IPyLexer
from sphinx import highlighting


def patch_IPythonConsoleLexer():
    """Patch lexers.IPythonConsoleLexer.ipytb_start.

    This extension uses IPyLexer to highlight the exception. An exception is
    recognized by the traceback. The start of a traceback is found using the
    ipytb_start regex. This regex has to be patched to also find exceptions
    without a preceding traceback.
    """
    import builtins
    import re

    from IPython.lib import lexers

    exceptions = [
        name
        for name, value in builtins.__dict__.items()
        if isinstance(value, type)
        and issubclass(value, Exception)
        and not issubclass(value, Warning)
    ]
    lexers.IPythonConsoleLexer.ipytb_start = re.compile(
        r"^(\^C)?(-+\n)|^(  File)(.*)(, line )(\d+\n)|^("
        + "|".join(exceptions)
        + r"): \S.*\n"
    )


patch_IPythonConsoleLexer()


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

ipy2 = IPyLexer(python3=False)
ipy3 = IPyLexer(python3=True)

highlighting.lexers["ipython"] = ipy2
highlighting.lexers["ipython2"] = ipy2
highlighting.lexers["ipython3"] = ipy3
