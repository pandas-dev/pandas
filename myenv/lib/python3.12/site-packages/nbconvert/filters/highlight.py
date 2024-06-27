"""
Module containing filter functions that allow code to be highlighted
from within Jinja templates.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

# pygments must not be imported at the module level
# because errors should be raised at runtime if it's actually needed,
# not import time, when it may not be needed.

from html import escape
from warnings import warn

from traitlets import Dict, observe

from nbconvert.utils.base import NbConvertBase

MULTILINE_OUTPUTS = ["text", "html", "svg", "latex", "javascript", "json"]

__all__ = ["Highlight2HTML", "Highlight2Latex"]


class Highlight2HTML(NbConvertBase):
    """Convert highlighted code to html."""

    extra_formatter_options = Dict(
        {},
        help="""
        Extra set of options to control how code is highlighted.

        Passed through to the pygments' HtmlFormatter class.
        See available list in https://pygments.org/docs/formatters/#HtmlFormatter
        """,
        config=True,
    )

    def __init__(self, pygments_lexer=None, **kwargs):
        """Initialize the converter."""
        self.pygments_lexer = pygments_lexer or "ipython3"
        super().__init__(**kwargs)

    @observe("default_language")
    def _default_language_changed(self, change):
        warn(
            "Setting default_language in config is deprecated as of 5.0, "
            "please use language_info metadata instead.",
            stacklevel=2,
        )
        self.pygments_lexer = change["new"]

    def __call__(self, source, language=None, metadata=None):
        """
        Return a syntax-highlighted version of the input source as html output.

        Parameters
        ----------
        source : str
            source of the cell to highlight
        language : str
            language to highlight the syntax of
        metadata : NotebookNode cell metadata
            metadata of the cell to highlight
        """
        from pygments.formatters import HtmlFormatter

        if not language:
            language = self.pygments_lexer

        return _pygments_highlight(
            source if len(source) > 0 else " ",
            # needed to help post processors:
            HtmlFormatter(
                cssclass=escape(f" highlight hl-{language}"), **self.extra_formatter_options
            ),
            language,
            metadata,
        )


class Highlight2Latex(NbConvertBase):
    """Convert highlighted code to latex."""

    extra_formatter_options = Dict(
        {},
        help="""
        Extra set of options to control how code is highlighted.

        Passed through to the pygments' LatexFormatter class.
        See available list in https://pygments.org/docs/formatters/#LatexFormatter
        """,
        config=True,
    )

    def __init__(self, pygments_lexer=None, **kwargs):
        """Initialize the converter."""
        self.pygments_lexer = pygments_lexer or "ipython3"
        super().__init__(**kwargs)

    @observe("default_language")
    def _default_language_changed(self, change):
        warn(
            "Setting default_language in config is deprecated as of 5.0, "
            "please use language_info metadata instead.",
            stacklevel=2,
        )
        self.pygments_lexer = change["new"]

    def __call__(self, source, language=None, metadata=None, strip_verbatim=False):
        """
        Return a syntax-highlighted version of the input source as latex output.

        Parameters
        ----------
        source : str
            source of the cell to highlight
        language : str
            language to highlight the syntax of
        metadata : NotebookNode cell metadata
            metadata of the cell to highlight
        strip_verbatim : bool
            remove the Verbatim environment that pygments provides by default
        """
        from pygments.formatters import LatexFormatter

        if not language:
            language = self.pygments_lexer

        latex = _pygments_highlight(
            source, LatexFormatter(**self.extra_formatter_options), language, metadata
        )
        if strip_verbatim:
            latex = latex.replace(r"\begin{Verbatim}[commandchars=\\\{\}]" + "\n", "")
            return latex.replace("\n\\end{Verbatim}\n", "")
        return latex


def _pygments_highlight(source, output_formatter, language="ipython", metadata=None):
    """
    Return a syntax-highlighted version of the input source

    Parameters
    ----------
    source : str
        source of the cell to highlight
    output_formatter : Pygments formatter
    language : str
        language to highlight the syntax of
    metadata : NotebookNode cell metadata
        metadata of the cell to highlight
    """
    from pygments import highlight
    from pygments.lexers import get_lexer_by_name
    from pygments.util import ClassNotFound

    # If the cell uses a magic extension language,
    # use the magic language instead.
    if language.startswith("ipython") and metadata and "magics_language" in metadata:
        language = metadata["magics_language"]

    lexer = None
    if language == "ipython2":
        try:
            from IPython.lib.lexers import IPythonLexer
        except ImportError:
            warn("IPython lexer unavailable, falling back on Python", stacklevel=2)
            language = "python"
        else:
            lexer = IPythonLexer()
    elif language == "ipython3":
        try:
            from IPython.lib.lexers import IPython3Lexer
        except ImportError:
            warn("IPython3 lexer unavailable, falling back on Python 3", stacklevel=2)
            language = "python3"
        else:
            lexer = IPython3Lexer()

    if lexer is None:
        try:
            lexer = get_lexer_by_name(language, stripall=True)
        except ClassNotFound:
            warn("No lexer found for language %r. Treating as plain text." % language, stacklevel=2)
            from pygments.lexers.special import TextLexer

            lexer = TextLexer()

    return highlight(source, lexer, output_formatter)
