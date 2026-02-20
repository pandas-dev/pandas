"""This preprocessor detect cells using a different language through
magic extensions such as `%%R` or `%%octave`. Cell's metadata is marked
so that the appropriate highlighter can be used in the `highlight`
filter.
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import re

from traitlets import Dict

from .base import Preprocessor


class HighlightMagicsPreprocessor(Preprocessor):
    """
    Detects and tags code cells that use a different languages than Python.
    """

    # list of magic language extensions and their associated pygment lexers
    default_languages = Dict(
        {
            "%%R": "r",
            "%%bash": "bash",
            "%%cython": "cython",
            "%%javascript": "javascript",
            "%%julia": "julia",
            "%%latex": "latex",
            "%%octave": "octave",
            "%%perl": "perl",
            "%%ruby": "ruby",
            "%%sh": "sh",
            "%%sql": "sql",
        }
    )

    # user defined language extensions
    languages = Dict(
        help=(
            "Syntax highlighting for magic's extension languages. "
            "Each item associates a language magic extension such as %%R, "
            "with a pygments lexer such as r."
        )
    ).tag(config=True)

    def __init__(self, config=None, **kw):
        """Public constructor"""

        super().__init__(config=config, **kw)

        # Update the default languages dict with the user configured ones
        self.default_languages.update(self.languages)

        # build a regular expression to catch language extensions and choose
        # an adequate pygments lexer
        any_language = "|".join(self.default_languages.keys())
        self.re_magic_language = re.compile(rf"^\s*({any_language})\s+")

    def which_magic_language(self, source):
        """
        When a cell uses another language through a magic extension,
        the other language is returned.
        If no language magic is detected, this function returns None.

        Parameters
        ----------
        source: str
            Source code of the cell to highlight
        """

        m = self.re_magic_language.match(source)

        if m:
            # By construction of the re, the matched language must be in the
            # languages dictionary
            return self.default_languages[m.group(1)]
        return None

    def preprocess_cell(self, cell, resources, cell_index):
        """
        Tags cells using a magic extension language

        Parameters
        ----------
        cell : NotebookNode cell
            Notebook cell being processed
        resources : dictionary
            Additional resources used in the conversion process.  Allows
            preprocessors to pass variables into the Jinja engine.
        cell_index : int
            Index of the cell being processed (see base.py)
        """

        # Only tag code cells
        if cell.cell_type == "code":
            magic_language = self.which_magic_language(cell.source)
            if magic_language:
                cell["metadata"]["magics_language"] = magic_language
        return cell, resources
