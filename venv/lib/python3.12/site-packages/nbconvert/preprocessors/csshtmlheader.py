"""Module that pre-processes the notebook for export to HTML."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import hashlib
import os

from jupyterlab_pygments import JupyterStyle  # type:ignore[import-untyped]
from pygments.style import Style
from traitlets import Type, Unicode, Union

from .base import Preprocessor

try:
    from notebook import DEFAULT_STATIC_FILES_PATH  # type:ignore[import-not-found]
except ImportError:
    DEFAULT_STATIC_FILES_PATH = None


class CSSHTMLHeaderPreprocessor(Preprocessor):
    """
    Preprocessor used to pre-process notebook for HTML output.  Adds IPython notebook
    front-end CSS and Pygments CSS to HTML output.
    """

    highlight_class = Unicode(".highlight", help="CSS highlight class identifier").tag(config=True)

    style = Union(
        [Unicode("default"), Type(klass=Style)],
        help="Name of the pygments style to use",
        default_value=JupyterStyle,
    ).tag(config=True)

    def __init__(self, *pargs, **kwargs):
        """Initialize the preprocessor."""
        Preprocessor.__init__(self, *pargs, **kwargs)
        self._default_css_hash = None

    def preprocess(self, nb, resources):
        """Fetch and add CSS to the resource dictionary

        Fetch CSS from IPython and Pygments to add at the beginning
        of the html files.  Add this css in resources in the
        "inlining.css" key

        Parameters
        ----------
        nb : NotebookNode
            Notebook being converted
        resources : dictionary
            Additional resources used in the conversion process.  Allows
            preprocessors to pass variables into the Jinja engine.
        """
        resources["inlining"] = {}
        resources["inlining"]["css"] = self._generate_header(resources)
        return nb, resources

    def _generate_header(self, resources):
        """
        Fills self.header with lines of CSS extracted from IPython
        and Pygments.
        """
        from pygments.formatters import HtmlFormatter  # noqa: PLC0415

        header = []

        formatter = HtmlFormatter(style=self.style)
        pygments_css = formatter.get_style_defs(self.highlight_class)
        header.append(pygments_css)

        # Load the user's custom CSS and IPython's default custom CSS.  If they
        # differ, assume the user has made modifications to his/her custom CSS
        # and that we should inline it in the nbconvert output.
        config_dir = resources["config_dir"]
        custom_css_filename = os.path.join(config_dir, "custom", "custom.css")
        if os.path.isfile(custom_css_filename):
            if DEFAULT_STATIC_FILES_PATH and self._default_css_hash is None:
                self._default_css_hash = self._hash(
                    os.path.join(DEFAULT_STATIC_FILES_PATH, "custom", "custom.css")
                )
            if self._hash(custom_css_filename) != self._default_css_hash:
                with open(custom_css_filename, encoding="utf-8") as f:
                    header.append(f.read())
        return header

    def _hash(self, filename):
        """Compute the hash of a file."""
        md5 = hashlib.md5()  # noqa: S324
        with open(filename, "rb") as f:
            md5.update(f.read())
        return md5.digest()
