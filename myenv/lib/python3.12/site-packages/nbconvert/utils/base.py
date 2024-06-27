"""Global configuration class."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import List, Unicode
from traitlets.config.configurable import LoggingConfigurable


class NbConvertBase(LoggingConfigurable):
    """Global configurable class for shared config

    Useful for display data priority that might be used by many transformers
    """

    display_data_priority = List(
        [
            "text/html",
            "application/pdf",
            "text/latex",
            "image/svg+xml",
            "image/png",
            "image/jpeg",
            "text/markdown",
            "text/plain",
        ],
        help="""
            An ordered list of preferred output type, the first
            encountered will usually be used when converting discarding
            the others.
            """,
    ).tag(config=True)

    default_language = Unicode(
        "ipython",
        help="Deprecated default highlight language as of 5.0, please use language_info metadata instead",
    ).tag(config=True)
