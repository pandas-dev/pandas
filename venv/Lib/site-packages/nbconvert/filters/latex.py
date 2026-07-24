"""Latex filters.

Module of useful filters for processing Latex within Jinja latex templates.
"""
# -----------------------------------------------------------------------------
# Copyright (c) 2013, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
import re

# -----------------------------------------------------------------------------
# Globals and constants
# -----------------------------------------------------------------------------

LATEX_RE_SUBS = ((re.compile(r"\.\.\.+"), r"{\\ldots}"),)

# Latex substitutions for escaping latex.
# see: http://stackoverflow.com/questions/16259923/how-can-i-escape-latex-special-characters-inside-django-templates

LATEX_SUBS = {
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\^{}",
    "\\": r"\textbackslash{}",
}


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

__all__ = ["escape_latex"]


def escape_latex(text):
    """
    Escape characters that may conflict with latex.

    Parameters
    ----------
    text : str
        Text containing characters that may conflict with Latex
    """
    text = "".join(LATEX_SUBS.get(c, c) for c in text)
    for pattern, replacement in LATEX_RE_SUBS:
        text = pattern.sub(replacement, text)

    return text
