# -*- coding: utf-8 -*-
"""
The IPython lexers are now a separate package, ipython-pygments-lexers.

Importing from here is deprecated and may break in the future.
"""
# -----------------------------------------------------------------------------
# Copyright (c) 2013, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from ipython_pygments_lexers import (
    IPythonLexer,
    IPython3Lexer,
    IPythonPartialTracebackLexer,
    IPythonTracebackLexer,
    IPythonConsoleLexer,
    IPyLexer,
)


__all__ = [
    "IPython3Lexer",
    "IPythonLexer",
    "IPythonPartialTracebackLexer",
    "IPythonTracebackLexer",
    "IPythonConsoleLexer",
    "IPyLexer",
]
