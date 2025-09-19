"""Sentinel class for constants with useful reprs"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations


class Sentinel:
    """Sentinel class for constants with useful reprs"""

    def __init__(self, name, module, docstring=None):
        """Initialize the sentinel."""
        self.name = name
        self.module = module
        if docstring:
            self.__doc__ = docstring

    def __repr__(self):
        """The string repr for the sentinel."""
        return str(self.module) + "." + self.name
