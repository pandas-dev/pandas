"""Backport of Python 3.11's contextlib.chdir."""

import os
from contextlib import AbstractContextManager


class chdir(AbstractContextManager):  # type:ignore[type-arg]
    """Non thread-safe context manager to change the current working directory."""

    def __init__(self, path):
        """Initialize the manager."""
        self.path = path
        self._old_cwd = []

    def __enter__(self):
        """Enter the context."""
        self._old_cwd.append(os.getcwd())
        os.chdir(self.path)

    def __exit__(self, *excinfo):
        """Exit the context."""
        os.chdir(self._old_cwd.pop())
