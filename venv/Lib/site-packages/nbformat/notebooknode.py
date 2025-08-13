"""NotebookNode - adding attribute access to dicts"""

from __future__ import annotations

from collections.abc import Mapping

from ._struct import Struct


class NotebookNode(Struct):
    """A dict-like node with attribute-access"""

    def __setitem__(self, key, value):
        """Set an item on the notebook."""
        if isinstance(value, Mapping) and not isinstance(value, NotebookNode):
            value = from_dict(value)
        super().__setitem__(key, value)

    def update(self, *args, **kwargs):
        """
        A dict-like update method based on CPython's MutableMapping `update`
        method.
        """
        if len(args) > 1:
            raise TypeError("update expected at most 1 arguments, got %d" % len(args))
        if args:
            other = args[0]
            if isinstance(other, Mapping):  # noqa: SIM114
                for key in other:
                    self[key] = other[key]
            elif hasattr(other, "keys"):
                for key in other:
                    self[key] = other[key]
            else:
                for key, value in other:
                    self[key] = value
        for key, value in kwargs.items():
            self[key] = value


def from_dict(d):
    """Convert dict to dict-like NotebookNode

    Recursively converts any dict in the container to a NotebookNode.
    This does not check that the contents of the dictionary make a valid
    notebook or part of a notebook.
    """
    if isinstance(d, dict):
        return NotebookNode({k: from_dict(v) for k, v in d.items()})
    if isinstance(d, (tuple, list)):
        return [from_dict(i) for i in d]
    return d
