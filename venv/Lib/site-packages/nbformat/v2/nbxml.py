"""REMOVED: Read and write notebook files as XML."""

from __future__ import annotations

REMOVED_MSG = """\
Reading notebooks as XML has been removed to harden security and avoid
possible denial-of-service attacks.

The XML notebook format was deprecated before the Jupyter (previously IPython)
Notebook was ever released. We are not aware of anyone using it, so we have
removed it.

If you were using this code, and you need to continue using it, feel free to
fork an earlier version of the nbformat package and maintain it yourself.
The issue which prompted this removal is:

https://github.com/jupyter/nbformat/issues/132
"""


def reads(s, **kwargs):
    """REMOVED"""
    raise Exception(REMOVED_MSG)


def read(fp, **kwargs):
    """REMOVED"""
    raise Exception(REMOVED_MSG)


def to_notebook(root, **kwargs):
    """REMOVED"""
    raise Exception(REMOVED_MSG)
