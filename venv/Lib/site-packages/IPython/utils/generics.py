"""Generic functions for extending IPython."""

from functools import singledispatch

from IPython.core.error import TryNext


@singledispatch
def inspect_object(obj):
    """Called when you do obj?"""
    raise TryNext


@singledispatch
def complete_object(obj, prev_completions):
    """Custom completer dispatching for python objects.

    Parameters
    ----------
    obj : object
        The object to complete.
    prev_completions : list
        List of attributes discovered so far.
    This should return the list of attributes in obj. If you only wish to
    add to the attributes already discovered normally, return
    own_attrs + prev_completions.
    """
    raise TryNext
