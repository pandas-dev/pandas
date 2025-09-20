"""Functionality related with node paths in a PyTables file.

Variables
=========

`__docformat`__
    The format of documentation strings in this module.

"""

import re
import keyword
import warnings

from .exceptions import NaturalNameWarning

__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""


_python_id_re = re.compile("^[a-zA-Z_][a-zA-Z0-9_]*$")
"""Python identifier regular expression."""

_reserved_id_re = re.compile("^_[cfgv]_")
"""PyTables reserved identifier regular expression.

- c: class variables
- f: class public methods
- g: class private methods
- v: instance variables
"""

_hidden_name_re = re.compile("^_[pi]_")
"""Nodes with a name *matching* this expression are considered hidden.

For instance, ``name`` would be visible while ``_i_name`` would not.
"""

_hidden_path_re = re.compile("/_[pi]_")
"""Nodes with a path *containing* this expression are considered hidden.

For instance, a node with a pathname like ``/a/b/c`` would be visible
while nodes with pathnames like ``/a/c/_i_x`` or ``/a/_p_x/y`` would
not.
"""

_warn_info = (
    "you will not be able to use natural naming to access this object; "
    "using ``getattr()`` will still work, though"
)
"""Warning printed when a name will not be reachable through natural naming"""


def check_attribute_name(name: str) -> None:
    """Check the validity of the `name` of an attribute in AttributeSet.

    If the name is not valid, a ``ValueError`` is raised.  If it is
    valid but it can not be used with natural naming, a
    `NaturalNameWarning` is issued.

    >>> warnings.simplefilter("ignore")
    >>> check_attribute_name('a')
    >>> check_attribute_name('a_b')
    >>> check_attribute_name('a:b')         # NaturalNameWarning
    >>> check_attribute_name('/a/b')        # NaturalNameWarning
    >>> check_attribute_name('/')           # NaturalNameWarning
    >>> check_attribute_name('.')           # NaturalNameWarning
    >>> check_attribute_name('__members__')
    Traceback (most recent call last):
     ...
    ValueError: ``__members__`` is not allowed as an object name
    >>> check_attribute_name(1)
    Traceback (most recent call last):
     ...
    TypeError: object name is not a string: 1
    >>> check_attribute_name('')
    Traceback (most recent call last):
     ...
    ValueError: the empty string is not allowed as an object name
    """
    if not isinstance(name, str):  # Python >= 2.3
        raise TypeError(f"object name is not a string: {name!r}")

    if name == "":
        raise ValueError("the empty string is not allowed as an object name")

    # Check whether `name` is a valid Python identifier.
    if not _python_id_re.match(name):
        warnings.warn(
            "object name is not a valid Python identifier: %r; "
            "it does not match the pattern ``%s``; %s"
            % (name, _python_id_re.pattern, _warn_info),
            NaturalNameWarning,
            stacklevel=2,
        )
        return

    # However, Python identifiers and keywords have the same form.
    if keyword.iskeyword(name):
        warnings.warn(
            f"object name is a Python keyword: {name!r}; {_warn_info}",
            NaturalNameWarning,
            stacklevel=2,
        )
        return

    # Still, names starting with reserved prefixes are not allowed.
    if _reserved_id_re.match(name):
        raise ValueError(
            "object name starts with a reserved prefix: %r; "
            "it matches the pattern ``%s``" % (name, _reserved_id_re.pattern)
        )

    # ``__members__`` is the only exception to that rule.
    if name == "__members__":
        raise ValueError("``__members__`` is not allowed as an object name")


def check_name_validity(name: str) -> None:
    """Check the validity of the `name` of a Node object.

    Validity of Node names is more limited than attribute names.

    If the name is not valid, a ``ValueError`` is raised.  If it is
    valid but it can not be used with natural naming, a
    `NaturalNameWarning` is issued.

    >>> warnings.simplefilter("ignore")
    >>> check_name_validity('a')
    >>> check_name_validity('a_b')
    >>> check_name_validity('a:b')          # NaturalNameWarning
    >>> check_name_validity('/a/b')
    Traceback (most recent call last):
     ...
    ValueError: the ``/`` character is not allowed in object names: '/a/b'
    >>> check_name_validity('.')
    Traceback (most recent call last):
     ...
    ValueError: ``.`` is not allowed as an object name
    >>> check_name_validity('')
    Traceback (most recent call last):
     ...
    ValueError: the empty string is not allowed as an object name

    """
    check_attribute_name(name)

    # Check whether `name` is a valid HDF5 name.
    # http://hdfgroup.org/HDF5/doc/UG/03_Model.html#Structure
    if name == ".":
        raise ValueError("``.`` is not allowed as an object name")
    elif "/" in name:
        raise ValueError(
            "the ``/`` character is not allowed " "in object names: %r" % name
        )


def join_path(parentpath: str, name: str) -> str:
    """Join a *canonical* `parentpath` with a *non-empty* `name`.

    .. versionchanged:: 3.0
       The *parentPath* parameter has been renamed into *parentpath*.

    >>> join_path('/', 'foo')
    '/foo'
    >>> join_path('/foo', 'bar')
    '/foo/bar'
    >>> join_path('/foo', '/foo2/bar')
    '/foo/foo2/bar'
    >>> join_path('/foo', '/')
    '/foo'

    """
    if name.startswith("./"):  # Support relative paths (mainly for links)
        name = name[2:]
    if parentpath == "/" and name.startswith("/"):
        pstr = "%s" % name
    elif parentpath == "/" or name.startswith("/"):
        pstr = f"{parentpath}{name}"
    else:
        pstr = f"{parentpath}/{name}"
    if pstr.endswith("/"):
        pstr = pstr[:-1]
    return pstr


def split_path(path: str) -> (str, str):
    """Split a *canonical* `path` into a parent path and a node name.

    The result is returned as a tuple.  The parent path does not
    include a trailing slash.

    >>> split_path('/')
    ('/', '')
    >>> split_path('/foo/bar')
    ('/foo', 'bar')

    """
    lastslash = path.rfind("/")
    ppath = path[:lastslash]
    name = path[lastslash + 1 :]

    if ppath == "":
        ppath = "/"

    return (ppath, name)


def isvisiblename(name: str) -> bool:
    """Return `True` if  `name` makes the named node visible."""
    return _hidden_name_re.match(name) is None


def isvisiblepath(path: str) -> bool:
    """Return `True` if `path` makes the named node visible."""
    return _hidden_path_re.search(path) is None


def _test() -> None:
    """Run ``doctest`` on this module."""
    import doctest

    doctest.testmod()


if __name__ == "__main__":
    _test()
