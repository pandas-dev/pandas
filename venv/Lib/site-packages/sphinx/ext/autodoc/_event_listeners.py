"""Some useful event listener factories for autodoc-process-docstring."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Protocol

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc._property_types import _AutodocObjType

    class _AutodocProcessDocstringListener(Protocol):
        # parameter names are non-normative
        def __call__(
            self,
            app: Sphinx,
            obj_type: _AutodocObjType,
            full_name: str,
            obj: Any,
            options: Any,
            docstring_lines: list[str],
            /,
        ) -> None: ...

    class _AutodocBeforeProcessSignatureListener(Protocol):  # NoQA: PYI046
        # parameter names are non-normative
        def __call__(
            self,
            app: Sphinx,
            obj: Any,
            is_bound_method: bool,
            /,
        ) -> None: ...

    class _AutodocProcessSignatureListener(Protocol):  # NoQA: PYI046
        # parameter names are non-normative
        # returns: (args, retann) | None
        def __call__(
            self,
            app: Sphinx,
            obj_type: _AutodocObjType,
            full_name: str,
            obj: Any,
            options: Any,
            args: str | None,
            retann: str | None,
            /,
        ) -> tuple[str | None, str | None] | None: ...

    class _AutodocProcessBasesListener(Protocol):  # NoQA: PYI046
        # parameter names are non-normative
        def __call__(
            self,
            app: Sphinx,
            full_name: str,
            obj: Any,
            _unused: None,  # previously: options
            obj_bases: list[type],
            /,
        ) -> None: ...

    class _AutodocSkipMemberListener(Protocol):  # NoQA: PYI046
        # parameter names are non-normative
        # returns: skip_member
        def __call__(
            self,
            app: Sphinx,
            obj_type: _AutodocObjType,
            member_name: str,
            member_obj: Any,
            skip: bool,
            options: Any,
            /,
        ) -> bool | None: ...


def cut_lines(
    pre: int, post: int = 0, what: Sequence[str] | None = None
) -> _AutodocProcessDocstringListener:
    """Return a listener that removes the first *pre* and last *post*
    lines of every docstring.  If *what* is a sequence of strings,
    only docstrings of a type in *what* will be processed.

    Use like this (e.g. in the ``setup()`` function of :file:`conf.py`)::

       from sphinx.ext.autodoc import cut_lines

       app.connect('autodoc-process-docstring', cut_lines(4, what={'module'}))

    This can (and should) be used in place of :confval:`automodule_skip_lines`.
    """
    if not what:
        what_unique: frozenset[str] = frozenset()
    elif isinstance(what, str):  # strongly discouraged
        what_unique = frozenset({what})
    else:
        what_unique = frozenset(what)

    def process(
        app: Sphinx,
        what_: _AutodocObjType,
        name: str,
        obj: Any,
        options: Any,
        lines: list[str],
    ) -> None:
        if what_unique and what_ not in what_unique:
            return
        del lines[:pre]
        if post:
            # remove one trailing blank line.
            if lines and not lines[-1]:
                lines.pop(-1)
            del lines[-post:]
        # make sure there is a blank line at the end
        if lines and lines[-1]:
            lines.append('')

    return process


def between(
    marker: str,
    what: Sequence[str] | None = None,
    keepempty: bool = False,
    exclude: bool = False,
) -> _AutodocProcessDocstringListener:
    """Return a listener that either keeps, or if *exclude* is True excludes,
    lines between lines that match the *marker* regular expression.  If no line
    matches, the resulting docstring would be empty, so no change will be made
    unless *keepempty* is true.

    If *what* is a sequence of strings, only docstrings of a type in *what* will
    be processed.
    """
    marker_re = re.compile(marker)

    def process(
        app: Sphinx,
        what_: _AutodocObjType,
        name: str,
        obj: Any,
        options: Any,
        lines: list[str],
    ) -> None:
        if what and what_ not in what:
            return
        deleted = 0
        delete = not exclude
        orig_lines = lines.copy()
        for i, line in enumerate(orig_lines):
            if delete:
                lines.pop(i - deleted)
                deleted += 1
            if marker_re.match(line):
                delete = not delete
                if delete:
                    lines.pop(i - deleted)
                    deleted += 1
        if not lines and not keepempty:
            lines[:] = orig_lines
        # make sure there is a blank line at the end
        if lines and lines[-1]:
            lines.append('')

    return process
