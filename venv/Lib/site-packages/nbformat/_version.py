"""The version information for nbformat."""

# Use "hatchling version xx.yy.zz" to handle version changes
from __future__ import annotations

import re
from importlib.metadata import version

__version__ = version("nbformat") or "0.0.0"

# matches tbump regex in pyproject.toml
_version_regex = re.compile(
    r"""
  (?P<major>\d+)
  \.
  (?P<minor>\d+)
  \.
  (?P<patch>\d+)
  (?P<pre>((a|b|rc)\d+))?
  (\.
    (?P<dev>dev\d*)
  )?
  """,
    re.VERBOSE,
)

_version_fields = _version_regex.match(__version__).groupdict()  # type:ignore[union-attr]
version_info = tuple(
    field
    for field in (
        int(_version_fields["major"]),
        int(_version_fields["minor"]),
        int(_version_fields["patch"]),
        _version_fields["pre"],
        _version_fields["dev"],
    )
    if field is not None
)
