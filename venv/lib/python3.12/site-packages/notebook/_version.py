"""Version info for notebook."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import re
from collections import namedtuple

# Use "hatch version xx.yy.zz" to handle version changes
__version__ = "7.4.5"

# PEP440 version parser
_version_regex = re.compile(
    r"""
  (?P<major>\d+)
  \.
  (?P<minor>\d+)
  \.
  (?P<micro>\d+)
  (?P<releaselevel>((a|b|rc|\.dev)))?
  (?P<serial>\d+)?
  """,
    re.VERBOSE,
)

_version_fields = _version_regex.match(__version__).groupdict()  # type:ignore[union-attr]

VersionInfo = namedtuple("VersionInfo", ["major", "minor", "micro", "releaselevel", "serial"])  # noqa: PYI024

version_info = VersionInfo(
    *[
        field
        for field in (
            int(_version_fields["major"]),
            int(_version_fields["minor"]),
            int(_version_fields["micro"]),
            _version_fields["releaselevel"] or "",
            _version_fields["serial"] or "",
        )
    ]
)
