# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from collections import namedtuple

VersionInfo = namedtuple("VersionInfo", ["major", "minor", "micro", "releaselevel", "serial"])

# DO NOT EDIT THIS DIRECTLY!  It is managed by bumpversion
version_info = VersionInfo(4, 5, 4, "final", 0)

_specifier_ = {"alpha": "a", "beta": "b", "candidate": "rc", "final": ""}

__version__ = "{}.{}.{}{}".format(
    version_info.major,
    version_info.minor,
    version_info.micro,
    (
        ""
        if version_info.releaselevel == "final"
        else _specifier_[version_info.releaselevel] + str(version_info.serial)
    ),
)
