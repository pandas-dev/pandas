from __future__ import annotations

import os

from setuptools_scm import _types as _t


def norm_real(path: _t.PathT) -> str:
    return os.path.normcase(os.path.realpath(path))
