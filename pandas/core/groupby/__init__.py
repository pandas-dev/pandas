from __future__ import annotations

from typing import Any

from pandas._core import groupby as groupby_
from pandas.core.common import _depr_core


def __getattr__(attr_name: str) -> Any:
    attr = getattr(groupby_, attr_name)
    _depr_core()
    return attr
