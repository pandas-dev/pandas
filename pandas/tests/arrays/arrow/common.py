from __future__ import annotations

import pytest

from pandas.compat import (
    is_platform_windows,
)
from pandas.compat.pyarrow import pa_version_under22p0

pa = pytest.importorskip("pyarrow")


def _require_timezone_database(request):
    if is_platform_windows() and pa_version_under22p0:
        mark = pytest.mark.xfail(
            raises=pa.ArrowInvalid,
            reason=(
                "TODO: Set ARROW_TIMEZONE_DATABASE environment variable "
                "on CI to path to the tzdata for pyarrow."
            ),
        )
        request.applymarker(mark)
