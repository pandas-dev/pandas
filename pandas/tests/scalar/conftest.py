from __future__ import annotations

import re

import pytest

from pandas._libs.tslibs import OutOfBoundsTimedelta


@pytest.fixture()
def timedelta_overflow() -> dict:
    """
    The expected message and exception when Timedelta ops overflow.
    """
    msg = re.escape(
        "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
    )
    return {"expected_exception": OutOfBoundsTimedelta, "match": msg}
