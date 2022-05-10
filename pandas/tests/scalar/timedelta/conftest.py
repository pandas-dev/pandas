import re

import pytest


@pytest.fixture(name="td_overflow_msg")
def fixture_td_overflow_msg() -> str:
    return re.escape(
        "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
    )
