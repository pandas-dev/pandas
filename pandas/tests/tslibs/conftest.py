import pytest

import pandas._testing as tm


@pytest.fixture(params=tm.utf8_locales, autouse=True)
def with_locale(request):
    with tm.set_locale(request.param):
        yield
