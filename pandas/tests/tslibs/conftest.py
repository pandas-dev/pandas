import pytest

import pandas._testing as tm


@pytest.fixture(params=tm.TESTING_LOCALES, autouse=True)
def with_locale(request):
    with tm.set_locale(request.param):
        yield
