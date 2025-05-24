import pytest

from pandas.tests.apply.common import MockEngineDecorator


@pytest.fixture(params=[None, MockEngineDecorator])
def engine(request):
    return request.param
