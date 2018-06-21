import pytest

CALL_COUNT = 0 


@pytest.fixture(scope="module")
def fixture(request, datapath):
    global CALL_COUNT
    CALL_COUNT += 1 

    return request.param


def pytest_generate_tests(metafunc):
    if "fixture" in metafunc.fixturenames:
        metafunc.parametrize("fixture", ["foo"], indirect=True, scope="module")


@pytest.mark.parametrize("param", ["bar", "zaz"])
def test_1(fixture, param):
    global CALL_COUNT
    assert CALL_COUNT == 1
