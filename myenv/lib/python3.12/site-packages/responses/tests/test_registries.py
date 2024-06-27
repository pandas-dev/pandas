import pytest
import requests
from requests.exceptions import ConnectionError

import responses
from responses import registries
from responses.registries import OrderedRegistry
from responses.tests.test_responses import assert_reset


def test_set_registry_not_empty():
    class CustomRegistry(registries.FirstMatchRegistry):
        pass

    @responses.activate
    def run():
        url = "http://fizzbuzz/foo"
        responses.add(method=responses.GET, url=url)
        with pytest.raises(AttributeError) as excinfo:
            responses.mock._set_registry(CustomRegistry)
        msg = str(excinfo.value)
        assert "Cannot replace Registry, current registry has responses" in msg

    run()
    assert_reset()


def test_set_registry():
    class CustomRegistry(registries.FirstMatchRegistry):
        pass

    @responses.activate(registry=CustomRegistry)
    def run_with_registry():
        assert type(responses.mock.get_registry()) == CustomRegistry

    @responses.activate
    def run():
        # test that registry does not leak to another test
        assert type(responses.mock.get_registry()) == registries.FirstMatchRegistry

    run_with_registry()
    run()
    assert_reset()


def test_set_registry_reversed():
    """See https://github.com/getsentry/responses/issues/563"""

    class CustomRegistry(registries.FirstMatchRegistry):
        pass

    @responses.activate
    def run():
        # test that registry does not leak to another test
        assert type(responses.mock.get_registry()) == registries.FirstMatchRegistry

    @responses.activate(registry=CustomRegistry)
    def run_with_registry():
        assert type(responses.mock.get_registry()) == CustomRegistry

    run()
    run_with_registry()
    assert_reset()


@pytest.mark.asyncio
async def test_registry_async():  # type: ignore[misc]
    class CustomRegistry(registries.FirstMatchRegistry):
        pass

    @responses.activate
    async def run():
        # test that registry does not leak to another test
        assert type(responses.mock.get_registry()) == registries.FirstMatchRegistry

    @responses.activate(registry=CustomRegistry)
    async def run_with_registry():
        assert type(responses.mock.get_registry()) == CustomRegistry

    await run()
    await run_with_registry()
    assert_reset()


def test_set_registry_context_manager():
    def run():
        class CustomRegistry(registries.FirstMatchRegistry):
            pass

        with responses.RequestsMock(
            assert_all_requests_are_fired=False, registry=CustomRegistry
        ) as rsps:
            assert type(rsps.get_registry()) == CustomRegistry
            assert type(responses.mock.get_registry()) == registries.FirstMatchRegistry

    run()
    assert_reset()


def test_registry_reset():
    def run():
        class CustomRegistry(registries.FirstMatchRegistry):
            pass

        with responses.RequestsMock(
            assert_all_requests_are_fired=False, registry=CustomRegistry
        ) as rsps:
            rsps.get_registry().reset()
            assert not rsps.registered()

    run()
    assert_reset()


class TestOrderedRegistry:
    def test_invocation_index(self):
        @responses.activate(registry=OrderedRegistry)
        def run():
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/foobar",
                status=666,
            )
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/foobar",
                status=667,
            )
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/foobar",
                status=668,
            )
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/foobar",
                status=669,
            )

            resp = requests.get("http://twitter.com/api/1/foobar")
            assert resp.status_code == 666
            resp = requests.get("http://twitter.com/api/1/foobar")
            assert resp.status_code == 667
            resp = requests.get("http://twitter.com/api/1/foobar")
            assert resp.status_code == 668
            resp = requests.get("http://twitter.com/api/1/foobar")
            assert resp.status_code == 669

        run()
        assert_reset()

    def test_not_match(self):
        @responses.activate(registry=OrderedRegistry)
        def run():
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/foobar",
                json={"msg": "not found"},
                status=667,
            )
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/barfoo",
                json={"msg": "not found"},
                status=404,
            )
            responses.add(
                responses.GET,
                "http://twitter.com/api/1/foobar",
                json={"msg": "OK"},
                status=200,
            )

            resp = requests.get("http://twitter.com/api/1/foobar")
            assert resp.status_code == 667

            with pytest.raises(ConnectionError) as excinfo:
                requests.get("http://twitter.com/api/1/foobar")

            msg = str(excinfo.value)
            assert (
                "- GET http://twitter.com/api/1/barfoo Next 'Response' in the "
                "order doesn't match due to the following reason: URL does not match"
            ) in msg

        run()
        assert_reset()

    def test_empty_registry(self):
        @responses.activate(registry=OrderedRegistry)
        def run():
            with pytest.raises(ConnectionError):
                requests.get("http://twitter.com/api/1/foobar")

        run()
        assert_reset()
