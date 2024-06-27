# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Test the kernels service API."""

import asyncio
import json
import os
from tempfile import TemporaryDirectory

import pytest
import tornado


def expected_http_error(error, expected_code, expected_message=None):
    """Check that the error matches the expected output error."""
    e = error.value
    if isinstance(e, tornado.web.HTTPError):
        if expected_code != e.status_code:
            return False
        if expected_message is not None and expected_message != str(e):
            return False
        return True
    elif any(
        [
            isinstance(e, tornado.httpclient.HTTPClientError),
            isinstance(e, tornado.httpclient.HTTPError),
        ]
    ):
        if expected_code != e.code:
            return False
        if expected_message:
            message = json.loads(e.response.body.decode())["message"]
            if expected_message != message:
                return False
        return True


@pytest.fixture
def build_api_tester(jp_serverapp, labapp, fetch_long):
    return BuildAPITester(labapp, fetch_long)


class BuildAPITester:
    """Wrapper for build REST API requests"""

    url = "lab/api/build"

    def __init__(self, labapp, fetch_long):
        self.labapp = labapp
        self.fetch = fetch_long

    async def _req(self, verb, path, body=None):
        return await self.fetch(self.url + path, method=verb, body=body)

    async def getStatus(self):
        return await self._req("GET", "")

    async def build(self):
        return await self._req("POST", "", json.dumps({}))

    async def clear(self):
        return await self._req("DELETE", "")


@pytest.mark.slow
class TestBuildAPI:
    def tempdir(self):
        td = TemporaryDirectory()
        self.tempdirs.append(td)
        return td.name

    def setUp(self):
        # Any TemporaryDirectory objects appended to this list will be cleaned
        # up at the end of the test run.
        self.tempdirs = []

        # TODO(@echarles) Move the cleanup in the fixture.
        @self.addCleanup
        def cleanup_tempdirs():
            for d in self.tempdirs:
                d.cleanup()

    #    @pytest.mark.timeout(timeout=30)
    #    @pytest.mark.gen_test(timeout=30)
    async def test_get_status(self, build_api_tester):
        """Make sure there are no kernels running at the start"""
        r = await build_api_tester.getStatus()
        res = r.body.decode()
        resp = json.loads(res)
        assert "status" in resp
        assert "message" in resp

    #    @pytest.mark.gen_test(timeout=30)
    # FIXME
    @pytest.mark.skipif(os.name == "nt", reason="Currently failing on windows")
    async def test_build(self, build_api_tester):
        r = await build_api_tester.build()
        assert r.code == 200

    #    @pytest.mark.gen_test(timeout=30)
    # FIXME
    @pytest.mark.skipif(os.name == "nt", reason="Currently failing on windows")
    async def test_clear(self, build_api_tester):
        with pytest.raises(tornado.httpclient.HTTPClientError) as e:
            r = await build_api_tester.clear()
            res = r.body.decode()
        assert expected_http_error(e, 500)

        loop = asyncio.get_event_loop()
        asyncio.ensure_future(build_api_tester.build(), loop=loop)  # noqa RUF006

        while True:
            r = await build_api_tester.getStatus()
            res = r.body.decode()
            resp = json.loads(res)
            if resp["status"] == "building":
                break

        r = await build_api_tester.clear()
        assert r.code == 204
