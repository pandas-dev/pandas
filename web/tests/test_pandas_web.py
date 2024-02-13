from unittest.mock import (  # noqa: TID251
    mock_open,
    patch,
)

import pytest
import requests

from web.pandas_web import Preprocessors


class MockResponse:
    def __init__(self, status_code: int, response: dict) -> None:
        self.status_code = status_code
        self._resp = response

    def json(self):
        return self._resp

    @staticmethod
    def raise_for_status() -> None:
        return


@pytest.fixture
def context() -> dict:
    return {
        "main": {"github_repo_url": "pandas-dev/pandas"},
        "target_path": "test_target_path",
    }


@pytest.fixture
def mock_response(monkeypatch, request) -> None:
    def mocked_resp(*args, **kwargs):
        status_code, response = request.param
        return MockResponse(status_code, response)

    monkeypatch.setattr(requests, "get", mocked_resp)


_releases_list = [
    {
        "prerelease": False,
        "published_at": "2024-01-19T03:34:05Z",
        "tag_name": "v1.5.6",
        "assets": None,
    },
    {
        "prerelease": False,
        "published_at": "2023-11-10T19:07:37Z",
        "tag_name": "v2.1.3",
        "assets": None,
    },
    {
        "prerelease": False,
        "published_at": "2023-08-30T13:24:32Z",
        "tag_name": "v2.1.0",
        "assets": None,
    },
    {
        "prerelease": False,
        "published_at": "2023-04-30T13:24:32Z",
        "tag_name": "v2.0.0",
        "assets": None,
    },
    {
        "prerelease": True,
        "published_at": "2023-01-19T03:34:05Z",
        "tag_name": "v1.5.3xd",
        "assets": None,
    },
    {
        "prerelease": False,
        "published_at": "2027-01-19T03:34:05Z",
        "tag_name": "v10.0.1",
        "assets": None,
    },
]


@pytest.mark.parametrize("mock_response", [(200, _releases_list)], indirect=True)
def test_web_preprocessor_creates_releases(mock_response, context) -> None:
    m = mock_open()
    with patch("builtins.open", m):
        context = Preprocessors.home_add_releases(context)
        release_versions = [release["name"] for release in context["releases"]]
        assert release_versions == ["10.0.1", "2.1.3", "2.0.0", "1.5.6"]
