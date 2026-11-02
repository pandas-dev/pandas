from unittest.mock import (  # noqa: TID251
    mock_open,
    patch,
)

import pytest
import requests

from web.pandas_web import (
    Preprocessors,
    main,
)


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


def _create_minimal_source(path) -> None:
    """Create the smallest source directory ``main`` can build offline."""
    path.mkdir(parents=True)
    (path / "config.yml").write_text(
        "main:\n"
        "  templates_path: templates\n"
        "  ignore: []\n"
        "  context_preprocessors: []\n",
        encoding="utf-8",
    )
    (path / "versions.json").write_text("{}", encoding="utf-8")


@pytest.mark.parametrize("target", ["site", "site/build", "."])
def test_web_main_overlapping_source_and_target(tmp_path, target) -> None:
    # GH#70082: an overlapping target is removed before rendering, which
    # would delete the source files.
    source = tmp_path / "site"
    source.mkdir()
    (source / "config.yml").touch()
    with pytest.raises(
        ValueError, match="Target path must not equal, contain, or be inside"
    ):
        main(source, tmp_path / target)
    assert (source / "config.yml").exists()


def test_web_main_disjoint_source_and_target(tmp_path) -> None:
    # Non-overlapping source and target must still build normally.
    source = tmp_path / "site"
    _create_minimal_source(source)
    target = tmp_path / "build"

    main(source, target)

    assert (source / "config.yml").exists()
    assert (source / "versions.json").exists()
    assert (target / "config.yml").exists()
    assert (target / "versions.json").exists()
