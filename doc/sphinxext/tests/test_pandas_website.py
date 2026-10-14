import datetime
import pathlib
import sys
from unittest.mock import (  # noqa: TID251
    mock_open,
    patch,
)

import pytest
import requests

sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
from pandas_website import Preprocessors


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


def test_web_preprocessor_blog_posts_from_front_matter(tmp_path: pathlib.Path) -> None:
    posts_path = tmp_path / "community" / "blog"
    posts_path.mkdir(parents=True)
    (posts_path / "my-post.md").write_text(
        "---\n"
        "title: 'My post: a test'\n"
        "date: 2026-01-21\n"
        "---\n"
        "\n"
        "# My post: a test\n"
        "\n"
        "Some *text* in the post.\n",
        encoding="utf-8",
    )
    context = {
        "source_path": tmp_path,
        "blog": {
            "posts_path": "community/blog",
            "author": "pandas team",
            "feed_name": "pandas blog",
            "feed": [],
            "num_posts": 10,
        },
    }
    context = Preprocessors.blog_add_posts(context)
    (post,) = context["blog"]["posts"]
    assert post["title"] == "My post: a test"
    assert post["published"] == datetime.datetime(2026, 1, 21)
    assert post["link"] == "community/blog/my-post.html"
    assert post["summary"].strip() == "Some text in the post."


def test_legacy_heading_ids() -> None:
    # ids of the headings in the pages previously generated with Python-Markdown,
    # kept so links to sections of the website pages keep working
    from pandas_website import _legacy_heading_id

    used_ids: set[str] = set()
    titles = [
        "Sponsors",
        "Highlights of pandas 3.0",
        "1. Dedicated string data type by default",
        "PDEP-7: Consistent copy/view semantics in pandas with Copy-on-Write",
        "Sponsors",
        "Sponsors",
    ]
    assert [_legacy_heading_id(title, used_ids) for title in titles] == [
        "sponsors",
        "highlights-of-pandas-30",
        "1-dedicated-string-data-type-by-default",
        "pdep-7-consistent-copyview-semantics-in-pandas-with-copy-on-write",
        "sponsors_1",
        "sponsors_2",
    ]
