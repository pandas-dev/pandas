from unittest.mock import (
    MagicMock,
    patch,
)

from web.pandas_web import Preprocessors


def test_home_add_releases():
    # Prepare test data
    releases = [
        {
            "tag_name": "v1.0.0",
            "prerelease": False,
            "published_at": "2021-04-01T00:00:00Z",
            "assets": [{"browser_download_url": "https://example.com/download"}],
        },
        {
            "tag_name": "v2.0.0",
            "prerelease": False,
            "published_at": "2022-01-01T00:00:00Z",
            "assets": [{"browser_download_url": "https://example.com/download"}],
        },
        {
            "tag_name": "v1.5.4",
            "prerelease": False,
            "published_at": "2023-08-01T00:00:00Z",
            "assets": [],
        },
        {
            "tag_name": "v1.5.3",
            "prerelease": False,
            "published_at": "2021-07-01T00:00:00Z",
            "assets": [],
        },
        {
            "tag_name": "v1.5.2",
            "prerelease": False,
            "published_at": "2021-06-01T00:00:00Z",
            "assets": [],
        },
    ]
    github_repo_url = "pandas-dev/pandas"
    context = {
        "main": {
            "github_repo_url": github_repo_url,
            "production_url": "https://example.com/",
        },
        "target_path": "/tmp",
        "releases": [],
    }

    # Mock the requests module to return the test data
    resp = MagicMock()
    resp.status_code = 200
    resp.json.return_value = releases
    with patch("requests.get", return_value=resp):
        # Call the function being tested
        Preprocessors.home_add_releases(context)

        # Assert that releases were correctly added to the context
        assert len(context["releases"]) == 3
        assert context["releases"][0]["name"] == "2.0.0"
        assert context["releases"][1]["name"] == "1.5.4"
        assert context["releases"][2]["name"] == "1.0.0"
