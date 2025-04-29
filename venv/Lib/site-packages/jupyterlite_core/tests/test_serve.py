"""Test that various serving options work"""

import json
import subprocess
import time

import pytest
from tornado import httpclient

from .conftest import CI, DARWIN, LINUX, PYPY

if CI and DARWIN:  # pragma: no cover
    pytest.skip("skipping flaky MacOS tests", allow_module_level=True)

if CI and LINUX and PYPY:  # pragma: no cover
    pytest.skip("skipping flaky Linux/PyPy tests", allow_module_level=True)


@pytest.mark.parametrize("base_url", [None, "/@foo/"])
def test_serve(an_empty_lite_dir, script_runner, base_url, an_unused_port):  # pragma: no cover
    """verify that serving kinda works"""
    args = ["jupyter", "lite", "serve", "--port", f"{an_unused_port}"]

    http_headers = {"x-foo": "bar"}
    extra_http_headers = {"x-baz": "boo"}
    all_headers = {**http_headers, **extra_http_headers}

    config = {
        "LiteBuildConfig": {
            "http_headers": http_headers,
            "extra_http_headers": extra_http_headers,
        }
    }

    config_path = an_empty_lite_dir / "jupyter_lite_config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    if base_url:
        args += ["--base-url", base_url]
    else:
        base_url = "/"

    url = f"http://127.0.0.1:{an_unused_port}{base_url}"

    server = subprocess.Popen(args, cwd=str(an_empty_lite_dir))  # noqa: S603
    time.sleep(2)

    app_urls = [""]
    for app in ["lab", "tree", "repl"]:
        app_urls += [
            f"{app}/",
            f"{app}/index.html",
        ]

    maybe_errors = [
        _fetch_without_errors(f"{url}{frag}", expect_headers=all_headers) for frag in app_urls
    ]

    errors = [e for e in maybe_errors if e]

    try:
        assert not errors
    finally:
        _fetch_without_errors(f"{url}shutdown")
        server.wait(timeout=10)


def _fetch_without_errors(url, retries=10, expect_headers=None):  # pragma: no cover
    retries = 10
    errors = []

    while retries:
        retries -= 1
        response = None
        try:
            client = httpclient.HTTPClient()
            response = client.fetch(url)
            assert b"jupyter-config-data" in response.body
            # it worked, eventually: clear errors
            errors = []
            break
        except Exception as error:  # pragma: no cover
            print(f"{error}: {retries} retries left...")
            time.sleep(0.5)
            errors = [error]

    if response and expect_headers:
        errors = []
        for header, value in expect_headers.items():
            try:
                assert response.headers[header] == value
            except Exception as err:  # pragma: no cover
                errors += [err]

    return errors
