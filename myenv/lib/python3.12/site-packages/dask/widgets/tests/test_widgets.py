from __future__ import annotations

import os.path

import pytest

jinja2 = pytest.importorskip("jinja2")

from dask.utils import format_bytes
from dask.widgets import FILTERS, TEMPLATE_PATHS, get_environment, get_template


@pytest.fixture(autouse=True)
def setup_testing():
    TEMPLATE_PATHS.append(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
    )
    FILTERS["custom_filter"] = lambda x: "baz"


def test_widgets():
    template = get_template("example.html.j2")
    assert isinstance(template, jinja2.Template)
    rendered = template.render(foo="bar")
    assert "Hello bar" in rendered


def test_environment():
    environment = get_environment()
    assert isinstance(environment, jinja2.Environment)


def test_unknown_template():
    with pytest.raises(jinja2.TemplateNotFound) as e:
        get_template("does_not_exist.html.j2")

        # The error should contain all the registered template directories to help the user
        # understand where jinja2 is looking. Including the one we registered in the fixture.
        assert os.path.dirname(os.path.abspath(__file__)) in str(e)


def test_filters():
    template = get_template("bytes.html.j2")
    assert format_bytes in FILTERS.values()
    assert format_bytes(2e9) in template.render(foo=2e9)

    template = get_template("custom_filter.html.j2")
    assert "baz" in template.render(foo=None)
