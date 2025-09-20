from __future__ import annotations

import datetime
import html
import os.path

from jinja2 import Environment, FileSystemLoader, Template
from jinja2.exceptions import TemplateNotFound

from dask.utils import format_bytes, format_time, format_time_ago, key_split, typename

FILTERS = {
    "datetime_from_timestamp": datetime.datetime.fromtimestamp,
    "format_bytes": format_bytes,
    "format_time": format_time,
    "format_time_ago": format_time_ago,
    "html_escape": html.escape,
    "key_split": key_split,
    "type": type,
    "typename": typename,
}

TEMPLATE_PATHS = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")]


def get_environment() -> Environment:
    loader = FileSystemLoader(TEMPLATE_PATHS)
    environment = Environment(loader=loader)
    environment.filters.update(FILTERS)

    return environment


def get_template(name: str) -> Template:
    try:
        return get_environment().get_template(name)
    except TemplateNotFound as e:
        raise TemplateNotFound(
            f"Unable to find {name} in dask.widgets.TEMPLATE_PATHS {TEMPLATE_PATHS}"
        ) from e
