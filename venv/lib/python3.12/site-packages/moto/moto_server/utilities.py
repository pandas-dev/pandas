import gzip
import json
from typing import Any, Dict
from urllib.parse import urlencode

from flask import current_app, request
from flask.testing import FlaskClient
from werkzeug.routing import BaseConverter


class RegexConverter(BaseConverter):
    # http://werkzeug.pocoo.org/docs/routing/#custom-converters

    part_isolating = False

    def __init__(self, url_map: Any, *items: Any):
        super().__init__(url_map)
        self.regex = items[0]


class AWSTestHelper(FlaskClient):
    def action_data(self, action_name: str, **kwargs: Any) -> str:
        """
        Method calls resource with action_name and returns data of response.
        """
        opts = {"Action": action_name}
        opts.update(kwargs)
        res = self.get(
            f"/?{urlencode(opts)}",
            headers={"Host": f"{self.application.service}.us-east-1.amazonaws.com"},  # type: ignore[attr-defined]
        )
        return res.data.decode("utf-8")

    def action_json(self, action_name: str, **kwargs: Any) -> Dict[str, Any]:
        """
        Method calls resource with action_name and returns object obtained via
        deserialization of output.
        """
        return json.loads(self.action_data(action_name, **kwargs))


def decompress_request_body() -> None:
    """Intended to be used as a `before_request` handler in the Moto Server Flask app."""
    if getattr(current_app, "service") == "s3":
        return
    if getattr(request, "content_encoding") == "gzip":
        request.stream = gzip.GzipFile(fileobj=request.stream)  # type: ignore[assignment]
        request.get_data(parse_form_data=True)
