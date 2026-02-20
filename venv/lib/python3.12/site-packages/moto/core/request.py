from __future__ import annotations

from typing import TYPE_CHECKING
from urllib.parse import urlparse

from werkzeug.wrappers import Request

from moto.core.utils import gzip_decompress
from moto.settings import MAX_FORM_MEMORY_SIZE
from moto.utilities.constants import APPLICATION_JSON, JSON_TYPES

if TYPE_CHECKING:
    from botocore.awsrequest import AWSPreparedRequest

    from moto.core.model import ServiceModel


def normalize_request(request: AWSPreparedRequest | Request) -> Request:
    if isinstance(request, Request):
        return request
    body = request.body
    # Request.from_values() does not automatically handle gzip-encoded bodies,
    # like the full WSGI server would, so we need to do it manually.
    if request.headers.get("Content-Encoding") == "gzip":
        body = gzip_decompress(body)  # type: ignore[arg-type]
    parsed_url = urlparse(request.url)
    Request.max_form_memory_size = MAX_FORM_MEMORY_SIZE
    normalized_request = Request.from_values(
        method=request.method,
        base_url=f"{parsed_url.scheme}://{parsed_url.netloc}",
        path=parsed_url.path,
        query_string=parsed_url.query,
        data=body,
        headers=[(k, v) for k, v in request.headers.items()],
    )
    return normalized_request


def determine_request_protocol(
    service_model: ServiceModel, content_type: str | None = None
) -> str:
    protocol = str(service_model.protocol)
    supported_protocols = service_model.metadata.get("protocols", [protocol])
    content_type = content_type if content_type is not None else ""
    if content_type in JSON_TYPES:
        protocol = "rest-json" if content_type == APPLICATION_JSON else "json"
    elif content_type.startswith("application/x-www-form-urlencoded"):
        protocol = "ec2" if "ec2" in supported_protocols else "query"
    if protocol not in supported_protocols:
        raise NotImplementedError(
            f"Unsupported protocol [{protocol}] for service {service_model.service_name}"
        )
    return protocol
