from __future__ import annotations

from aiohttp import hdrs, web
from aiohttp.typedefs import Handler

from ..exposition import _bake_output
from ..registry import Collector, REGISTRY


def make_aiohttp_handler(
    registry: Collector = REGISTRY,
    disable_compression: bool = False,
) -> Handler:
    """Create a aiohttp handler which serves the metrics from a registry."""

    async def prometheus_handler(request: web.Request) -> web.Response:
        # Prepare parameters
        params = {key: request.query.getall(key) for key in request.query.keys()}
        accept_header = ",".join(request.headers.getall(hdrs.ACCEPT, []))
        accept_encoding_header = ""
        # Bake output
        status, headers, output = _bake_output(
            registry,
            accept_header,
            accept_encoding_header,
            params,
            # use AIOHTTP's compression
            disable_compression=True,
        )
        response = web.Response(
            status=int(status.split(" ")[0]),
            headers=headers,
            body=output,
        )
        if not disable_compression:
            response.enable_compression()
        return response

    return prometheus_handler
