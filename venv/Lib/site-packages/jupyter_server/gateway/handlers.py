"""Gateway API handlers."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import asyncio
import logging
import mimetypes
import os
import random
import warnings
from typing import Any, Optional, cast

from jupyter_client.session import Session
from tornado import web
from tornado.concurrent import Future
from tornado.escape import json_decode, url_escape, utf8
from tornado.httpclient import HTTPRequest
from tornado.ioloop import IOLoop, PeriodicCallback
from tornado.websocket import WebSocketHandler, websocket_connect
from traitlets.config.configurable import LoggingConfigurable

from ..base.handlers import APIHandler, JupyterHandler
from ..utils import url_path_join
from .gateway_client import GatewayClient

warnings.warn(
    "The jupyter_server.gateway.handlers module is deprecated and will not be supported in Jupyter Server 3.0",
    DeprecationWarning,
    stacklevel=2,
)


# Keepalive ping interval (default: 30 seconds)
GATEWAY_WS_PING_INTERVAL_SECS = int(os.getenv("GATEWAY_WS_PING_INTERVAL_SECS", "30"))


class WebSocketChannelsHandler(WebSocketHandler, JupyterHandler):
    """Gateway web socket channels handler."""

    session = None
    gateway = None
    kernel_id = None
    ping_callback = None

    def check_origin(self, origin=None):
        """Check origin for the socket."""
        return JupyterHandler.check_origin(self, origin)

    def set_default_headers(self):
        """Undo the set_default_headers in JupyterHandler which doesn't make sense for websockets"""

    def get_compression_options(self):
        """Get the compression options for the socket."""
        # use deflate compress websocket
        return {}

    def authenticate(self):
        """Run before finishing the GET request

        Extend this method to add logic that should fire before
        the websocket finishes completing.
        """
        # authenticate the request before opening the websocket
        if self.current_user is None:
            self.log.warning("Couldn't authenticate WebSocket connection")
            raise web.HTTPError(403)

        if self.get_argument("session_id", None):
            assert self.session is not None
            self.session.session = self.get_argument("session_id")  # type:ignore[unreachable]
        else:
            self.log.warning("No session ID specified")

    def initialize(self):
        """Initialize the socket."""
        self.log.debug("Initializing websocket connection %s", self.request.path)
        self.session = Session(config=self.config)
        self.gateway = GatewayWebSocketClient(gateway_url=GatewayClient.instance().url)

    async def get(self, kernel_id, *args, **kwargs):
        """Get the socket."""
        self.authenticate()
        self.kernel_id = kernel_id
        kwargs["kernel_id"] = kernel_id
        await super().get(*args, **kwargs)

    def send_ping(self):
        """Send a ping to the socket."""
        if self.ws_connection is None and self.ping_callback is not None:
            self.ping_callback.stop()  # type:ignore[unreachable]
            return

        self.ping(b"")

    def open(self, kernel_id, *args, **kwargs):
        """Handle web socket connection open to notebook server and delegate to gateway web socket handler"""
        self.ping_callback = PeriodicCallback(self.send_ping, GATEWAY_WS_PING_INTERVAL_SECS * 1000)
        self.ping_callback.start()

        assert self.gateway is not None
        self.gateway.on_open(
            kernel_id=kernel_id,
            message_callback=self.write_message,
            compression_options=self.get_compression_options(),
        )

    def on_message(self, message):
        """Forward message to gateway web socket handler."""
        assert self.gateway is not None
        self.gateway.on_message(message)

    def write_message(self, message, binary=False):
        """Send message back to notebook client.  This is called via callback from self.gateway._read_messages."""
        if self.ws_connection:  # prevent WebSocketClosedError
            if isinstance(message, bytes):
                binary = True
            super().write_message(message, binary=binary)
        elif self.log.isEnabledFor(logging.DEBUG):
            msg_summary = WebSocketChannelsHandler._get_message_summary(json_decode(utf8(message)))
            self.log.debug(
                f"Notebook client closed websocket connection - message dropped: {msg_summary}"
            )

    def on_close(self):
        """Handle a closing socket."""
        self.log.debug("Closing websocket connection %s", self.request.path)
        assert self.gateway is not None
        self.gateway.on_close()
        super().on_close()

    @staticmethod
    def _get_message_summary(message):
        """Get a summary of a message."""
        summary = []
        message_type = message["msg_type"]
        summary.append(f"type: {message_type}")

        if message_type == "status":
            summary.append(", state: {}".format(message["content"]["execution_state"]))
        elif message_type == "error":
            summary.append(
                ", {}:{}:{}".format(
                    message["content"]["ename"],
                    message["content"]["evalue"],
                    message["content"]["traceback"],
                )
            )
        else:
            summary.append(", ...")  # don't display potentially sensitive data

        return "".join(summary)


class GatewayWebSocketClient(LoggingConfigurable):
    """Proxy web socket connection to a kernel/enterprise gateway."""

    def __init__(self, **kwargs):
        """Initialize the gateway web socket client."""
        super().__init__()
        self.kernel_id = None
        self.ws = None
        self.ws_future: Future[Any] = Future()
        self.disconnected = False
        self.retry = 0

    async def _connect(self, kernel_id, message_callback):
        """Connect to the socket."""
        # websocket is initialized before connection
        self.ws = None
        self.kernel_id = kernel_id
        client = GatewayClient.instance()
        assert client.ws_url is not None

        ws_url = url_path_join(
            client.ws_url,
            client.kernels_endpoint,
            url_escape(kernel_id),
            "channels",
        )
        self.log.info(f"Connecting to {ws_url}")
        kwargs: dict[str, Any] = {}
        kwargs = client.load_connection_args(**kwargs)

        request = HTTPRequest(ws_url, **kwargs)
        self.ws_future = cast("Future[Any]", websocket_connect(request))
        self.ws_future.add_done_callback(self._connection_done)

        loop = IOLoop.current()
        loop.add_future(self.ws_future, lambda future: self._read_messages(message_callback))

    def _connection_done(self, fut):
        """Handle a finished connection."""
        if (
            not self.disconnected and fut.exception() is None
        ):  # prevent concurrent.futures._base.CancelledError
            self.ws = fut.result()
            self.retry = 0
            self.log.debug(f"Connection is ready: ws: {self.ws}")
        else:
            self.log.warning(
                "Websocket connection has been closed via client disconnect or due to error.  "
                f"Kernel with ID '{self.kernel_id}' may not be terminated on GatewayClient: {GatewayClient.instance().url}"
            )

    def _disconnect(self):
        """Handle a disconnect."""
        self.disconnected = True
        if self.ws is not None:
            # Close connection
            self.ws.close()
        elif not self.ws_future.done():
            # Cancel pending connection.  Since future.cancel() is a noop on tornado, we'll track cancellation locally
            self.ws_future.cancel()
            self.log.debug(f"_disconnect: future cancelled, disconnected: {self.disconnected}")

    async def _read_messages(self, callback):
        """Read messages from gateway server."""
        while self.ws is not None:
            message = None
            if not self.disconnected:
                try:
                    message = await self.ws.read_message()
                except Exception as e:
                    self.log.error(
                        f"Exception reading message from websocket: {e}"
                    )  # , exc_info=True)
                if message is None:
                    if not self.disconnected:
                        self.log.warning(f"Lost connection to Gateway: {self.kernel_id}")
                    break
                callback(
                    message
                )  # pass back to notebook client (see self.on_open and WebSocketChannelsHandler.open)
            else:  # ws cancelled - stop reading
                break

        # NOTE(esevan): if websocket is not disconnected by client, try to reconnect.
        if not self.disconnected and self.retry < GatewayClient.instance().gateway_retry_max:
            jitter = random.randint(10, 100) * 0.01  # noqa: S311
            retry_interval = (
                min(
                    GatewayClient.instance().gateway_retry_interval * (2**self.retry),
                    GatewayClient.instance().gateway_retry_interval_max,
                )
                + jitter
            )
            self.retry += 1
            self.log.info(
                "Attempting to re-establish the connection to Gateway in %s secs (%s/%s): %s",
                retry_interval,
                self.retry,
                GatewayClient.instance().gateway_retry_max,
                self.kernel_id,
            )
            await asyncio.sleep(retry_interval)
            loop = IOLoop.current()
            loop.spawn_callback(self._connect, self.kernel_id, callback)

    def on_open(self, kernel_id, message_callback, **kwargs):
        """Web socket connection open against gateway server."""
        loop = IOLoop.current()
        loop.spawn_callback(self._connect, kernel_id, message_callback)

    def on_message(self, message):
        """Send message to gateway server."""
        if self.ws is None:
            loop = IOLoop.current()
            loop.add_future(self.ws_future, lambda future: self._write_message(message))
        else:
            self._write_message(message)

    def _write_message(self, message):
        """Send message to gateway server."""
        try:
            if not self.disconnected and self.ws is not None:
                self.ws.write_message(message)
        except Exception as e:
            self.log.error(f"Exception writing message to websocket: {e}")  # , exc_info=True)

    def on_close(self):
        """Web socket closed event."""
        self._disconnect()


class GatewayResourceHandler(APIHandler):
    """Retrieves resources for specific kernelspec definitions from kernel/enterprise gateway."""

    @web.authenticated
    async def get(self, kernel_name, path, include_body=True):
        """Get a gateway resource by name and path."""
        mimetype: Optional[str] = None
        ksm = self.kernel_spec_manager
        kernel_spec_res = await ksm.get_kernel_spec_resource(  # type:ignore[attr-defined]
            kernel_name, path
        )
        if kernel_spec_res is None:
            self.log.warning(
                f"Kernelspec resource '{path}' for '{kernel_name}' not found.  Gateway may not support"
                " resource serving."
            )
        else:
            mimetype = mimetypes.guess_type(path)[0] or "text/plain"
        self.finish(kernel_spec_res, set_content_type=mimetype)


from ..services.kernels.handlers import _kernel_id_regex
from ..services.kernelspecs.handlers import kernel_name_regex

default_handlers = [
    (r"/api/kernels/%s/channels" % _kernel_id_regex, WebSocketChannelsHandler),
    (r"/kernelspecs/%s/(?P<path>.*)" % kernel_name_regex, GatewayResourceHandler),
]
