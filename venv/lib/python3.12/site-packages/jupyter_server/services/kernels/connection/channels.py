"""An implementation of a kernel connection."""

from __future__ import annotations

import asyncio
import json
import time
import typing as t
import weakref
from concurrent.futures import Future
from textwrap import dedent

from jupyter_client import protocol_version as client_protocol_version  # type:ignore[attr-defined]
from tornado import gen, web
from tornado.ioloop import IOLoop
from tornado.websocket import WebSocketClosedError
from traitlets import Any, Bool, Dict, Float, Instance, Int, List, Unicode, default

try:
    from jupyter_client.jsonutil import json_default
except ImportError:
    from jupyter_client.jsonutil import date_default as json_default

from jupyter_core.utils import ensure_async

from jupyter_server.transutils import _i18n

from ..websocket import KernelWebsocketHandler
from .abc import KernelWebsocketConnectionABC
from .base import (
    BaseKernelWebsocketConnection,
    deserialize_binary_message,
    deserialize_msg_from_ws_v1,
    serialize_binary_message,
    serialize_msg_to_ws_v1,
)


def _ensure_future(f):
    """Wrap a concurrent future as an asyncio future if there is a running loop."""
    try:
        asyncio.get_running_loop()
        return asyncio.wrap_future(f)
    except RuntimeError:
        return f


class ZMQChannelsWebsocketConnection(BaseKernelWebsocketConnection):
    """A Jupyter Server Websocket Connection"""

    limit_rate = Bool(
        True,
        config=True,
        help=_i18n(
            "Whether to limit the rate of IOPub messages (default: True). "
            "If True, use iopub_msg_rate_limit, iopub_data_rate_limit and/or rate_limit_window "
            "to tune the rate."
        ),
    )

    iopub_msg_rate_limit = Float(
        1000,
        config=True,
        help=_i18n(
            """(msgs/sec)
        Maximum rate at which messages can be sent on iopub before they are
        limited."""
        ),
    )

    iopub_data_rate_limit = Float(
        1000000,
        config=True,
        help=_i18n(
            """(bytes/sec)
        Maximum rate at which stream output can be sent on iopub before they are
        limited."""
        ),
    )

    rate_limit_window = Float(
        3,
        config=True,
        help=_i18n(
            """(sec) Time window used to
        check the message and data rate limits."""
        ),
    )

    websocket_handler = Instance(KernelWebsocketHandler)

    @property
    def write_message(self):
        """Alias to the websocket handler's write_message method."""
        return self.websocket_handler.write_message

    # class-level registry of open sessions
    # allows checking for conflict on session-id,
    # which is used as a zmq identity and must be unique.
    _open_sessions: dict[str, KernelWebsocketHandler] = {}
    _open_sockets: t.MutableSet[ZMQChannelsWebsocketConnection] = weakref.WeakSet()

    _kernel_info_future: Future[t.Any]
    _close_future: Future[t.Any]

    channels = Dict({})
    kernel_info_channel = Any(allow_none=True)

    _kernel_info_future = Instance(klass=Future)  # type:ignore[assignment]

    @default("_kernel_info_future")
    def _default_kernel_info_future(self):
        """The default kernel info future."""
        return Future()

    _close_future = Instance(klass=Future)  # type:ignore[assignment]

    @default("_close_future")
    def _default_close_future(self):
        """The default close future."""
        return Future()

    session_key = Unicode("")

    _iopub_window_msg_count = Int()
    _iopub_window_byte_count = Int()
    _iopub_msgs_exceeded = Bool(False)
    _iopub_data_exceeded = Bool(False)
    # Queue of (time stamp, byte count)
    # Allows you to specify that the byte count should be lowered
    # by a delta amount at some point in the future.
    _iopub_window_byte_queue: List[t.Any] = List([])

    @classmethod
    async def close_all(cls):
        """Tornado does not provide a way to close open sockets, so add one."""
        for connection in list(cls._open_sockets):
            connection.disconnect()
            await _ensure_future(connection._close_future)

    @property
    def subprotocol(self):
        """The sub protocol."""
        try:
            protocol = self.websocket_handler.selected_subprotocol
        except Exception:
            protocol = None
        return protocol

    def create_stream(self):
        """Create a stream."""
        identity = self.session.bsession
        for channel in ("iopub", "shell", "control", "stdin"):
            meth = getattr(self.kernel_manager, "connect_" + channel)
            self.channels[channel] = stream = meth(identity=identity)
            stream.channel = channel

    def nudge(self):
        """Nudge the zmq connections with kernel_info_requests
        Returns a Future that will resolve when we have received
        a shell or control reply and at least one iopub message,
        ensuring that zmq subscriptions are established,
        sockets are fully connected, and kernel is responsive.
        Keeps retrying kernel_info_request until these are both received.
        """
        # Do not nudge busy kernels as kernel info requests sent to shell are
        # queued behind execution requests.
        # nudging in this case would cause a potentially very long wait
        # before connections are opened,
        # plus it is *very* unlikely that a busy kernel will not finish
        # establishing its zmq subscriptions before processing the next request.
        if getattr(self.kernel_manager, "execution_state", None) == "busy":
            self.log.debug("Nudge: not nudging busy kernel %s", self.kernel_id)
            f: Future[t.Any] = Future()
            f.set_result(None)
            return _ensure_future(f)
        # Use a transient shell channel to prevent leaking
        # shell responses to the front-end.
        shell_channel = self.kernel_manager.connect_shell()
        # Use a transient control channel to prevent leaking
        # control responses to the front-end.
        control_channel = self.kernel_manager.connect_control()
        # The IOPub used by the client, whose subscriptions we are verifying.
        iopub_channel = self.channels["iopub"]

        info_future: Future[t.Any] = Future()
        iopub_future: Future[t.Any] = Future()
        both_done = gen.multi([info_future, iopub_future])

        def finish(_=None):
            """Ensure all futures are resolved
            which in turn triggers cleanup
            """
            for f in (info_future, iopub_future):
                if not f.done():
                    f.set_result(None)

        def cleanup(_=None):
            """Common cleanup"""
            loop.remove_timeout(nudge_handle)
            iopub_channel.stop_on_recv()
            if not shell_channel.closed():
                shell_channel.close()
            if not control_channel.closed():
                control_channel.close()

        # trigger cleanup when both message futures are resolved
        both_done.add_done_callback(cleanup)

        def on_shell_reply(msg):
            """Handle nudge shell replies."""
            self.log.debug("Nudge: shell info reply received: %s", self.kernel_id)
            if not info_future.done():
                self.log.debug("Nudge: resolving shell future: %s", self.kernel_id)
                info_future.set_result(None)

        def on_control_reply(msg):
            """Handle nudge control replies."""
            self.log.debug("Nudge: control info reply received: %s", self.kernel_id)
            if not info_future.done():
                self.log.debug("Nudge: resolving control future: %s", self.kernel_id)
                info_future.set_result(None)

        def on_iopub(msg):
            """Handle nudge iopub replies."""
            self.log.debug("Nudge: IOPub received: %s", self.kernel_id)
            if not iopub_future.done():
                iopub_channel.stop_on_recv()
                self.log.debug("Nudge: resolving iopub future: %s", self.kernel_id)
                iopub_future.set_result(None)

        iopub_channel.on_recv(on_iopub)
        shell_channel.on_recv(on_shell_reply)
        control_channel.on_recv(on_control_reply)
        loop = IOLoop.current()

        # Nudge the kernel with kernel info requests until we get an IOPub message
        def nudge(count):
            """Nudge the kernel."""
            count += 1
            # check for stopped kernel
            if self.kernel_id not in self.multi_kernel_manager:
                self.log.debug("Nudge: cancelling on stopped kernel: %s", self.kernel_id)
                finish()
                return

            # check for closed zmq socket
            if shell_channel.closed():
                self.log.debug("Nudge: cancelling on closed zmq socket: %s", self.kernel_id)
                finish()
                return

            # check for closed zmq socket
            if control_channel.closed():
                self.log.debug("Nudge: cancelling on closed zmq socket: %s", self.kernel_id)
                finish()
                return

            if not both_done.done():
                log = self.log.warning if count % 10 == 0 else self.log.debug
                log(f"Nudge: attempt {count} on kernel {self.kernel_id}")
                self.session.send(shell_channel, "kernel_info_request")
                self.session.send(control_channel, "kernel_info_request")
                nonlocal nudge_handle  # type: ignore[misc]
                nudge_handle = loop.call_later(0.5, nudge, count)

        nudge_handle = loop.call_later(0, nudge, count=0)

        # resolve with a timeout if we get no response
        future = gen.with_timeout(loop.time() + self.kernel_info_timeout, both_done)
        # ensure we have no dangling resources or unresolved Futures in case of timeout
        future.add_done_callback(finish)
        return _ensure_future(future)

    async def _register_session(self):
        """Ensure we aren't creating a duplicate session.

        If a previous identical session is still open, close it to avoid collisions.
        This is likely due to a client reconnecting from a lost network connection,
        where the socket on our side has not been cleaned up yet.
        """
        self.session_key = f"{self.kernel_id}:{self.session.session}"
        stale_handler = self._open_sessions.get(self.session_key)
        if stale_handler:
            self.log.warning("Replacing stale connection: %s", self.session_key)
            stale_handler.close()
        if (
            self.kernel_id in self.multi_kernel_manager
        ):  # only update open sessions if kernel is actively managed
            self._open_sessions[self.session_key] = t.cast(
                KernelWebsocketHandler, self.websocket_handler
            )

    async def prepare(self):
        """Prepare a kernel connection."""
        # check session collision:
        await self._register_session()
        # then request kernel info, waiting up to a certain time before giving up.
        # We don't want to wait forever, because browsers don't take it well when
        # servers never respond to websocket connection requests.

        if hasattr(self.kernel_manager, "ready"):
            ready = self.kernel_manager.ready
            if not isinstance(ready, asyncio.Future):
                ready = asyncio.wrap_future(ready)
            try:
                await ready
            except Exception as e:
                self.kernel_manager.execution_state = "dead"
                self.kernel_manager.reason = str(e)
                raise web.HTTPError(500, str(e)) from e

        t0 = time.time()
        while not await ensure_async(self.kernel_manager.is_alive()):
            await asyncio.sleep(0.1)
            if (time.time() - t0) > self.multi_kernel_manager.kernel_info_timeout:
                msg = "Kernel never reached an 'alive' state."
                raise TimeoutError(msg)

        self.session.key = self.kernel_manager.session.key
        future = self.request_kernel_info()

        def give_up():
            """Don't wait forever for the kernel to reply"""
            if future.done():
                return
            self.log.warning("Timeout waiting for kernel_info reply from %s", self.kernel_id)
            future.set_result({})

        loop = IOLoop.current()
        loop.add_timeout(loop.time() + self.kernel_info_timeout, give_up)
        # actually wait for it
        await asyncio.wrap_future(future)

    def connect(self):
        """Handle a connection."""
        self.multi_kernel_manager.notify_connect(self.kernel_id)

        # on new connections, flush the message buffer
        buffer_info = self.multi_kernel_manager.get_buffer(self.kernel_id, self.session_key)
        if buffer_info and buffer_info["session_key"] == self.session_key:
            self.log.info("Restoring connection for %s", self.session_key)
            if self.multi_kernel_manager.ports_changed(self.kernel_id):
                # If the kernel's ports have changed (some restarts trigger this)
                # then reset the channels so nudge() is using the correct iopub channel
                self.create_stream()
            else:
                # The kernel's ports have not changed; use the channels captured in the buffer
                self.channels = buffer_info["channels"]

            connected = self.nudge()

            def replay(value):
                replay_buffer = buffer_info["buffer"]
                if replay_buffer:
                    self.log.info("Replaying %s buffered messages", len(replay_buffer))
                    for channel, msg_list in replay_buffer:
                        stream = self.channels[channel]
                        self.handle_outgoing_message(stream, msg_list)

            connected.add_done_callback(replay)
        else:
            try:
                self.create_stream()
                connected = self.nudge()
            except web.HTTPError as e:
                # Do not log error if the kernel is already shutdown,
                # as it's normal that it's not responding
                try:
                    self.multi_kernel_manager.get_kernel(self.kernel_id)
                    self.log.error("Error opening stream: %s", e)
                except KeyError:
                    pass
                # WebSockets don't respond to traditional error codes so we
                # close the connection.
                for stream in self.channels.values():
                    if not stream.closed():
                        stream.close()
                self.disconnect()
                return None

        self.multi_kernel_manager.add_restart_callback(self.kernel_id, self.on_kernel_restarted)
        self.multi_kernel_manager.add_restart_callback(
            self.kernel_id, self.on_restart_failed, "dead"
        )

        def subscribe(value):
            for stream in self.channels.values():
                stream.on_recv_stream(self.handle_outgoing_message)

        connected.add_done_callback(subscribe)
        ZMQChannelsWebsocketConnection._open_sockets.add(self)
        return connected

    def close(self):
        """Close the connection."""
        return self.disconnect()

    def disconnect(self):
        """Handle a disconnect."""
        self.log.debug("Websocket closed %s", self.session_key)
        # unregister myself as an open session (only if it's really me)
        if self._open_sessions.get(self.session_key) is self.websocket_handler:
            self._open_sessions.pop(self.session_key)

        if self.kernel_id in self.multi_kernel_manager:
            self.multi_kernel_manager.notify_disconnect(self.kernel_id)
            self.multi_kernel_manager.remove_restart_callback(
                self.kernel_id,
                self.on_kernel_restarted,
            )
            self.multi_kernel_manager.remove_restart_callback(
                self.kernel_id,
                self.on_restart_failed,
                "dead",
            )

            # start buffering instead of closing if this was the last connection
            if (
                self.kernel_id in self.multi_kernel_manager._kernel_connections
                and self.multi_kernel_manager._kernel_connections[self.kernel_id] == 0
            ):
                self.multi_kernel_manager.start_buffering(
                    self.kernel_id, self.session_key, self.channels
                )
                ZMQChannelsWebsocketConnection._open_sockets.remove(self)
                self._close_future.set_result(None)
                return

        # This method can be called twice, once by self.kernel_died and once
        # from the WebSocket close event. If the WebSocket connection is
        # closed before the ZMQ streams are setup, they could be None.
        for stream in self.channels.values():
            if stream is not None and not stream.closed():
                stream.on_recv(None)
                stream.close()

        self.channels = {}
        try:
            ZMQChannelsWebsocketConnection._open_sockets.remove(self)
            self._close_future.set_result(None)
        except Exception:
            pass

    def handle_incoming_message(self, incoming_msg: str) -> None:
        """Handle incoming messages from Websocket to ZMQ Sockets."""
        ws_msg = incoming_msg
        if not self.channels:
            # already closed, ignore the message
            self.log.debug("Received message on closed websocket %r", ws_msg)
            return

        if self.subprotocol == "v1.kernel.websocket.jupyter.org":
            channel, msg_list = deserialize_msg_from_ws_v1(ws_msg)
            msg = {
                "header": None,
            }
        else:
            if isinstance(ws_msg, bytes):  # type:ignore[unreachable]
                msg = deserialize_binary_message(ws_msg)  # type:ignore[unreachable]
            else:
                msg = json.loads(ws_msg)
            msg_list = []
            channel = msg.pop("channel", None)

        if channel is None:
            self.log.warning("No channel specified, assuming shell: %s", msg)
            channel = "shell"
        if channel not in self.channels:
            self.log.warning("No such channel: %r", channel)
            return
        am = self.multi_kernel_manager.allowed_message_types
        ignore_msg = False
        if am:
            msg["header"] = self.get_part("header", msg["header"], msg_list)
            assert msg["header"] is not None
            if msg["header"]["msg_type"] not in am:  # type:ignore[unreachable]
                self.log.warning(
                    'Received message of type "%s", which is not allowed. Ignoring.'
                    % msg["header"]["msg_type"]
                )
                ignore_msg = True
        if not ignore_msg:
            stream = self.channels[channel]
            if self.subprotocol == "v1.kernel.websocket.jupyter.org":
                self.session.send_raw(stream, msg_list)
            else:
                self.session.send(stream, msg)

    def handle_outgoing_message(self, stream: str, outgoing_msg: list[t.Any]) -> None:
        """Handle the outgoing messages from ZMQ sockets to Websocket."""
        msg_list = outgoing_msg
        _, fed_msg_list = self.session.feed_identities(msg_list)

        if self.subprotocol == "v1.kernel.websocket.jupyter.org":
            msg = {"header": None, "parent_header": None, "content": None}
        else:
            msg = self.session.deserialize(fed_msg_list)

        if isinstance(stream, str):
            stream = self.channels[stream]

        channel = getattr(stream, "channel", None)
        parts = fed_msg_list[1:]

        self._on_error(channel, msg, parts)

        if self._limit_rate(channel, msg, parts):
            return

        if self.subprotocol == "v1.kernel.websocket.jupyter.org":
            self._on_zmq_reply(stream, parts)
        else:
            self._on_zmq_reply(stream, msg)

    def get_part(self, field, value, msg_list):
        """Get a part of a message."""
        if value is None:
            field2idx = {
                "header": 0,
                "parent_header": 1,
                "content": 3,
            }
            value = self.session.unpack(msg_list[field2idx[field]])
        return value

    def _reserialize_reply(self, msg_or_list, channel=None):
        """Reserialize a reply message using JSON.

        msg_or_list can be an already-deserialized msg dict or the zmq buffer list.
        If it is the zmq list, it will be deserialized with self.session.

        This takes the msg list from the ZMQ socket and serializes the result for the websocket.
        This method should be used by self._on_zmq_reply to build messages that can
        be sent back to the browser.

        """
        if isinstance(msg_or_list, dict):
            # already unpacked
            msg = msg_or_list
        else:
            _, msg_list = self.session.feed_identities(msg_or_list)
            msg = self.session.deserialize(msg_list)
        if channel:
            msg["channel"] = channel
        if msg["buffers"]:
            buf = serialize_binary_message(msg)
            return buf
        else:
            return json.dumps(msg, default=json_default)

    def _on_zmq_reply(self, stream, msg_list):
        """Handle a zmq reply."""
        # Sometimes this gets triggered when the on_close method is scheduled in the
        # eventloop but hasn't been called.
        if stream.closed():
            self.log.warning("zmq message arrived on closed channel")
            self.disconnect()
            return
        channel = getattr(stream, "channel", None)
        if self.subprotocol == "v1.kernel.websocket.jupyter.org":
            bin_msg = serialize_msg_to_ws_v1(msg_list, channel)
            self.write_message(bin_msg, binary=True)
        else:
            try:
                msg = self._reserialize_reply(msg_list, channel=channel)
            except Exception:
                self.log.critical("Malformed message: %r" % msg_list, exc_info=True)
            else:
                try:
                    self.write_message(msg, binary=isinstance(msg, bytes))
                except WebSocketClosedError as e:
                    self.log.warning(str(e))

    def request_kernel_info(self):
        """send a request for kernel_info"""
        try:
            # check for previous request
            future = self.kernel_manager._kernel_info_future
        except AttributeError:
            self.log.debug("Requesting kernel info from %s", self.kernel_id)
            # Create a kernel_info channel to query the kernel protocol version.
            # This channel will be closed after the kernel_info reply is received.
            if self.kernel_info_channel is None:
                self.kernel_info_channel = self.multi_kernel_manager.connect_shell(self.kernel_id)
            assert self.kernel_info_channel is not None
            self.kernel_info_channel.on_recv(self._handle_kernel_info_reply)
            self.session.send(self.kernel_info_channel, "kernel_info_request")
            # store the future on the kernel, so only one request is sent
            self.kernel_manager._kernel_info_future = self._kernel_info_future
        else:
            if not future.done():
                self.log.debug("Waiting for pending kernel_info request")
            future.add_done_callback(lambda f: self._finish_kernel_info(f.result()))
        return _ensure_future(self._kernel_info_future)

    def _handle_kernel_info_reply(self, msg):
        """process the kernel_info_reply

        enabling msg spec adaptation, if necessary
        """
        idents, msg = self.session.feed_identities(msg)
        try:
            msg = self.session.deserialize(msg)
        except BaseException:
            self.log.error("Bad kernel_info reply", exc_info=True)
            self._kernel_info_future.set_result({})
            return
        else:
            info = msg["content"]
            self.log.debug("Received kernel info: %s", info)
            if msg["msg_type"] != "kernel_info_reply" or "protocol_version" not in info:
                self.log.error("Kernel info request failed, assuming current %s", info)
                info = {}
            self._finish_kernel_info(info)

        # close the kernel_info channel, we don't need it anymore
        if self.kernel_info_channel:
            self.kernel_info_channel.close()
        self.kernel_info_channel = None

    def _finish_kernel_info(self, info):
        """Finish handling kernel_info reply

        Set up protocol adaptation, if needed,
        and signal that connection can continue.
        """
        protocol_version = info.get("protocol_version", client_protocol_version)
        if protocol_version != client_protocol_version:
            self.session.adapt_version = int(protocol_version.split(".")[0])
            self.log.info(
                f"Adapting from protocol version {protocol_version} (kernel {self.kernel_id}) to {client_protocol_version} (client)."
            )
        if not self._kernel_info_future.done():
            self._kernel_info_future.set_result(info)

    def write_stderr(self, error_message, parent_header):
        """Write a message to stderr."""
        self.log.warning(error_message)
        err_msg = self.session.msg(
            "stream",
            content={"text": error_message + "\n", "name": "stderr"},
            parent=parent_header,
        )
        if self.subprotocol == "v1.kernel.websocket.jupyter.org":
            bin_msg = serialize_msg_to_ws_v1(err_msg, "iopub", self.session.pack)
            self.write_message(bin_msg, binary=True)
        else:
            err_msg["channel"] = "iopub"
            self.write_message(json.dumps(err_msg, default=json_default))

    def _limit_rate(self, channel, msg, msg_list):
        """Limit the message rate on a channel."""
        if not (self.limit_rate and channel == "iopub"):
            return False

        msg["header"] = self.get_part("header", msg["header"], msg_list)

        msg_type = msg["header"]["msg_type"]
        if msg_type == "status":
            msg["content"] = self.get_part("content", msg["content"], msg_list)
            if msg["content"].get("execution_state") == "idle":
                # reset rate limit counter on status=idle,
                # to avoid 'Run All' hitting limits prematurely.
                self._iopub_window_byte_queue = []
                self._iopub_window_msg_count = 0
                self._iopub_window_byte_count = 0
                self._iopub_msgs_exceeded = False
                self._iopub_data_exceeded = False

        if msg_type not in {"status", "comm_open", "execute_input"}:
            # Remove the counts queued for removal.
            now = IOLoop.current().time()
            while len(self._iopub_window_byte_queue) > 0:
                queued = self._iopub_window_byte_queue[0]
                if now >= queued[0]:
                    self._iopub_window_byte_count -= queued[1]
                    self._iopub_window_msg_count -= 1
                    del self._iopub_window_byte_queue[0]
                else:
                    # This part of the queue hasn't be reached yet, so we can
                    # abort the loop.
                    break

            # Increment the bytes and message count
            self._iopub_window_msg_count += 1
            byte_count = sum(len(x) for x in msg_list) if msg_type == "stream" else 0
            self._iopub_window_byte_count += byte_count

            # Queue a removal of the byte and message count for a time in the
            # future, when we are no longer interested in it.
            self._iopub_window_byte_queue.append((now + self.rate_limit_window, byte_count))

            # Check the limits, set the limit flags, and reset the
            # message and data counts.
            msg_rate = float(self._iopub_window_msg_count) / self.rate_limit_window
            data_rate = float(self._iopub_window_byte_count) / self.rate_limit_window

            # Check the msg rate
            if self.iopub_msg_rate_limit > 0 and msg_rate > self.iopub_msg_rate_limit:
                if not self._iopub_msgs_exceeded:
                    self._iopub_msgs_exceeded = True
                    msg["parent_header"] = self.get_part(
                        "parent_header", msg["parent_header"], msg_list
                    )
                    self.write_stderr(
                        dedent(
                            f"""\
                    IOPub message rate exceeded.
                    The Jupyter server will temporarily stop sending output
                    to the client in order to avoid crashing it.
                    To change this limit, set the config variable
                    `--ServerApp.iopub_msg_rate_limit`.

                    Current values:
                    ServerApp.iopub_msg_rate_limit={self.iopub_msg_rate_limit} (msgs/sec)
                    ServerApp.rate_limit_window={self.rate_limit_window} (secs)
                    """
                        ),
                        msg["parent_header"],
                    )
            # resume once we've got some headroom below the limit
            elif self._iopub_msgs_exceeded and msg_rate < (0.8 * self.iopub_msg_rate_limit):
                self._iopub_msgs_exceeded = False
                if not self._iopub_data_exceeded:
                    self.log.warning("iopub messages resumed")

            # Check the data rate
            if self.iopub_data_rate_limit > 0 and data_rate > self.iopub_data_rate_limit:
                if not self._iopub_data_exceeded:
                    self._iopub_data_exceeded = True
                    msg["parent_header"] = self.get_part(
                        "parent_header", msg["parent_header"], msg_list
                    )
                    self.write_stderr(
                        dedent(
                            f"""\
                    IOPub data rate exceeded.
                    The Jupyter server will temporarily stop sending output
                    to the client in order to avoid crashing it.
                    To change this limit, set the config variable
                    `--ServerApp.iopub_data_rate_limit`.

                    Current values:
                    ServerApp.iopub_data_rate_limit={self.iopub_data_rate_limit} (bytes/sec)
                    ServerApp.rate_limit_window={self.rate_limit_window} (secs)
                    """
                        ),
                        msg["parent_header"],
                    )
            # resume once we've got some headroom below the limit
            elif self._iopub_data_exceeded and data_rate < (0.8 * self.iopub_data_rate_limit):
                self._iopub_data_exceeded = False
                if not self._iopub_msgs_exceeded:
                    self.log.warning("iopub messages resumed")

            # If either of the limit flags are set, do not send the message.
            if self._iopub_msgs_exceeded or self._iopub_data_exceeded:
                # we didn't send it, remove the current message from the calculus
                self._iopub_window_msg_count -= 1
                self._iopub_window_byte_count -= byte_count
                self._iopub_window_byte_queue.pop(-1)
                return True

            return False

    def _send_status_message(self, status):
        """Send a status message."""
        iopub = self.channels.get("iopub", None)
        if iopub and not iopub.closed():
            # flush IOPub before sending a restarting/dead status message
            # ensures proper ordering on the IOPub channel
            # that all messages from the stopped kernel have been delivered
            iopub.flush()
        msg = self.session.msg("status", {"execution_state": status})
        if self.subprotocol == "v1.kernel.websocket.jupyter.org":
            bin_msg = serialize_msg_to_ws_v1(msg, "iopub", self.session.pack)
            self.write_message(bin_msg, binary=True)
        else:
            msg["channel"] = "iopub"
            self.write_message(json.dumps(msg, default=json_default))

    def on_kernel_restarted(self):
        """Handle a kernel restart."""
        self.log.warning("kernel %s restarted", self.kernel_id)
        self._send_status_message("restarting")

    def on_restart_failed(self):
        """Handle a kernel restart failure."""
        self.log.error("kernel %s restarted failed!", self.kernel_id)
        self._send_status_message("dead")

    def _on_error(self, channel, msg, msg_list):
        """Handle an error message."""
        if self.multi_kernel_manager.allow_tracebacks:
            return

        if channel == "iopub":
            msg["header"] = self.get_part("header", msg["header"], msg_list)
            if msg["header"]["msg_type"] == "error":
                msg["content"] = self.get_part("content", msg["content"], msg_list)
                msg["content"]["ename"] = "ExecutionError"
                msg["content"]["evalue"] = "Execution error"
                msg["content"]["traceback"] = [self.kernel_manager.traceback_replacement_message]
                if self.subprotocol == "v1.kernel.websocket.jupyter.org":
                    msg_list[3] = self.session.pack(msg["content"])


KernelWebsocketConnectionABC.register(ZMQChannelsWebsocketConnection)
