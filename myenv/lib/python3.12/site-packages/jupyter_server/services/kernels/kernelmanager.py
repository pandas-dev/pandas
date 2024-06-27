"""A MultiKernelManager for use in the Jupyter server

- raises HTTPErrors
- creates REST API models
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import asyncio
import os
import pathlib  # noqa: TCH003
import typing as t
import warnings
from collections import defaultdict
from datetime import datetime, timedelta
from functools import partial, wraps

from jupyter_client.ioloop.manager import AsyncIOLoopKernelManager
from jupyter_client.multikernelmanager import AsyncMultiKernelManager, MultiKernelManager
from jupyter_client.session import Session
from jupyter_core.paths import exists
from jupyter_core.utils import ensure_async
from jupyter_events import EventLogger
from jupyter_events.schema_registry import SchemaRegistryException
from overrides import overrides
from tornado import web
from tornado.concurrent import Future
from tornado.ioloop import IOLoop, PeriodicCallback
from traitlets import (
    Any,
    Bool,
    Dict,
    Float,
    Instance,
    Integer,
    List,
    TraitError,
    Unicode,
    default,
    validate,
)

from jupyter_server import DEFAULT_EVENTS_SCHEMA_PATH
from jupyter_server._tz import isoformat, utcnow
from jupyter_server.prometheus.metrics import KERNEL_CURRENTLY_RUNNING_TOTAL
from jupyter_server.utils import ApiPath, import_item, to_os_path


class MappingKernelManager(MultiKernelManager):
    """A KernelManager that handles
    - File mapping
    - HTTP error handling
    - Kernel message filtering
    """

    @default("kernel_manager_class")
    def _default_kernel_manager_class(self):
        return "jupyter_client.ioloop.IOLoopKernelManager"

    kernel_argv = List(Unicode())

    root_dir = Unicode(config=True)

    _kernel_connections = Dict()

    _kernel_ports: dict[str, list[int]] = Dict()  # type: ignore[assignment]

    _culler_callback = None

    _initialized_culler = False

    @default("root_dir")
    def _default_root_dir(self):
        if not self.parent:
            return os.getcwd()
        return self.parent.root_dir

    @validate("root_dir")
    def _update_root_dir(self, proposal):
        """Do a bit of validation of the root dir."""
        value = proposal["value"]
        if not os.path.isabs(value):
            # If we receive a non-absolute path, make it absolute.
            value = os.path.abspath(value)
        if not exists(value) or not os.path.isdir(value):
            raise TraitError("kernel root dir %r is not a directory" % value)
        return value

    cull_idle_timeout = Integer(
        0,
        config=True,
        help="""Timeout (in seconds) after which a kernel is considered idle and ready to be culled.
        Values of 0 or lower disable culling. Very short timeouts may result in kernels being culled
        for users with poor network connections.""",
    )

    cull_interval_default = 300  # 5 minutes
    cull_interval = Integer(
        cull_interval_default,
        config=True,
        help="""The interval (in seconds) on which to check for idle kernels exceeding the cull timeout value.""",
    )

    cull_connected = Bool(
        False,
        config=True,
        help="""Whether to consider culling kernels which have one or more connections.
        Only effective if cull_idle_timeout > 0.""",
    )

    cull_busy = Bool(
        False,
        config=True,
        help="""Whether to consider culling kernels which are busy.
        Only effective if cull_idle_timeout > 0.""",
    )

    buffer_offline_messages = Bool(
        True,
        config=True,
        help="""Whether messages from kernels whose frontends have disconnected should be buffered in-memory.

        When True (default), messages are buffered and replayed on reconnect,
        avoiding lost messages due to interrupted connectivity.

        Disable if long-running kernels will produce too much output while
        no frontends are connected.
        """,
    )

    kernel_info_timeout = Float(
        60,
        config=True,
        help="""Timeout for giving up on a kernel (in seconds).

        On starting and restarting kernels, we check whether the
        kernel is running and responsive by sending kernel_info_requests.
        This sets the timeout in seconds for how long the kernel can take
        before being presumed dead.
        This affects the MappingKernelManager (which handles kernel restarts)
        and the ZMQChannelsHandler (which handles the startup).
        """,
    )

    _kernel_buffers = Any()

    @default("_kernel_buffers")
    def _default_kernel_buffers(self):
        return defaultdict(lambda: {"buffer": [], "session_key": "", "channels": {}})

    last_kernel_activity = Instance(
        datetime,
        help="The last activity on any kernel, including shutting down a kernel",
    )

    def __init__(self, **kwargs):
        """Initialize a kernel manager."""
        self.pinned_superclass = MultiKernelManager
        self._pending_kernel_tasks = {}
        self.pinned_superclass.__init__(self, **kwargs)
        self.last_kernel_activity = utcnow()

    allowed_message_types = List(
        trait=Unicode(),
        config=True,
        help="""White list of allowed kernel message types.
        When the list is empty, all message types are allowed.
        """,
    )

    allow_tracebacks = Bool(
        True, config=True, help=("Whether to send tracebacks to clients on exceptions.")
    )

    traceback_replacement_message = Unicode(
        "An exception occurred at runtime, which is not shown due to security reasons.",
        config=True,
        help=("Message to print when allow_tracebacks is False, and an exception occurs"),
    )

    # -------------------------------------------------------------------------
    # Methods for managing kernels and sessions
    # -------------------------------------------------------------------------

    def _handle_kernel_died(self, kernel_id):
        """notice that a kernel died"""
        self.log.warning("Kernel %s died, removing from map.", kernel_id)
        self.remove_kernel(kernel_id)

    def cwd_for_path(self, path, **kwargs):
        """Turn API path into absolute OS path."""
        os_path = to_os_path(path, self.root_dir)
        # in the case of documents and kernels not being on the same filesystem,
        # walk up to root_dir if the paths don't exist
        while not os.path.isdir(os_path) and os_path != self.root_dir:
            os_path = os.path.dirname(os_path)
        return os_path

    async def _remove_kernel_when_ready(self, kernel_id, kernel_awaitable):
        """Remove a kernel when it is ready."""
        await super()._remove_kernel_when_ready(kernel_id, kernel_awaitable)
        self._kernel_connections.pop(kernel_id, None)
        self._kernel_ports.pop(kernel_id, None)

    # TODO: DEC 2022: Revise the type-ignore once the signatures have been changed upstream
    # https://github.com/jupyter/jupyter_client/pull/905
    async def _async_start_kernel(  # type:ignore[override]
        self, *, kernel_id: str | None = None, path: ApiPath | None = None, **kwargs: str
    ) -> str:
        """Start a kernel for a session and return its kernel_id.

        Parameters
        ----------
        kernel_id : uuid (str)
            The uuid to associate the new kernel with. If this
            is not None, this kernel will be persistent whenever it is
            requested.
        path : API path
            The API path (unicode, '/' delimited) for the cwd.
            Will be transformed to an OS path relative to root_dir.
        kernel_name : str
            The name identifying which kernel spec to launch. This is ignored if
            an existing kernel is returned, but it may be checked in the future.
        """
        if kernel_id is None or kernel_id not in self:
            if path is not None:
                kwargs["cwd"] = self.cwd_for_path(path, env=kwargs.get("env", {}))
            if kernel_id is not None:
                assert kernel_id is not None, "Never Fail, but necessary for mypy "
                kwargs["kernel_id"] = kernel_id
            kernel_id = await self.pinned_superclass._async_start_kernel(self, **kwargs)
            self._kernel_connections[kernel_id] = 0
            task = asyncio.create_task(self._finish_kernel_start(kernel_id))
            if not getattr(self, "use_pending_kernels", None):
                await task
            else:
                self._pending_kernel_tasks[kernel_id] = task
            # add busy/activity markers:
            kernel = self.get_kernel(kernel_id)
            kernel.execution_state = "starting"  # type:ignore[attr-defined]
            kernel.reason = ""  # type:ignore[attr-defined]
            kernel.last_activity = utcnow()  # type:ignore[attr-defined]
            self.log.info("Kernel started: %s", kernel_id)
            self.log.debug("Kernel args: %r", kwargs)

            # Increase the metric of number of kernels running
            # for the relevant kernel type by 1
            KERNEL_CURRENTLY_RUNNING_TOTAL.labels(type=self._kernels[kernel_id].kernel_name).inc()

        else:
            self.log.info("Using existing kernel: %s", kernel_id)

        # Initialize culling if not already
        if not self._initialized_culler:
            self.initialize_culler()
        assert kernel_id is not None
        return kernel_id

    # see https://github.com/jupyter-server/jupyter_server/issues/1165
    # this assignment is technically incorrect, but might need a change of API
    # in jupyter_client.
    start_kernel = _async_start_kernel  # type:ignore[assignment]

    async def _finish_kernel_start(self, kernel_id):
        """Handle a kernel that finishes starting."""
        km = self.get_kernel(kernel_id)
        if hasattr(km, "ready"):
            ready = km.ready
            if not isinstance(ready, asyncio.Future):
                ready = asyncio.wrap_future(ready)
            try:
                await ready
            except Exception:
                self.log.exception("Error waiting for kernel manager ready")
                return

        self._kernel_ports[kernel_id] = km.ports
        self.start_watching_activity(kernel_id)
        # register callback for failed auto-restart
        self.add_restart_callback(
            kernel_id,
            lambda: self._handle_kernel_died(kernel_id),
            "dead",
        )

    def ports_changed(self, kernel_id):
        """Used by ZMQChannelsHandler to determine how to coordinate nudge and replays.

        Ports are captured when starting a kernel (via MappingKernelManager).  Ports
        are considered changed (following restarts) if the referenced KernelManager
        is using a set of ports different from those captured at startup.  If changes
        are detected, the captured set is updated and a value of True is returned.

        NOTE: Use is exclusive to ZMQChannelsHandler because this object is a singleton
        instance while ZMQChannelsHandler instances are per WebSocket connection that
        can vary per kernel lifetime.
        """
        changed_ports = self._get_changed_ports(kernel_id)
        if changed_ports:
            # If changed, update captured ports and return True, else return False.
            self.log.debug("Port change detected for kernel: %s", kernel_id)
            self._kernel_ports[kernel_id] = changed_ports
            return True
        return False

    def _get_changed_ports(self, kernel_id):
        """Internal method to test if a kernel's ports have changed and, if so, return their values.

        This method does NOT update the captured ports for the kernel as that can only be done
        by ZMQChannelsHandler, but instead returns the new list of ports if they are different
        than those captured at startup.  This enables the ability to conditionally restart
        activity monitoring immediately following a kernel's restart (if ports have changed).
        """
        # Get current ports and return comparison with ports captured at startup.
        km = self.get_kernel(kernel_id)
        assert isinstance(km.ports, list)
        assert isinstance(self._kernel_ports[kernel_id], list)
        if km.ports != self._kernel_ports[kernel_id]:
            return km.ports
        return None

    def start_buffering(self, kernel_id, session_key, channels):
        """Start buffering messages for a kernel

        Parameters
        ----------
        kernel_id : str
            The id of the kernel to stop buffering.
        session_key : str
            The session_key, if any, that should get the buffer.
            If the session_key matches the current buffered session_key,
            the buffer will be returned.
        channels : dict({'channel': ZMQStream})
            The zmq channels whose messages should be buffered.
        """

        if not self.buffer_offline_messages:
            for stream in channels.values():
                stream.close()
            return

        self.log.info("Starting buffering for %s", session_key)
        self._check_kernel_id(kernel_id)
        # clear previous buffering state
        self.stop_buffering(kernel_id)
        buffer_info = self._kernel_buffers[kernel_id]
        # record the session key because only one session can buffer
        buffer_info["session_key"] = session_key
        # TODO: the buffer should likely be a memory bounded queue, we're starting with a list to keep it simple
        buffer_info["buffer"] = []
        buffer_info["channels"] = channels

        # forward any future messages to the internal buffer
        def buffer_msg(channel, msg_parts):
            self.log.debug("Buffering msg on %s:%s", kernel_id, channel)
            buffer_info["buffer"].append((channel, msg_parts))

        for channel, stream in channels.items():
            stream.on_recv(partial(buffer_msg, channel))

    def get_buffer(self, kernel_id, session_key):
        """Get the buffer for a given kernel

        Parameters
        ----------
        kernel_id : str
            The id of the kernel to stop buffering.
        session_key : str, optional
            The session_key, if any, that should get the buffer.
            If the session_key matches the current buffered session_key,
            the buffer will be returned.
        """
        self.log.debug("Getting buffer for %s", kernel_id)
        if kernel_id not in self._kernel_buffers:
            return None

        buffer_info = self._kernel_buffers[kernel_id]
        if buffer_info["session_key"] == session_key:
            # remove buffer
            self._kernel_buffers.pop(kernel_id)
            # only return buffer_info if it's a match
            return buffer_info
        else:
            self.stop_buffering(kernel_id)

    def stop_buffering(self, kernel_id):
        """Stop buffering kernel messages

        Parameters
        ----------
        kernel_id : str
            The id of the kernel to stop buffering.
        """
        self.log.debug("Clearing buffer for %s", kernel_id)
        self._check_kernel_id(kernel_id)

        if kernel_id not in self._kernel_buffers:
            return
        buffer_info = self._kernel_buffers.pop(kernel_id)
        # close buffering streams
        for stream in buffer_info["channels"].values():
            if not stream.socket.closed:
                stream.on_recv(None)
                stream.close()

        msg_buffer = buffer_info["buffer"]
        if msg_buffer:
            self.log.info(
                "Discarding %s buffered messages for %s",
                len(msg_buffer),
                buffer_info["session_key"],
            )

    async def _async_shutdown_kernel(self, kernel_id, now=False, restart=False):
        """Shutdown a kernel by kernel_id"""
        self._check_kernel_id(kernel_id)

        # Decrease the metric of number of kernels
        # running for the relevant kernel type by 1
        KERNEL_CURRENTLY_RUNNING_TOTAL.labels(type=self._kernels[kernel_id].kernel_name).dec()

        if kernel_id in self._pending_kernel_tasks:
            task = self._pending_kernel_tasks.pop(kernel_id)
            task.cancel()

        self.stop_watching_activity(kernel_id)
        self.stop_buffering(kernel_id)

        return await self.pinned_superclass._async_shutdown_kernel(
            self, kernel_id, now=now, restart=restart
        )

    shutdown_kernel = _async_shutdown_kernel

    async def _async_restart_kernel(self, kernel_id, now=False):
        """Restart a kernel by kernel_id"""
        self._check_kernel_id(kernel_id)
        await self.pinned_superclass._async_restart_kernel(self, kernel_id, now=now)
        kernel = self.get_kernel(kernel_id)
        # return a Future that will resolve when the kernel has successfully restarted
        channel = kernel.connect_shell()
        future: Future[Any] = Future()

        def finish():
            """Common cleanup when restart finishes/fails for any reason."""
            if not channel.closed():  # type:ignore[operator]
                channel.close()
            loop.remove_timeout(timeout)
            kernel.remove_restart_callback(on_restart_failed, "dead")
            kernel._pending_restart_cleanup = None  # type:ignore[attr-defined]

        def on_reply(msg):
            self.log.debug("Kernel info reply received: %s", kernel_id)
            finish()
            if not future.done():
                future.set_result(msg)

        def on_timeout():
            self.log.warning("Timeout waiting for kernel_info_reply: %s", kernel_id)
            finish()
            if not future.done():
                future.set_exception(TimeoutError("Timeout waiting for restart"))

        def on_restart_failed():
            self.log.warning("Restarting kernel failed: %s", kernel_id)
            finish()
            if not future.done():
                future.set_exception(RuntimeError("Restart failed"))

        kernel.add_restart_callback(on_restart_failed, "dead")
        kernel._pending_restart_cleanup = finish  # type:ignore[attr-defined]
        kernel.session.send(channel, "kernel_info_request")
        channel.on_recv(on_reply)  # type:ignore[operator]
        loop = IOLoop.current()
        timeout = loop.add_timeout(loop.time() + self.kernel_info_timeout, on_timeout)
        # Re-establish activity watching if ports have changed...
        if self._get_changed_ports(kernel_id) is not None:
            self.stop_watching_activity(kernel_id)
            self.start_watching_activity(kernel_id)
        return future

    restart_kernel = _async_restart_kernel

    def notify_connect(self, kernel_id):
        """Notice a new connection to a kernel"""
        if kernel_id in self._kernel_connections:
            self._kernel_connections[kernel_id] += 1

    def notify_disconnect(self, kernel_id):
        """Notice a disconnection from a kernel"""
        if kernel_id in self._kernel_connections:
            self._kernel_connections[kernel_id] -= 1

    def kernel_model(self, kernel_id):
        """Return a JSON-safe dict representing a kernel

        For use in representing kernels in the JSON APIs.
        """
        self._check_kernel_id(kernel_id)
        kernel = self._kernels[kernel_id]

        model = {
            "id": kernel_id,
            "name": kernel.kernel_name,
            "last_activity": isoformat(kernel.last_activity),
            "execution_state": kernel.execution_state,
            "connections": self._kernel_connections.get(kernel_id, 0),
        }
        if getattr(kernel, "reason", None):
            model["reason"] = kernel.reason
        return model

    def list_kernels(self):
        """Returns a list of kernel_id's of kernels running."""
        kernels = []
        kernel_ids = self.pinned_superclass.list_kernel_ids(self)
        for kernel_id in kernel_ids:
            try:
                model = self.kernel_model(kernel_id)
                kernels.append(model)
            except (web.HTTPError, KeyError):
                # Probably due to a (now) non-existent kernel, continue building the list
                pass
        return kernels

    # override _check_kernel_id to raise 404 instead of KeyError
    def _check_kernel_id(self, kernel_id):
        """Check a that a kernel_id exists and raise 404 if not."""
        if kernel_id not in self:
            raise web.HTTPError(404, "Kernel does not exist: %s" % kernel_id)

    # monitoring activity:

    def start_watching_activity(self, kernel_id):
        """Start watching IOPub messages on a kernel for activity.

        - update last_activity on every message
        - record execution_state from status messages
        """
        kernel = self._kernels[kernel_id]
        # add busy/activity markers:
        kernel.execution_state = "starting"
        kernel.reason = ""
        kernel.last_activity = utcnow()
        kernel._activity_stream = kernel.connect_iopub()
        session = Session(
            config=kernel.session.config,
            key=kernel.session.key,
        )

        def record_activity(msg_list):
            """Record an IOPub message arriving from a kernel"""
            self.last_kernel_activity = kernel.last_activity = utcnow()

            idents, fed_msg_list = session.feed_identities(msg_list)
            msg = session.deserialize(fed_msg_list, content=False)

            msg_type = msg["header"]["msg_type"]
            if msg_type == "status":
                msg = session.deserialize(fed_msg_list)
                kernel.execution_state = msg["content"]["execution_state"]
                self.log.debug(
                    "activity on %s: %s (%s)",
                    kernel_id,
                    msg_type,
                    kernel.execution_state,
                )
            else:
                self.log.debug("activity on %s: %s", kernel_id, msg_type)

        kernel._activity_stream.on_recv(record_activity)

    def stop_watching_activity(self, kernel_id):
        """Stop watching IOPub messages on a kernel for activity."""
        kernel = self._kernels[kernel_id]
        if getattr(kernel, "_activity_stream", None):
            if not kernel._activity_stream.socket.closed:
                kernel._activity_stream.close()
            kernel._activity_stream = None
        if getattr(kernel, "_pending_restart_cleanup", None):
            kernel._pending_restart_cleanup()

    def initialize_culler(self):
        """Start idle culler if 'cull_idle_timeout' is greater than zero.

        Regardless of that value, set flag that we've been here.
        """
        if (
            not self._initialized_culler
            and self.cull_idle_timeout > 0
            and self._culler_callback is None
        ):
            _ = IOLoop.current()
            if self.cull_interval <= 0:  # handle case where user set invalid value
                self.log.warning(
                    "Invalid value for 'cull_interval' detected (%s) - using default value (%s).",
                    self.cull_interval,
                    self.cull_interval_default,
                )
                self.cull_interval = self.cull_interval_default
            self._culler_callback = PeriodicCallback(self.cull_kernels, 1000 * self.cull_interval)
            self.log.info(
                "Culling kernels with idle durations > %s seconds at %s second intervals ...",
                self.cull_idle_timeout,
                self.cull_interval,
            )
            if self.cull_busy:
                self.log.info("Culling kernels even if busy")
            if self.cull_connected:
                self.log.info("Culling kernels even with connected clients")
            self._culler_callback.start()

        self._initialized_culler = True

    async def cull_kernels(self):
        """Handle culling kernels."""
        self.log.debug(
            "Polling every %s seconds for kernels idle > %s seconds...",
            self.cull_interval,
            self.cull_idle_timeout,
        )
        """Create a separate list of kernels to avoid conflicting updates while iterating"""
        for kernel_id in list(self._kernels):
            try:
                await self.cull_kernel_if_idle(kernel_id)
            except Exception as e:
                self.log.exception(
                    "The following exception was encountered while checking the idle duration of kernel %s: %s",
                    kernel_id,
                    e,
                )

    async def cull_kernel_if_idle(self, kernel_id):
        """Cull a kernel if it is idle."""
        kernel = self._kernels[kernel_id]

        if getattr(kernel, "execution_state", None) == "dead":
            self.log.warning(
                "Culling '%s' dead kernel '%s' (%s).",
                kernel.execution_state,
                kernel.kernel_name,
                kernel_id,
            )
            await ensure_async(self.shutdown_kernel(kernel_id))
            return

        kernel_spec_metadata = kernel.kernel_spec.metadata
        cull_idle_timeout = kernel_spec_metadata.get("cull_idle_timeout", self.cull_idle_timeout)

        if hasattr(
            kernel, "last_activity"
        ):  # last_activity is monkey-patched, so ensure that has occurred
            self.log.debug(
                "kernel_id=%s, kernel_name=%s, last_activity=%s",
                kernel_id,
                kernel.kernel_name,
                kernel.last_activity,
            )
            dt_now = utcnow()
            dt_idle = dt_now - kernel.last_activity
            # Compute idle properties
            is_idle_time = dt_idle > timedelta(seconds=cull_idle_timeout)
            is_idle_execute = self.cull_busy or (kernel.execution_state != "busy")
            connections = self._kernel_connections.get(kernel_id, 0)
            is_idle_connected = self.cull_connected or not connections
            # Cull the kernel if all three criteria are met
            if is_idle_time and is_idle_execute and is_idle_connected:
                idle_duration = int(dt_idle.total_seconds())
                self.log.warning(
                    "Culling '%s' kernel '%s' (%s) with %d connections due to %s seconds of inactivity.",
                    kernel.execution_state,
                    kernel.kernel_name,
                    kernel_id,
                    connections,
                    idle_duration,
                )
                await ensure_async(self.shutdown_kernel(kernel_id))


# AsyncMappingKernelManager inherits as much as possible from MappingKernelManager,
# overriding only what is different.
class AsyncMappingKernelManager(MappingKernelManager, AsyncMultiKernelManager):  # type:ignore[misc]
    """An asynchronous mapping kernel manager."""

    @default("kernel_manager_class")
    def _default_kernel_manager_class(self):
        return "jupyter_server.services.kernels.kernelmanager.ServerKernelManager"

    @validate("kernel_manager_class")
    def _validate_kernel_manager_class(self, proposal):
        """A validator for the kernel manager class."""
        km_class_value = proposal.value
        km_class = import_item(km_class_value)
        if not issubclass(km_class, ServerKernelManager):
            warnings.warn(
                f"KernelManager class '{km_class}' is not a subclass of 'ServerKernelManager'.  Custom "
                "KernelManager classes should derive from 'ServerKernelManager' beginning with jupyter-server 2.0 "
                "or risk missing functionality.  Continuing...",
                FutureWarning,
                stacklevel=3,
            )
        return km_class_value

    def __init__(self, **kwargs):
        """Initialize an async mapping kernel manager."""
        self.pinned_superclass = MultiKernelManager
        self._pending_kernel_tasks = {}
        self.pinned_superclass.__init__(self, **kwargs)
        self.last_kernel_activity = utcnow()


def emit_kernel_action_event(success_msg: str = "") -> t.Callable[..., t.Any]:
    """Decorate kernel action methods to
    begin emitting jupyter kernel action events.

    Parameters
    ----------
    success_msg: str
        A formattable string that's passed to the message field of
        the emitted event when the action succeeds. You can include
        the kernel_id, kernel_name, or action in the message using
        a formatted string argument,
        e.g. "{kernel_id} succeeded to {action}."

    error_msg: str
        A formattable string that's passed to the message field of
        the emitted event when the action fails. You can include
        the kernel_id, kernel_name, or action in the message using
        a formatted string argument,
        e.g. "{kernel_id} failed to {action}."
    """

    def wrap_method(method):
        @wraps(method)
        async def wrapped_method(self, *args, **kwargs):
            """"""
            # Get the method name from the
            action = method.__name__.replace("_kernel", "")
            # If the method succeeds, emit a success event.
            try:
                out = await method(self, *args, **kwargs)
                data = {
                    "kernel_name": self.kernel_name,
                    "action": action,
                    "status": "success",
                    "msg": success_msg.format(
                        kernel_id=self.kernel_id, kernel_name=self.kernel_name, action=action
                    ),
                }
                if self.kernel_id:
                    data["kernel_id"] = self.kernel_id
                self.emit(
                    schema_id="https://events.jupyter.org/jupyter_server/kernel_actions/v1",
                    data=data,
                )
                return out
            # If the method fails, emit a failed event.
            except Exception as err:
                data = {
                    "kernel_name": self.kernel_name,
                    "action": action,
                    "status": "error",
                    "msg": str(err),
                }
                if self.kernel_id:
                    data["kernel_id"] = self.kernel_id
                # If the exception is an HTTPError (usually via a gateway request)
                # log the status_code and HTTPError log_message.
                if isinstance(err, web.HTTPError):
                    msg = err.log_message or ""
                    data["status_code"] = err.status_code
                    data["msg"] = msg
                self.emit(
                    schema_id="https://events.jupyter.org/jupyter_server/kernel_actions/v1",
                    data=data,
                )
                raise err

        return wrapped_method

    return wrap_method


class ServerKernelManager(AsyncIOLoopKernelManager):
    """A server-specific kernel manager."""

    # Define activity-related attributes:
    execution_state = Unicode(
        None, allow_none=True, help="The current execution state of the kernel"
    )
    reason = Unicode("", help="The reason for the last failure against the kernel")

    last_activity = Instance(datetime, help="The last activity on the kernel")

    # A list of pathlib objects, each pointing at an event
    # schema to register with this kernel manager's eventlogger.
    # This trait should not be overridden.
    @property
    def core_event_schema_paths(self) -> list[pathlib.Path]:
        return [DEFAULT_EVENTS_SCHEMA_PATH / "kernel_actions" / "v1.yaml"]

    # This trait is intended for subclasses to override and define
    # custom event schemas.
    extra_event_schema_paths: List[str] = List(
        default_value=[],
        help="""
        A list of pathlib.Path objects pointing at to register with
        the kernel manager's eventlogger.
        """,
    ).tag(config=True)

    event_logger = Instance(EventLogger)

    @default("event_logger")
    def _default_event_logger(self):
        """Initialize the logger and ensure all required events are present."""
        if (
            self.parent is not None
            and self.parent.parent is not None
            and hasattr(self.parent.parent, "event_logger")
        ):
            logger = self.parent.parent.event_logger
        else:
            # If parent does not have an event logger, create one.
            logger = EventLogger()
        # Ensure that all the expected schemas are registered. If not, register them.
        schemas = self.core_event_schema_paths + self.extra_event_schema_paths
        for schema_path in schemas:
            # Try registering the event.
            try:
                logger.register_event_schema(schema_path)
            # Pass if it already exists.
            except SchemaRegistryException:
                pass
        return logger

    def emit(self, schema_id, data):
        """Emit an event from the kernel manager."""
        self.event_logger.emit(schema_id=schema_id, data=data)

    @overrides
    @emit_kernel_action_event(
        success_msg="Kernel {kernel_id} was started.",
    )
    async def start_kernel(self, *args, **kwargs):
        return await super().start_kernel(*args, **kwargs)

    @overrides
    @emit_kernel_action_event(
        success_msg="Kernel {kernel_id} was shutdown.",
    )
    async def shutdown_kernel(self, *args, **kwargs):
        return await super().shutdown_kernel(*args, **kwargs)

    @overrides
    @emit_kernel_action_event(
        success_msg="Kernel {kernel_id} was restarted.",
    )
    async def restart_kernel(self, *args, **kwargs):
        return await super().restart_kernel(*args, **kwargs)

    @overrides
    @emit_kernel_action_event(
        success_msg="Kernel {kernel_id} was interrupted.",
    )
    async def interrupt_kernel(self, *args, **kwargs):
        return await super().interrupt_kernel(*args, **kwargs)
