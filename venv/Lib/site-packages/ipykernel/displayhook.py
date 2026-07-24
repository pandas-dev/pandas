"""Replacements for sys.displayhook that publish over ZMQ."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import builtins
import sys
import threading
import typing as t
from contextvars import ContextVar

from IPython.core.displayhook import DisplayHook
from jupyter_client.session import Session, extract_header
from traitlets import Any, Instance, default

from ipykernel.jsonutil import encode_images, json_clean


class ZMQDisplayHook:
    """A simple displayhook that publishes the object's repr over a ZeroMQ
    socket."""

    topic = b"execute_result"

    def __init__(self, session, pub_socket):
        """Initialize the hook."""
        self.session = session
        self.pub_socket = pub_socket

        self._parent_header: ContextVar[dict[str, Any]] = ContextVar("parent_header")
        self._parent_header.set({})
        self._parent_header_global = {}

    def get_execution_count(self):
        """This method is replaced in kernelapp"""
        return 0

    def __call__(self, obj):
        """Handle a hook call."""
        if obj is None:
            return

        builtins._ = obj  # type:ignore[attr-defined]
        sys.stdout.flush()
        sys.stderr.flush()
        contents = {
            "execution_count": self.get_execution_count(),
            "data": {"text/plain": repr(obj)},
            "metadata": {},
        }
        self.session.send(
            self.pub_socket,
            "execute_result",
            contents,
            parent=self.parent_header,
            ident=self.topic,
        )

    @property
    def parent_header(self):
        try:
            return self._parent_header.get()
        except LookupError:
            return self._parent_header_global

    def set_parent(self, parent):
        """Set the parent header."""
        parent_header = extract_header(parent)
        self._parent_header.set(parent_header)
        self._parent_header_global = parent_header


class ZMQShellDisplayHook(DisplayHook):
    """A displayhook subclass that publishes data using ZeroMQ. This is intended
    to work with an InteractiveShell instance. It sends a dict of different
    representations of the object."""

    topic = None

    session = Instance(Session, allow_none=True)
    pub_socket = Any(allow_none=True)
    _parent_header: ContextVar[dict[str, Any]]
    _thread_local = Any()
    msg: dict[str, t.Any] | None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent_header = ContextVar("parent_header")
        self._parent_header.set({})

    @default("_thread_local")
    def _default_thread_local(self):
        return threading.local()

    @property
    def _hooks(self):
        if not hasattr(self._thread_local, "hooks"):
            self._thread_local.hooks = []
        return self._thread_local.hooks

    def register_hook(self, hook):
        """Register a transform hook on the execute_result message.

        Mirrors ``ZMQDisplayPublisher.register_hook``. Each hook receives the
        outbound message dict and must return either a (possibly mutated)
        message dict to continue the chain, or ``None`` to suppress the send.
        """
        self._hooks.append(hook)

    def unregister_hook(self, hook):
        """Remove a previously registered hook. Returns True on success."""
        try:
            self._hooks.remove(hook)
            return True
        except ValueError:
            return False

    @property
    def parent_header(self):
        try:
            return self._parent_header.get()
        except LookupError:
            return self._parent_header_global

    def set_parent(self, parent):
        """Set the parent header."""
        parent_header = extract_header(parent)
        self._parent_header.set(parent_header)
        self._parent_header_global = parent_header

    def start_displayhook(self):
        """Start the display hook."""
        if self.session:
            self.msg = self.session.msg(
                "execute_result",
                {
                    "data": {},
                    "metadata": {},
                },
                parent=self.parent_header,
            )

    def write_output_prompt(self):
        """Write the output prompt."""
        if self.msg:
            self.msg["content"]["execution_count"] = self.prompt_count

    def write_format_data(self, format_dict, md_dict=None):
        """Write format data to the message."""
        if self.msg:
            self.msg["content"]["data"] = json_clean(encode_images(format_dict))
            self.msg["content"]["metadata"] = md_dict

    def finish_displayhook(self):
        """Finish up all displayhook activities.

        Runs the registered hook chain before ``session.send``. Each hook
        either returns a message (to continue) or ``None`` (to suppress the
        send). This mirrors the transform pipeline on
        ``ZMQDisplayPublisher.publish`` so a single hook implementation can
        attach to both the ``display_data`` and ``execute_result`` paths.
        """
        sys.stdout.flush()
        sys.stderr.flush()
        if self.msg and self.msg["content"]["data"] and self.session:
            msg = self.msg
            for hook in self._hooks:
                msg = hook(msg)
                if msg is None:
                    self.msg = None
                    return
            self.session.send(self.pub_socket, msg, ident=self.topic)
        self.msg = None
