""" A session for managing a language server process
"""

import asyncio
import atexit
import os
import string
import subprocess
from datetime import datetime, timezone

from tornado.ioloop import IOLoop
from tornado.queues import Queue
from tornado.websocket import WebSocketHandler
from traitlets import Bunch, Instance, Set, Unicode, UseEnum, observe
from traitlets.config import LoggingConfigurable

from . import stdio
from .schema import LANGUAGE_SERVER_SPEC
from .specs.utils import censored_spec
from .trait_types import Schema
from .types import SessionStatus


class LanguageServerSession(LoggingConfigurable):
    """Manage a session for a connection to a language server"""

    language_server = Unicode(help="the language server implementation name")
    spec = Schema(LANGUAGE_SERVER_SPEC)

    # run-time specifics
    process = Instance(
        subprocess.Popen, help="the language server subprocess", allow_none=True
    )
    writer = Instance(stdio.LspStdIoWriter, help="the JSON-RPC writer", allow_none=True)
    reader = Instance(stdio.LspStdIoReader, help="the JSON-RPC reader", allow_none=True)
    from_lsp = Instance(
        Queue, help="a queue for string messages from the server", allow_none=True
    )
    to_lsp = Instance(
        Queue, help="a queue for string message to the server", allow_none=True
    )
    handlers = Set(
        trait=Instance(WebSocketHandler),
        default_value=[],
        help="the currently subscribed websockets",
    )
    status = UseEnum(SessionStatus, default_value=SessionStatus.NOT_STARTED)
    last_handler_message_at = Instance(datetime, allow_none=True)
    last_server_message_at = Instance(datetime, allow_none=True)

    _tasks = None

    _skip_serialize = ["argv", "debug_argv"]

    def __init__(self, *args, **kwargs):
        """set up the required traitlets and exit behavior for a session"""
        super().__init__(*args, **kwargs)
        atexit.register(self.stop)

    def __repr__(self):  # pragma: no cover
        return (
            "<LanguageServerSession(" "language_server={language_server}, argv={argv})>"
        ).format(language_server=self.language_server, **self.spec)

    def to_json(self):
        return dict(
            handler_count=len(self.handlers),
            status=self.status.value,
            last_server_message_at=(
                self.last_server_message_at.isoformat()
                if self.last_server_message_at
                else None
            ),
            last_handler_message_at=(
                self.last_handler_message_at.isoformat()
                if self.last_handler_message_at
                else None
            ),
            spec=censored_spec(self.spec),
        )

    def initialize(self):
        """(re)initialize a language server session"""
        self.stop()
        self.status = SessionStatus.STARTING
        self.init_queues()
        self.init_process()
        self.init_writer()
        self.init_reader()

        loop = asyncio.get_event_loop()
        self._tasks = [
            loop.create_task(coro())
            for coro in [self._read_lsp, self._write_lsp, self._broadcast_from_lsp]
        ]

        self.status = SessionStatus.STARTED

    def stop(self):
        """clean up all of the state of the session"""

        self.status = SessionStatus.STOPPING

        if self.process:
            self.process.terminate()
            self.process = None
        if self.reader:
            self.reader.close()
            self.reader = None
        if self.writer:
            self.writer.close()
            self.writer = None

        if self._tasks:
            [task.cancel() for task in self._tasks]

        self.status = SessionStatus.STOPPED

    @observe("handlers")
    def _on_handlers(self, change: Bunch):
        """re-initialize if someone starts listening, or stop if nobody is"""
        if change["new"] and not self.process:
            self.initialize()
        elif not change["new"] and self.process:
            self.stop()

    def write(self, message):
        """wrapper around the write queue to keep it mostly internal"""
        self.last_handler_message_at = self.now()
        IOLoop.current().add_callback(self.to_lsp.put_nowait, message)

    def now(self):
        return datetime.now(timezone.utc)

    def init_process(self):
        """start the language server subprocess"""
        self.process = subprocess.Popen(
            self.spec["argv"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            env=self.substitute_env(self.spec.get("env", {}), os.environ),
            bufsize=0,
        )

    def init_queues(self):
        """create the queues"""
        self.from_lsp = Queue()
        self.to_lsp = Queue()

    def init_reader(self):
        """create the stdout reader (from the language server)"""
        self.reader = stdio.LspStdIoReader(
            stream=self.process.stdout, queue=self.from_lsp, parent=self
        )

    def init_writer(self):
        """create the stdin writer (to the language server)"""
        self.writer = stdio.LspStdIoWriter(
            stream=self.process.stdin, queue=self.to_lsp, parent=self
        )

    def substitute_env(self, env, base):
        final_env = base.copy()

        for key, value in env.items():
            final_env.update({key: string.Template(value).safe_substitute(base)})

        return final_env

    async def _read_lsp(self):
        await self.reader.read()

    async def _write_lsp(self):
        await self.writer.write()

    async def _broadcast_from_lsp(self):
        """loop for reading messages from the queue of messages from the language
        server
        """
        async for message in self.from_lsp:
            self.last_server_message_at = self.now()
            await self.parent.on_server_message(message, self)
            self.from_lsp.task_done()
