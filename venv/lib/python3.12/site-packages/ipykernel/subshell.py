"""A thread for a subshell."""

import asyncio
from typing import Any

import zmq

from .socket_pair import SocketPair
from .thread import BaseThread


class SubshellThread(BaseThread):
    """A thread for a subshell.

    .. versionadded:: 7
    """

    def __init__(
        self,
        subshell_id: str,
        context: zmq.Context[Any],
        **kwargs,
    ):
        """Initialize the thread."""
        super().__init__(name=f"subshell-{subshell_id}", **kwargs)

        self.shell_channel_to_subshell = SocketPair(context, subshell_id)
        self.subshell_to_shell_channel = SocketPair(context, subshell_id + "-reverse")

        # When aborting flag is set, execute_request messages to this subshell will be aborted.
        self.aborting = False

        self.asyncio_lock = asyncio.Lock()

    def run(self) -> None:
        """Run the thread."""
        try:
            super().run()
        finally:
            self.shell_channel_to_subshell.close()
            self.subshell_to_shell_channel.close()
