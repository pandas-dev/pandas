"""An output widget mimic."""
from __future__ import annotations

from typing import Any

from jupyter_client.client import KernelClient
from nbformat import NotebookNode
from nbformat.v4 import output_from_msg

from .jsonutil import json_clean


class OutputWidget:
    """This class mimics a front end output widget"""

    def __init__(
        self, comm_id: str, state: dict[str, Any], kernel_client: KernelClient, executor: Any
    ) -> None:
        """Initialize the widget."""
        self.comm_id: str = comm_id
        self.state: dict[str, Any] = state
        self.kernel_client: KernelClient = kernel_client
        self.executor = executor
        self.topic: bytes = ("comm-%s" % self.comm_id).encode("ascii")
        self.outputs: list[NotebookNode] = self.state["outputs"]
        self.clear_before_next_output: bool = False

    def clear_output(self, outs: list[NotebookNode], msg: dict[str, Any], cell_index: int) -> None:
        """Clear output."""
        self.parent_header = msg["parent_header"]
        content = msg["content"]
        if content.get("wait"):
            self.clear_before_next_output = True
        else:
            self.outputs = []
            # sync back the state to the kernel
            self.sync_state()
            if hasattr(self.executor, "widget_state"):
                # sync the state to the nbconvert state as well, since that is used for testing
                self.executor.widget_state[self.comm_id]["outputs"] = self.outputs

    def sync_state(self) -> None:
        """Sync state."""
        state = {"outputs": self.outputs}
        msg = {"method": "update", "state": state, "buffer_paths": []}
        self.send(msg)

    def _publish_msg(
        self,
        msg_type: str,
        data: dict[str, Any] | None = None,
        metadata: dict[str, Any] | None = None,
        buffers: list[Any] | None = None,
        **keys: Any,
    ) -> None:
        """Helper for sending a comm message on IOPub"""
        data = {} if data is None else data
        metadata = {} if metadata is None else metadata
        content = json_clean(dict(data=data, comm_id=self.comm_id, **keys))
        msg = self.kernel_client.session.msg(
            msg_type, content=content, parent=self.parent_header, metadata=metadata
        )
        self.kernel_client.shell_channel.send(msg)

    def send(
        self,
        data: dict[str, Any] | None = None,
        metadata: dict[str, Any] | None = None,
        buffers: list[Any] | None = None,
    ) -> None:
        """Send a comm message."""
        self._publish_msg("comm_msg", data=data, metadata=metadata, buffers=buffers)

    def output(
        self, outs: list[NotebookNode], msg: dict[str, Any], display_id: str, cell_index: int
    ) -> None:
        """Handle output."""
        if self.clear_before_next_output:
            self.outputs = []
            self.clear_before_next_output = False
        self.parent_header = msg["parent_header"]
        output = output_from_msg(msg)  # type:ignore[no-untyped-call]

        if self.outputs:
            # try to coalesce/merge output text
            last_output = self.outputs[-1]
            if (
                last_output["output_type"] == "stream"
                and output["output_type"] == "stream"
                and last_output["name"] == output["name"]
            ):
                last_output["text"] += output["text"]
            else:
                self.outputs.append(output)
        else:
            self.outputs.append(output)
        self.sync_state()
        if hasattr(self.executor, "widget_state"):
            # sync the state to the nbconvert state as well, since that is used for testing
            self.executor.widget_state[self.comm_id]["outputs"] = self.outputs

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the state."""
        if "msg_id" in state:
            msg_id = state.get("msg_id")
            if msg_id:
                self.executor.register_output_hook(msg_id, self)
                self.msg_id = msg_id
            else:
                self.executor.remove_output_hook(self.msg_id, self)
                self.msg_id = msg_id

    def handle_msg(self, msg: dict[str, Any]) -> None:
        """Handle a message."""
        content = msg["content"]
        comm_id = content["comm_id"]
        if comm_id != self.comm_id:
            raise AssertionError("Mismatched comm id")
        data = content["data"]
        if "state" in data:
            self.set_state(data["state"])
