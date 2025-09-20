"""Exceptions for nbclient."""
from __future__ import annotations

from typing import Any

from nbformat import NotebookNode


class CellControlSignal(Exception):  # noqa
    """
    A custom exception used to indicate that the exception is used for cell
    control actions (not the best model, but it's needed to cover existing
    behavior without major refactors).
    """

    pass


class CellTimeoutError(TimeoutError, CellControlSignal):
    """
    A custom exception to capture when a cell has timed out during execution.
    """

    @classmethod
    def error_from_timeout_and_cell(
        cls, msg: str, timeout: int, cell: NotebookNode
    ) -> CellTimeoutError:
        """Create an error from a timeout on a cell."""
        if cell and cell.source:
            src_by_lines = cell.source.strip().split("\n")
            src = (
                cell.source
                if len(src_by_lines) < 11
                else f"{src_by_lines[:5]}\n...\n{src_by_lines[-5:]}"
            )
        else:
            src = "Cell contents not found."
        return cls(timeout_err_msg.format(timeout=timeout, msg=msg, cell_contents=src))


class DeadKernelError(RuntimeError):
    """A dead kernel error."""

    pass


class CellExecutionComplete(CellControlSignal):
    """
    Used as a control signal for cell execution across execute_cell and
    process_message function calls. Raised when all execution requests
    are completed and no further messages are expected from the kernel
    over zeromq channels.
    """

    pass


class CellExecutionError(CellControlSignal):
    """
    Custom exception to propagate exceptions that are raised during
    notebook execution to the caller. This is mostly useful when
    using nbconvert as a library, since it allows to deal with
    failures gracefully.
    """

    def __init__(self, traceback: str, ename: str, evalue: str) -> None:
        """Initialize the error."""
        super().__init__(traceback)
        self.traceback = traceback
        self.ename = ename
        self.evalue = evalue

    def __reduce__(self) -> tuple[Any]:
        """Reduce implementation."""
        return type(self), (self.traceback, self.ename, self.evalue)  # type:ignore[return-value]

    def __str__(self) -> str:
        """Str repr."""
        if self.traceback:
            return self.traceback
        else:
            return f"{self.ename}: {self.evalue}"

    @classmethod
    def from_cell_and_msg(cls, cell: NotebookNode, msg: dict[str, Any]) -> CellExecutionError:
        """Instantiate from a code cell object and a message contents
        (message is either execute_reply or error)
        """

        # collect stream outputs for our error message
        stream_outputs: list[str] = []
        for output in cell.outputs:
            if output["output_type"] == "stream":
                stream_outputs.append(
                    stream_output_msg.format(name=output["name"], text=output["text"].rstrip())
                )
        if stream_outputs:
            # add blank line before, trailing separator
            # if there is any stream output to display
            stream_outputs.insert(0, "")
            stream_outputs.append("------------------")
        stream_output: str = "\n".join(stream_outputs)

        tb = "\n".join(msg.get("traceback", []) or [])
        return cls(
            exec_err_msg.format(
                cell=cell,
                stream_output=stream_output,
                traceback=tb,
            ),
            ename=msg.get("ename", "<Error>"),
            evalue=msg.get("evalue", ""),
        )


stream_output_msg: str = """\
----- {name} -----
{text}"""

exec_err_msg: str = """\
An error occurred while executing the following cell:
------------------
{cell.source}
------------------
{stream_output}

{traceback}
"""


timeout_err_msg: str = """\
A cell timed out while it was being executed, after {timeout} seconds.
The message was: {msg}.
Here is a preview of the cell contents:
-------------------
{cell_contents}
-------------------
"""
