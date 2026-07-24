"""A Jupyter-aware wrapper for the yarn package manager."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from shutil import which

HERE = str(Path(__file__).resolve().parent)
YARN_PATH = str(Path(__file__).resolve().parent / "yarn.js")


def _which_node_js(env: dict[str, str] | None = None) -> str:
    """Get the full path to Node.js executable.

    Parameters
    ----------
    env: dict, optional
        The environment variables, defaults to `os.environ`.

    """
    env = env or os.environ  # type:ignore[assignment]
    path = env.get("PATH") or os.defpath  # type:ignore[union-attr]
    command_with_path = which("node", path=path)

    # Allow nodejs as an alias to node.
    if not command_with_path:
        command_with_path = which("nodejs", path=path)

    if not command_with_path:
        msg = (
            "Please install Node.js and npm before continuing installation. "
            "You may be able to install Node.js from your package manager, "
            "from conda, or directly from the Node.js website "
            "(https://nodejs.org)."
        )
        raise ValueError(msg)
    # Return the path as found on PATH, without resolving symlinks: version
    # managers such as Volta install `node` as a symlink to a dispatcher
    # binary that selects the tool from the name it is invoked under, so the
    # dispatch only works when the command keeps the `node` name.
    return command_with_path


def _execvp_node(argv: list[str]) -> None:
    """Execute node, using Popen on Windows due to Python bug #9148."""
    cmd = _which_node_js()
    if os.name == "nt":
        import signal  # noqa: PLC0415

        p = subprocess.Popen([cmd, *argv[1:]])  # noqa: S603
        # Don't raise KeyboardInterrupt in the parent process.
        # Set this after spawning, to avoid subprocess inheriting handler.
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        p.wait()
        sys.exit(p.returncode)
    else:
        os.execvp(cmd, argv)  # noqa: S606


def main(argv: list[str] | None = None) -> None:
    """Run node and return the result."""
    # Make sure node is available.
    argv = argv or sys.argv[1:]
    # Disable Yarn telemetry by default for jlpm
    os.environ.setdefault("YARN_ENABLE_TELEMETRY", "0")
    _execvp_node(["node", YARN_PATH, *argv])
