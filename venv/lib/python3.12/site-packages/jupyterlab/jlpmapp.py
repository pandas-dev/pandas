"""A Jupyter-aware wrapper for the yarn package manager"""

import os

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import sys

from jupyterlab_server.process import subprocess, which

HERE = os.path.dirname(os.path.abspath(__file__))
YARN_PATH = os.path.join(HERE, "staging", "yarn.js")


def execvp(cmd, argv):
    """Execvp, except on Windows where it uses Popen.

    The first argument, by convention, should point to the filename
    associated with the file being executed.

    Python provides execvp on Windows, but its behavior is problematic
    (Python bug#9148).
    """
    cmd = which(cmd)
    if os.name == "nt":
        import signal
        import sys

        p = subprocess.Popen([cmd] + argv[1:])
        # Don't raise KeyboardInterrupt in the parent process.
        # Set this after spawning, to avoid subprocess inheriting handler.
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        p.wait()
        sys.exit(p.returncode)
    else:
        os.execvp(cmd, argv)  # noqa S606


def main(argv=None):
    """Run node and return the result."""
    # Make sure node is available.
    argv = argv or sys.argv[1:]
    execvp("node", ["node", YARN_PATH, *argv])
