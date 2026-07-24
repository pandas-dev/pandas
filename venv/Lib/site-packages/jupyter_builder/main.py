"""Jupyter Builder CLI entry point."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import sys

from jupyter_core.application import JupyterApp

from jupyter_builder.extension_commands.build import BuildLabExtensionApp
from jupyter_builder.extension_commands.develop import DevelopLabExtensionApp
from jupyter_builder.extension_commands.watch import WatchLabExtensionApp

from ._version import __version__

_EXAMPLES = """
jupyter-builder build                       # build a prebuilt labextension
jupyter-builder develop                     # develop a prebuilt labextension
jupyter-builder watch                       # watch a prebuilt labextension
"""


class BuilderApp(JupyterApp):
    """Base jupyter-builder command entry point."""

    name = "jupyter builder"
    version = __version__
    description = "Build JupyterLab extensions"
    examples = _EXAMPLES

    subcommands = {
        "develop": (DevelopLabExtensionApp, "Develop labextension(s)"),
        "build": (BuildLabExtensionApp, "Build labextension"),
        "watch": (WatchLabExtensionApp, "Watch labextension"),
    }

    def start(self) -> None:
        """Perform the App's functions as configured."""
        super().start()

        # The above should have called a subcommand and raised NoStart; if we
        # get here, it didn't, so we should self.log.info a message.
        subcmds = ", ".join(sorted(self.subcommands))
        self.exit(f"Please supply at least one subcommand: {subcmds}")


main = BuilderApp.launch_instance

if __name__ == "__main__":
    sys.exit(main())
