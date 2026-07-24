"""nbclient cli."""
from __future__ import annotations

import logging
import sys
import typing
from pathlib import Path
from textwrap import dedent

import nbformat
from jupyter_core.application import JupyterApp
from traitlets import Bool, Integer, List, Unicode, default
from traitlets.config import catch_config_error

from nbclient import __version__

from .client import NotebookClient

# mypy: disable-error-code="no-untyped-call"

nbclient_aliases: dict[str, str] = {
    "timeout": "NbClientApp.timeout",
    "startup_timeout": "NbClientApp.startup_timeout",
    "kernel_name": "NbClientApp.kernel_name",
    "output": "NbClientApp.output_base",
}

nbclient_flags: dict[str, typing.Any] = {
    "allow-errors": (
        {
            "NbClientApp": {
                "allow_errors": True,
            },
        },
        "Errors are ignored and execution is continued until the end of the notebook.",
    ),
    "inplace": (
        {
            "NbClientApp": {
                "inplace": True,
            },
        },
        "Overwrite input notebook with executed results.",
    ),
}


class NbClientApp(JupyterApp):
    """
    An application used to execute notebook files (``*.ipynb``)
    """

    version = Unicode(__version__)
    name = "jupyter-execute"
    aliases = nbclient_aliases
    flags = nbclient_flags

    description = "An application used to execute notebook files (*.ipynb)"
    notebooks = List(Unicode(), help="Path of notebooks to convert").tag(config=True)
    timeout = Integer(
        None,
        allow_none=True,
        help=dedent(
            """
            The time to wait (in seconds) for output from executions.
            If a cell execution takes longer, a TimeoutError is raised.
            ``-1`` will disable the timeout.
            """
        ),
    ).tag(config=True)
    startup_timeout = Integer(
        60,
        help=dedent(
            """
            The time to wait (in seconds) for the kernel to start.
            If kernel startup takes longer, a RuntimeError is
            raised.
            """
        ),
    ).tag(config=True)
    allow_errors = Bool(
        False,
        help=dedent(
            """
            When a cell raises an error the default behavior is that
            execution is stopped and a :py:class:`nbclient.exceptions.CellExecutionError`
            is raised.
            If this flag is provided, errors are ignored and execution
            is continued until the end of the notebook.
            """
        ),
    ).tag(config=True)
    skip_cells_with_tag = Unicode(
        "skip-execution",
        help=dedent(
            """
            Name of the cell tag to use to denote a cell that should be skipped.
            """
        ),
    ).tag(config=True)
    kernel_name = Unicode(
        "",
        help=dedent(
            """
            Name of kernel to use to execute the cells.
            If not set, use the kernel_spec embedded in the notebook.
            """
        ),
    ).tag(config=True)
    inplace = Bool(
        False,
        help=dedent(
            """
            Default is execute notebook without writing the newly executed notebook.
            If this flag is provided, the newly generated notebook will
            overwrite the input notebook.
            """
        ),
    ).tag(config=True)
    output_base = Unicode(
        None,
        allow_none=True,
        help=dedent(
            """
            Write executed notebook to this file base name.
            Supports pattern replacements ``'{notebook_name}'``,
            the name of the input notebook file without extension.
            Note that output is always relative to the parent directory of the
            input notebook.
            """
        ),
    ).tag(config=True)

    @default("log_level")
    def _log_level_default(self) -> int:
        return logging.INFO

    @catch_config_error
    def initialize(self, argv: list[str] | None = None) -> None:
        """Initialize the app."""
        super().initialize(argv)

        # Get notebooks to run
        self.notebooks = self.get_notebooks()

        # If there are none, throw an error
        if not self.notebooks:
            sys.exit(-1)

        # If output, must have single notebook
        if len(self.notebooks) > 1 and self.output_base is not None:
            if "{notebook_name}" not in self.output_base:
                msg = (
                    "If passing multiple notebooks with `--output=output` option, "
                    "output string must contain {notebook_name}"
                )
                raise ValueError(msg)

        # Loop and run them one by one
        for path in self.notebooks:
            self.run_notebook(path)

    def get_notebooks(self) -> list[str]:
        """Get the notebooks for the app."""
        # If notebooks were provided from the command line, use those
        if self.extra_args:
            notebooks = self.extra_args
        # If not, look to the class attribute
        else:
            notebooks = self.notebooks

        # Return what we got.
        return notebooks

    def run_notebook(self, notebook_path: str) -> None:
        """Run a notebook by path."""
        # Log it
        self.log.info(f"Executing {notebook_path}")

        input_path = Path(notebook_path).with_suffix(".ipynb")

        # Get its parent directory so we can add it to the $PATH
        path = input_path.parent.absolute()

        # Optional output of executed notebook
        if self.inplace:
            output_path = input_path
        elif self.output_base:
            output_path = input_path.parent.joinpath(
                self.output_base.format(notebook_name=input_path.with_suffix("").name)
            ).with_suffix(".ipynb")
        else:
            output_path = None

        if output_path and not output_path.parent.is_dir():
            msg = f"Cannot write to directory={output_path.parent} that does not exist"
            raise ValueError(msg)

        # Open up the notebook we're going to run
        with input_path.open() as f:
            nb = nbformat.read(f, as_version=4)

        # Configure nbclient to run the notebook
        client = NotebookClient(
            nb,
            timeout=self.timeout,
            startup_timeout=self.startup_timeout,
            skip_cells_with_tag=self.skip_cells_with_tag,
            allow_errors=self.allow_errors,
            kernel_name=self.kernel_name,
            resources={"metadata": {"path": path}},
        )

        # Run it
        client.execute()

        # Save it
        if output_path:
            self.log.info(f"Save executed results to {output_path}")
            nbformat.write(nb, output_path)


main = NbClientApp.launch_instance
