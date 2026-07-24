# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""A license reporting CLI

Mostly ready-to-use, the downstream must provide the location of the application's
static resources. Licenses from an app's federated_extensions will also be discovered
as configured in `labextensions_path` and `extra_labextensions_path`.

    from traitlets import default
    from jupyterlab_server import LicensesApp

    class MyLicensesApp(LicensesApp):
        version = "0.1.0"

        @default("static_dir")
        def _default_static_dir(self):
            return "my-static/"

    class MyApp(JupyterApp, LabConfig):
        ...
        subcommands = dict(
            licenses=(MyLicensesApp, MyLicensesApp.description.splitlines()[0])
        )

"""
from typing import Any

from jupyter_core.application import JupyterApp, base_aliases, base_flags
from traitlets import Bool, Enum, Instance, Unicode

from ._version import __version__
from .config import LabConfig
from .licenses_handler import LicensesManager


class LicensesApp(JupyterApp, LabConfig):
    """A license management app."""

    version = __version__

    description = """
    Report frontend licenses
    """

    static_dir = Unicode("", config=True, help="The static directory from which to show licenses")

    full_text = Bool(False, config=True, help="Also print out full license text (if available)")

    report_format = Enum(
        ["markdown", "json", "csv"], "markdown", config=True, help="Reporter format"
    )

    bundles_pattern = Unicode(".*", config=True, help="A regular expression of bundles to print")

    licenses_manager = Instance(LicensesManager)

    aliases = {
        **base_aliases,
        "bundles": "LicensesApp.bundles_pattern",
        "report-format": "LicensesApp.report_format",
    }

    flags = {
        **base_flags,
        "full-text": (
            {"LicensesApp": {"full_text": True}},
            "Print out full license text (if available)",
        ),
        "json": (
            {"LicensesApp": {"report_format": "json"}},
            "Print out report as JSON (implies --full-text)",
        ),
        "csv": (
            {"LicensesApp": {"report_format": "csv"}},
            "Print out report as CSV (implies --full-text)",
        ),
    }

    def initialize(self, *args: Any, **kwargs: Any) -> None:
        """Initialize the app."""
        super().initialize(*args, **kwargs)
        self.init_licenses_manager()

    def init_licenses_manager(self) -> None:
        """Initialize the license manager."""
        self.licenses_manager = LicensesManager(
            parent=self,
        )

    def start(self) -> None:
        """Start the app."""
        report = self.licenses_manager.report(
            report_format=self.report_format,
            full_text=self.full_text,
            bundles_pattern=self.bundles_pattern,
        )[0]
        print(report)
        self.exit(0)
