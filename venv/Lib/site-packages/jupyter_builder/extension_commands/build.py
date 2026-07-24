"""Build command for JupyterLab extensions."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import annotations

from pathlib import Path

from traitlets import Bool, Unicode

from jupyter_builder.base_extension_app import BaseExtensionApp
from jupyter_builder.federated_extensions import build_labextension

HERE = str(Path(__file__).resolve().parent)


class BuildLabExtensionApp(BaseExtensionApp):
    """Application for building a prebuilt JupyterLab extension."""

    description = "Build labextension"

    static_url = Unicode("", config=True, help="Sets the url for static assets when building")

    development = Bool(False, config=True, help="Build in development mode")

    source_map = Bool(False, config=True, help="Generate source maps")

    core_package_file = Unicode(
        "",
        config=True,
        help="Path to the core application package definition file",
    )

    core_path = Unicode(
        "",
        config=True,
        help="(Deprecated) Directory containing core application package.json file",
    )

    core_version = Unicode(
        "",
        config=True,
        help=(
            "Version of JupyterLab core to use when building, e.g. 'X.Y.Z' or 'vX.Y.Z' "
            "(ignored if core-package-file is set). Building fails if the version cannot "
            "be resolved."
        ),
    )

    aliases = {  # noqa: RUF012
        "static-url": "BuildLabExtensionApp.static_url",
        "development": "BuildLabExtensionApp.development",
        "source-map": "BuildLabExtensionApp.source_map",
        "core-package-file": "BuildLabExtensionApp.core_package_file",
        "core-path": "BuildLabExtensionApp.core_path",
        "core-version": "BuildLabExtensionApp.core_version",
    }

    def run_task(self) -> None:
        """Build the labextension in the configured path."""
        self.extra_args = self.extra_args or [str(Path.cwd())]
        build_labextension(
            self.extra_args[0],
            logger=self.log,
            development=self.development,
            static_url=self.static_url or None,
            source_map=self.source_map,
            core_version=self.core_version or None,
            core_package_file=self.core_package_file or None,
            core_path=self.core_path or None,
        )


def main() -> None:
    """Run the build labextension app."""
    app = BuildLabExtensionApp()
    app.initialize()
    app.start()


if __name__ == "__main__":
    main()
