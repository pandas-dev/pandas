"""A qt exporter."""

import os
import sys
import tempfile
import time

from traitlets import default

from .html import HTMLExporter


class QtExporter(HTMLExporter):
    """A qt exporter."""

    paginate = None
    format = ""

    @default("file_extension")
    def _file_extension_default(self):
        return ".html"

    def _check_launch_reqs(self):
        if sys.platform.startswith("win") and self.format == "png":
            msg = "Exporting to PNG using Qt is currently not supported on Windows."
            raise RuntimeError(msg)
        from .qt_screenshot import QT_INSTALLED  # noqa: PLC0415

        if not QT_INSTALLED:
            msg = (
                f"PyQtWebEngine is not installed to support Qt {self.format.upper()} conversion. "
                f"Please install `nbconvert[qt{self.format}]` to enable."
            )
            raise RuntimeError(msg)
        from .qt_screenshot import QtScreenshot  # noqa: PLC0415

        return QtScreenshot

    def _run_pyqtwebengine(self, html):
        ext = ".html"
        temp_file = tempfile.NamedTemporaryFile(  # noqa: SIM115
            suffix=ext, delete=False
        )
        filename = f"{temp_file.name[: -len(ext)]}.{self.format}"
        with temp_file:
            temp_file.write(html.encode("utf-8"))
        try:
            QtScreenshot = self._check_launch_reqs()
            s = QtScreenshot()
            s.capture(f"file://{temp_file.name}", filename, self.paginate)
        finally:
            # Ensure the file is deleted even if pyqtwebengine raises an exception
            os.unlink(temp_file.name)
        # Prefer Qt's in-memory bytes, but fall back to reading the file on disk
        data = getattr(s, "data", b"")

        if (not data) and os.path.exists(filename):
            deadline = time.time() + 5.0
            while time.time() < deadline:
                try:
                    if os.path.getsize(filename) > 0:
                        break
                except OSError:
                    pass
                time.sleep(0.05)

            if os.path.exists(filename) and os.path.getsize(filename) > 0:
                with open(filename, "rb") as f:
                    data = f.read()

        # Best-effort cleanup of the generated output file
        try:
            if os.path.exists(filename):
                os.unlink(filename)
        except OSError:
            pass

        return data

    def from_notebook_node(self, nb, resources=None, **kw):
        """Convert from notebook node."""
        self._check_launch_reqs()
        html, resources = super().from_notebook_node(nb, resources=resources, **kw)

        self.log.info("Building %s", self.format.upper())
        data = self._run_pyqtwebengine(html)
        self.log.info("%s successfully created", self.format.upper())

        # convert output extension
        # the writer above required it to be html
        resources["output_extension"] = f".{self.format}"

        return data, resources
