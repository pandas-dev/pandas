"""A qt exporter."""

import os
import sys
import tempfile

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
        from .qt_screenshot import QT_INSTALLED

        if not QT_INSTALLED:
            msg = (
                f"PyQtWebEngine is not installed to support Qt {self.format.upper()} conversion. "
                f"Please install `nbconvert[qt{self.format}]` to enable."
            )
            raise RuntimeError(msg)
        from .qt_screenshot import QtScreenshot

        return QtScreenshot

    def _run_pyqtwebengine(self, html):
        ext = ".html"
        temp_file = tempfile.NamedTemporaryFile(suffix=ext, delete=False)
        filename = f"{temp_file.name[:-len(ext)]}.{self.format}"
        with temp_file:
            temp_file.write(html.encode("utf-8"))
        try:
            QtScreenshot = self._check_launch_reqs()
            s = QtScreenshot()
            s.capture(f"file://{temp_file.name}", filename, self.paginate)
        finally:
            # Ensure the file is deleted even if pyqtwebengine raises an exception
            os.unlink(temp_file.name)
        return s.data

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
