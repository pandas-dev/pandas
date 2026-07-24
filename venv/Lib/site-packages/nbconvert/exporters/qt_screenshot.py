"""A qt screenshot exporter."""

import os

try:
    from PyQt5 import QtCore  # type:ignore[import-not-found]
    from PyQt5.QtGui import QPageLayout, QPageSize  # type:ignore[import-not-found]
    from PyQt5.QtWebEngineWidgets import (  # type:ignore[import-not-found]
        QWebEngineSettings,
        QWebEngineView,
    )
    from PyQt5.QtWidgets import QApplication  # type:ignore[import-not-found]

    QT_INSTALLED = True
except ModuleNotFoundError:
    QT_INSTALLED = False


if QT_INSTALLED:
    APP = None
    if not QApplication.instance():
        APP = QApplication([])

    class QtScreenshot(QWebEngineView):  # type:ignore[misc]
        """A qt screenshot exporter."""

        def __init__(self):
            """Initialize the exporter."""
            super().__init__()
            self.app = APP

        def capture(self, url, output_file, paginate):
            """Capture the screenshot."""
            self.output_file = output_file
            self.paginate = paginate
            self.load(QtCore.QUrl(url))
            self.loadFinished.connect(self.on_loaded)
            # Create hidden view without scrollbars
            self.setAttribute(QtCore.Qt.WA_DontShowOnScreen)
            self.page().settings().setAttribute(QWebEngineSettings.ShowScrollBars, False)
            self.data = b""
            if output_file.endswith(".pdf"):
                self.export = self.export_pdf

                def cleanup(*args):
                    """Cleanup the app."""
                    self.get_data()
                    self.app.quit()  # type:ignore[union-attr]

                self.page().pdfPrintingFinished.connect(cleanup)
            elif output_file.endswith(".png"):
                self.export = self.export_png
            else:
                msg = f"Export file extension not supported: {output_file}"
                raise RuntimeError(msg)
            self.show()
            self.app.exec()  # type:ignore[union-attr]

        def on_loaded(self):
            """Handle app load."""
            self.size = self.page().contentsSize().toSize()
            self.resize(self.size)
            # Wait for resize
            QtCore.QTimer.singleShot(1000, self.export)

        def export_pdf(self):
            """Export to pdf."""
            if self.paginate:
                page_size = QPageSize(QPageSize.A4)
                page_layout = QPageLayout(page_size, QPageLayout.Portrait, QtCore.QMarginsF())
            else:
                factor = 0.75
                page_size = QPageSize(
                    QtCore.QSizeF(self.size.width() * factor, self.size.height() * factor),
                    QPageSize.Point,
                )
                page_layout = QPageLayout(page_size, QPageLayout.Portrait, QtCore.QMarginsF())

            self.page().printToPdf(self.output_file, pageLayout=page_layout)

        def export_png(self):
            """Export to png."""
            self.grab().save(self.output_file, "PNG")
            self.get_data()
            self.app.quit()  # type:ignore[union-attr]

        def get_data(self):
            """Get output data."""
            if os.path.exists(self.output_file):
                with open(self.output_file, "rb") as f:
                    self.data = f.read()
                os.unlink(self.output_file)
