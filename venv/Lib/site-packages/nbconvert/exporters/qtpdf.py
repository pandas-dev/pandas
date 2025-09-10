"""Export to PDF via a headless browser"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import Bool

from .qt_exporter import QtExporter


class QtPDFExporter(QtExporter):
    """Writer designed to write to PDF files.

    This inherits from :class:`HTMLExporter`. It creates the HTML using the
    template machinery, and then uses pyqtwebengine to create a pdf.
    """

    export_from_notebook = "PDF via HTML"
    format = "pdf"

    paginate = Bool(  # type:ignore[assignment]
        True,
        help="""
        Split generated notebook into multiple pages.

        If False, a PDF with one long page will be generated.

        Set to True to match behavior of LaTeX based PDF generator
        """,
    ).tag(config=True)
