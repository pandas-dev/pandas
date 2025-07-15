"""ASCIIDoc Exporter class"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import default
from traitlets.config import Config

from .templateexporter import TemplateExporter


class ASCIIDocExporter(TemplateExporter):
    """
    Exports to an ASCIIDoc document (.asciidoc)
    """

    @default("file_extension")
    def _file_extension_default(self):
        return ".asciidoc"

    @default("template_name")
    def _template_name_default(self):
        return "asciidoc"

    output_mimetype = "text/asciidoc"
    export_from_notebook = "AsciiDoc"

    @default("raw_mimetypes")
    def _raw_mimetypes_default(self):
        return ["text/asciidoc/", "text/markdown", "text/html", ""]

    @property
    def default_config(self):
        c = Config(
            {
                "NbConvertBase": {
                    "display_data_priority": [
                        "text/html",
                        "text/markdown",
                        "image/svg+xml",
                        "image/png",
                        "image/jpeg",
                        "text/plain",
                        "text/latex",
                    ]
                },
                "ExtractOutputPreprocessor": {"enabled": True},
                "HighlightMagicsPreprocessor": {"enabled": True},
            }
        )
        if super().default_config:
            c2 = super().default_config.copy()
            c2.merge(c)
            c = c2
        return c
