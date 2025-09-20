"""Python script Exporter class"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import default

from .templateexporter import TemplateExporter


class PythonExporter(TemplateExporter):
    """
    Exports a Python code file.
    Note that the file produced will have a shebang of '#!/usr/bin/env python'
    regardless of the actual python version used in the notebook.
    """

    @default("file_extension")
    def _file_extension_default(self):
        return ".py"

    @default("template_name")
    def _template_name_default(self):
        return "python"

    output_mimetype = "text/x-python"
