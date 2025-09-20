"""NotebookExporter class"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import nbformat
from traitlets import Enum, default

from .exporter import Exporter


class NotebookExporter(Exporter):
    """Exports to an IPython notebook.

    This is useful when you want to use nbconvert's preprocessors to operate on
    a notebook (e.g. to execute it) and then write it back to a notebook file.
    """

    nbformat_version = Enum(
        list(nbformat.versions),
        default_value=nbformat.current_nbformat,
        help="""The nbformat version to write.
        Use this to downgrade notebooks.
        """,
    ).tag(config=True)

    @default("file_extension")
    def _file_extension_default(self):
        return ".ipynb"

    output_mimetype = "application/json"
    export_from_notebook = "Notebook"

    def from_notebook_node(self, nb, resources=None, **kw):
        """Convert from notebook node."""
        nb_copy, resources = super().from_notebook_node(nb, resources, **kw)
        if self.nbformat_version != nb_copy.nbformat:
            resources["output_suffix"] = ".v%i" % self.nbformat_version
        else:
            resources["output_suffix"] = ".nbconvert"
        output = nbformat.writes(nb_copy, version=self.nbformat_version)
        if not output.endswith("\n"):
            output = output + "\n"
        return output, resources
