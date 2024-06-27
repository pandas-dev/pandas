"""Module containing a preprocessor that removes the outputs from code cells"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import Set

from .base import Preprocessor


class ClearOutputPreprocessor(Preprocessor):
    """
    Removes the output from all code cells in a notebook.
    """

    remove_metadata_fields = Set({"collapsed", "scrolled"}).tag(config=True)

    def preprocess_cell(self, cell, resources, cell_index):
        """
        Apply a transformation on each cell. See base.py for details.
        """
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
            # Remove metadata associated with output
            if "metadata" in cell:
                for field in self.remove_metadata_fields:
                    cell.metadata.pop(field, None)
        return cell, resources
