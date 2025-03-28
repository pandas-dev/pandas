"""Module containing a preprocessor that converts outputs in the notebook from
one format to another.

Converts all of the outputs in a notebook from one format to another.
"""
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import Unicode

from .base import Preprocessor


class ConvertFiguresPreprocessor(Preprocessor):
    """
    Converts all of the outputs in a notebook from one format to another.
    """

    from_format = Unicode(help="Format the converter accepts").tag(config=True)
    to_format = Unicode(help="Format the converter writes").tag(config=True)

    def __init__(self, **kw):
        """
        Public constructor
        """
        super().__init__(**kw)

    def convert_figure(self, data_format, data):
        """Convert the figure."""
        raise NotImplementedError()

    def preprocess_cell(self, cell, resources, cell_index):
        """
        Apply a transformation on each cell,

        See base.py
        """

        # Loop through all of the datatypes of the outputs in the cell.
        for output in cell.get("outputs", []):
            if (
                output.output_type in {"execute_result", "display_data"}
                and self.from_format in output.data
                and self.to_format not in output.data
            ):
                output.data[self.to_format] = self.convert_figure(
                    self.from_format, output.data[self.from_format]
                )

        return cell, resources
