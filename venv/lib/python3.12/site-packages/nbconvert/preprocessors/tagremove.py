"""
Module containing a preprocessor that removes cells if they match
one or more regular expression.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from traitlets import Set, Unicode

from .base import Preprocessor


class TagRemovePreprocessor(Preprocessor):
    """
    Removes inputs, outputs, or cells from a notebook that
    have tags that designate they are to be removed prior to exporting
    the notebook.

    remove_cell_tags
        removes cells tagged with these values

    remove_all_outputs_tags
        removes entire output areas on cells
        tagged with these values

    remove_single_output_tags
        removes individual output objects on
        outputs tagged with these values

    remove_input_tags
        removes inputs tagged with these values
    """

    remove_cell_tags: set[str] = Set(  # type:ignore[assignment]
        Unicode(),
        default_value=[],
        help=(
            "Tags indicating which cells are to be removed,matches tags in ``cell.metadata.tags``."
        ),
    ).tag(config=True)
    remove_all_outputs_tags: set[str] = Set(  # type:ignore[assignment]
        Unicode(),
        default_value=[],
        help=(
            "Tags indicating cells for which the outputs are to be removed,"
            "matches tags in ``cell.metadata.tags``."
        ),
    ).tag(config=True)
    remove_single_output_tags: set[str] = Set(  # type:ignore[assignment]
        Unicode(),
        default_value=[],
        help=(
            "Tags indicating which individual outputs are to be removed,"
            "matches output *i* tags in ``cell.outputs[i].metadata.tags``."
        ),
    ).tag(config=True)
    remove_input_tags: set[str] = Set(  # type:ignore[assignment]
        Unicode(),
        default_value=[],
        help=(
            "Tags indicating cells for which input is to be removed,"
            "matches tags in ``cell.metadata.tags``."
        ),
    ).tag(config=True)
    remove_metadata_fields: set[str] = Set({"collapsed", "scrolled"}).tag(config=True)  # type:ignore[assignment]

    def check_cell_conditions(self, cell, resources, index):
        """
        Checks that a cell has a tag that is to be removed

        Returns: Boolean.
        True means cell should *not* be removed.
        """

        # Return true if any of the tags in the cell are removable.
        return not self.remove_cell_tags.intersection(cell.get("metadata", {}).get("tags", []))

    def preprocess(self, nb, resources):
        """
        Preprocessing to apply to each notebook. See base.py for details.
        """
        # Skip preprocessing if the list of patterns is empty
        if not any(
            [
                self.remove_cell_tags,
                self.remove_all_outputs_tags,
                self.remove_single_output_tags,
                self.remove_input_tags,
            ]
        ):
            return nb, resources

        # Filter out cells that meet the conditions
        nb.cells = [
            self.preprocess_cell(cell, resources, index)[0]
            for index, cell in enumerate(nb.cells)
            if self.check_cell_conditions(cell, resources, index)
        ]

        return nb, resources

    def preprocess_cell(self, cell, resources, cell_index):
        """
        Apply a transformation on each cell. See base.py for details.
        """

        if (
            self.remove_all_outputs_tags.intersection(cell.get("metadata", {}).get("tags", []))
            and cell.cell_type == "code"
        ):
            cell.outputs = []
            cell.execution_count = None
            # Remove metadata associated with output
            if "metadata" in cell:
                for field in self.remove_metadata_fields:
                    cell.metadata.pop(field, None)

        if self.remove_input_tags.intersection(cell.get("metadata", {}).get("tags", [])):
            cell.metadata["transient"] = {"remove_source": True}

        if cell.get("outputs", []):
            cell.outputs = [
                output
                for output_index, output in enumerate(cell.outputs)
                if self.check_output_conditions(output, resources, cell_index, output_index)
            ]
        return cell, resources

    def check_output_conditions(self, output, resources, cell_index, output_index):
        """
        Checks that an output has a tag that indicates removal.

        Returns: Boolean.
        True means output should *not* be removed.
        """
        return not self.remove_single_output_tags.intersection(
            output.get("metadata", {}).get("tags", [])
        )
