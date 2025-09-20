"""Module containing a preprocessor that removes metadata from code cells"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import Bool, Set

from .base import Preprocessor


class ClearMetadataPreprocessor(Preprocessor):
    """
    Removes all the metadata from all code cells in a notebook.
    """

    clear_cell_metadata = Bool(
        True,
        help=("Flag to choose if cell metadata is to be cleared in addition to notebook metadata."),
    ).tag(config=True)
    clear_notebook_metadata = Bool(
        True,
        help=("Flag to choose if notebook metadata is to be cleared in addition to cell metadata."),
    ).tag(config=True)
    preserve_nb_metadata_mask = Set(
        [("language_info", "name")],
        help=(
            "Indicates the key paths to preserve when deleting metadata "
            "across both cells and notebook metadata fields. Tuples of "
            "keys can be passed to preserved specific nested values"
        ),
    ).tag(config=True)
    preserve_cell_metadata_mask = Set(
        help=(
            "Indicates the key paths to preserve when deleting metadata "
            "across both cells and notebook metadata fields. Tuples of "
            "keys can be passed to preserved specific nested values"
        )
    ).tag(config=True)

    def current_key(self, mask_key):
        """Get the current key for a mask key."""
        if isinstance(mask_key, str):
            return mask_key
        if len(mask_key) == 0:
            # Safeguard
            return None
        return mask_key[0]

    def current_mask(self, mask):
        """Get the current mask for a mask."""
        return {self.current_key(k) for k in mask if self.current_key(k) is not None}

    def nested_masks(self, mask):
        """Get the nested masks for a mask."""
        return {
            self.current_key(k[0]): k[1:]
            for k in mask
            if k and not isinstance(k, str) and len(k) > 1
        }

    def nested_filter(self, items, mask):
        """Get the nested filter for items given a mask."""
        keep_current = self.current_mask(mask)
        keep_nested_lookup = self.nested_masks(mask)
        for k, v in items:
            keep_nested = keep_nested_lookup.get(k)
            if k in keep_current:
                if keep_nested is not None:
                    if isinstance(v, dict):
                        yield k, dict(self.nested_filter(v.items(), keep_nested))
                else:
                    yield k, v

    def preprocess_cell(self, cell, resources, cell_index):
        """
        All the code cells are returned with an empty metadata field.
        """
        if self.clear_cell_metadata and cell.cell_type == "code":  # noqa: SIM102
            # Remove metadata
            if "metadata" in cell:
                cell.metadata = dict(
                    self.nested_filter(cell.metadata.items(), self.preserve_cell_metadata_mask)
                )
        return cell, resources

    def preprocess(self, nb, resources):
        """
        Preprocessing to apply on each notebook.

        Must return modified nb, resources.

        Parameters
        ----------
        nb : NotebookNode
            Notebook being converted
        resources : dictionary
            Additional resources used in the conversion process.  Allows
            preprocessors to pass variables into the Jinja engine.
        """
        nb, resources = super().preprocess(nb, resources)
        if self.clear_notebook_metadata and "metadata" in nb:
            nb.metadata = dict(
                self.nested_filter(nb.metadata.items(), self.preserve_nb_metadata_mask)
            )
        return nb, resources
