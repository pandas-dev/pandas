"""
Module that extracts attachments from notebooks into their own files
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import os
from base64 import b64decode

from traitlets import Bool, Unicode

from .base import Preprocessor


class ExtractAttachmentsPreprocessor(Preprocessor):
    """
    Extracts attachments from all (markdown and raw) cells in a notebook.
    The extracted attachments are stored in a directory ('attachments' by default).
    https://nbformat.readthedocs.io/en/latest/format_description.html#cell-attachments
    """

    attachments_directory_template = Unicode(
        "{notebook_name}_attachments",
        help="Directory to place attachments if use_separate_dir is True",
    ).tag(config=True)

    use_separate_dir = Bool(
        False,
        help="Whether to use output_files_dir (which ExtractOutput also uses) or "
        "create a separate directory for attachments",
    ).tag(config=True)

    def __init__(self, **kw):
        """
        Public constructor
        """
        super().__init__(**kw)
        # directory path,
        self.path_name = ""  # will be set in self.preprocess, needs resources
        # Where extracted attachments are stored in resources
        self.resources_item_key = (
            "attachments"  # Here as a default, in case someone doesn't want to call preprocess
        )

    # Add condition and configurability here
    def preprocess(self, nb, resources):
        """
        Determine some settings and apply preprocessor to notebook
        """
        if self.use_separate_dir:
            self.path_name = self.attachments_directory_template.format(
                notebook_name=resources["unique_key"]
            )
            # Initialize resources for attachments
            resources["attachment_files_dir"] = self.path_name
            resources["attachments"] = {}
            self.resources_item_key = "attachments"
        else:
            # Use same resources as ExtractOutput
            self.path_name = resources["output_files_dir"]
            self.resources_item_key = "outputs"

        # Make sure key exists
        if not isinstance(resources[self.resources_item_key], dict):
            resources[self.resources_item_key] = {}

        nb, resources = super().preprocess(nb, resources)
        return nb, resources

    def preprocess_cell(self, cell, resources, index):
        """
        Extract attachments to individual files and
        change references to them.
        E.g.
        '![image.png](attachment:021fdd80.png)'
        becomes
        '![image.png]({path_name}/021fdd80.png)'
        Assumes self.path_name and self.resources_item_key is set properly (usually in preprocess).
        """
        if "attachments" in cell:
            for fname in cell.attachments:
                self.log.debug("Encountered attachment %s", fname)

                # Add file for writer

                # Right now I don't know of a situation where there would be multiple
                # mime types under same filename, and I can't index into it without the mimetype.
                # So I only read the first one.
                for mimetype in cell.attachments[fname]:
                    # convert to bytes and decode
                    data = cell.attachments[fname][mimetype].encode("utf-8")
                    decoded = b64decode(data)
                    break

                # FilesWriter wants path to be in attachment filename here
                new_filename = os.path.join(self.path_name, fname)
                resources[self.resources_item_key][new_filename] = decoded

                # Edit the reference to the attachment

                # os.path.join on windows uses "\\" separator,
                # but files like markdown still want "/"
                if os.path.sep != "/":
                    new_filename = new_filename.replace(os.path.sep, "/")
                cell.source = cell.source.replace("attachment:" + fname, new_filename)

        return cell, resources
