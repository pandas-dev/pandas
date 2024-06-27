"""Contains writer for writing nbconvert output to filesystem."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import errno
import glob
import os
from pathlib import Path

from traitlets import Unicode, observe

from nbconvert.utils.io import link_or_copy

from .base import WriterBase


class FilesWriter(WriterBase):
    """Consumes nbconvert output and produces files."""

    build_directory = Unicode(
        "",
        help="""Directory to write output(s) to. Defaults
                              to output to the directory of each notebook. To recover
                              previous default behaviour (outputting to the current
                              working directory) use . as the flag value.""",
    ).tag(config=True)

    relpath = Unicode(
        help="""When copying files that the notebook depends on, copy them in
        relation to this path, such that the destination filename will be
        os.path.relpath(filename, relpath). If FilesWriter is operating on a
        notebook that already exists elsewhere on disk, then the default will be
        the directory containing that notebook."""
    ).tag(config=True)

    # Make sure that the output directory exists.
    @observe("build_directory")
    def _build_directory_changed(self, change):
        new = change["new"]
        if new:
            self._makedir(new)

    def __init__(self, **kw):
        """Initialize the writer."""
        super().__init__(**kw)
        self._build_directory_changed({"new": self.build_directory})

    def _makedir(self, path, mode=0o755):
        """ensure that a directory exists

        If it doesn't exist, try to create it and protect against a race condition
        if another process is doing the same.

        The default permissions are 755, which differ from os.makedirs default of 777.
        """
        if not os.path.exists(path):
            self.log.info("Making directory %s", path)
            try:
                os.makedirs(path, mode=mode)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        elif not os.path.isdir(path):
            raise OSError("%r exists but is not a directory" % path)

    def _write_items(self, items, build_dir):
        """Write a dict containing filename->binary data"""
        for filename, data in items:
            # Determine where to write the file to
            dest = os.path.join(build_dir, filename)
            path = os.path.dirname(dest)
            self._makedir(path)

            # Write file
            self.log.debug("Writing %i bytes to %s", len(data), dest)
            with open(dest, "wb") as f:
                f.write(data)

    def write(self, output, resources, notebook_name=None, **kw):
        """
        Consume and write Jinja output to the file system.  Output directory
        is set via the 'build_directory' variable of this instance (a
        configurable).

        See base for more...
        """

        # Verify that a notebook name is provided.
        if notebook_name is None:
            msg = "notebook_name"
            raise TypeError(msg)

        # Pull the extension and subdir from the resources dict.
        output_extension = resources.get("output_extension", None)

        # Get the relative path for copying files
        resource_path = resources.get("metadata", {}).get("path", "")
        relpath = self.relpath or resource_path
        build_directory = self.build_directory or resource_path

        # Write the extracted outputs to the destination directory.
        # NOTE: WE WRITE EVERYTHING AS-IF IT'S BINARY.  THE EXTRACT FIG
        # PREPROCESSOR SHOULD HANDLE UNIX/WINDOWS LINE ENDINGS...

        items = resources.get("outputs", {}).items()
        if items:
            self.log.info(
                "Support files will be in %s",
                os.path.join(resources.get("output_files_dir", ""), ""),
            )
            self._write_items(items, build_directory)

        # Write the extracted attachments
        # if ExtractAttachmentsOutput specified a separate directory
        attachments = resources.get("attachments", {}).items()
        if attachments:
            self.log.info(
                "Attachments will be in %s",
                os.path.join(resources.get("attachment_files_dir", ""), ""),
            )
            self._write_items(attachments, build_directory)

        # Copy referenced files to output directory
        if build_directory:
            for filename in self.files:
                # Copy files that match search pattern
                for matching_filename in glob.glob(filename):
                    # compute the relative path for the filename
                    if relpath != "":
                        dest_filename = os.path.relpath(matching_filename, relpath)
                    else:
                        dest_filename = matching_filename

                    # Make sure folder exists.
                    dest = os.path.join(build_directory, dest_filename)
                    path = os.path.dirname(dest)
                    self._makedir(path)

                    # Copy if destination is different.
                    if os.path.normpath(dest) != os.path.normpath(matching_filename):
                        self.log.info("Copying %s -> %s", matching_filename, dest)
                        link_or_copy(matching_filename, dest)

        # Determine where to write conversion results.
        dest = notebook_name + output_extension if output_extension is not None else notebook_name
        dest_path = Path(build_directory) / dest

        # Write conversion results.
        self.log.info("Writing %i bytes to %s", len(output), dest_path)
        if isinstance(output, str):
            with open(dest_path, "w", encoding="utf-8") as f:
                f.write(output)
        else:
            with open(dest_path, "wb") as f:
                f.write(output)

        return dest_path
