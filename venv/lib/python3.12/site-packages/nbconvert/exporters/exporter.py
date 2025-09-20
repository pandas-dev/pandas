"""This module defines a base Exporter class. For Jinja template-based export,
see templateexporter.py.
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import collections
import copy
import datetime
import os
import sys
import typing as t

import nbformat
from nbformat import NotebookNode, validator
from traitlets import Bool, HasTraits, List, TraitError, Unicode
from traitlets.config import Config
from traitlets.config.configurable import LoggingConfigurable
from traitlets.utils.importstring import import_item


class ResourcesDict(collections.defaultdict):  # type:ignore[type-arg]
    """A default dict for resources."""

    def __missing__(self, key):
        """Handle missing value."""
        return ""


class FilenameExtension(Unicode):  # type:ignore[type-arg]
    """A trait for filename extensions."""

    default_value = ""
    info_text = "a filename extension, beginning with a dot"

    def validate(self, obj, value):
        """Validate the file name."""
        # cast to proper unicode
        value = super().validate(obj, value)

        # check that it starts with a dot
        if value and not value.startswith("."):
            msg = "FileExtension trait '{}' does not begin with a dot: {!r}"
            raise TraitError(msg.format(self.name, value))

        return value


class Exporter(LoggingConfigurable):
    """
    Class containing methods that sequentially run a list of preprocessors on a
    NotebookNode object and then return the modified NotebookNode object and
    accompanying resources dict.
    """

    enabled = Bool(True, help="Disable this exporter (and any exporters inherited from it).").tag(
        config=True
    )

    file_extension = FilenameExtension(
        help="Extension of the file that should be written to disk"
    ).tag(config=True)

    optimistic_validation = Bool(
        False,
        help="Reduces the number of validation steps so that it only occurs after all preprocesors have run.",
    ).tag(config=True)

    # MIME type of the result file, for HTTP response headers.
    # This is *not* a traitlet, because we want to be able to access it from
    # the class, not just on instances.
    output_mimetype = ""

    # Should this converter be accessible from the notebook front-end?
    # If so, should be a friendly name to display (and possibly translated).
    export_from_notebook: str = None  # type:ignore[assignment]

    # Configurability, allows the user to easily add filters and preprocessors.
    preprocessors: List[t.Any] = List(
        help="""List of preprocessors, by name or namespace, to enable."""
    ).tag(config=True)

    _preprocessors: List[t.Any] = List()

    default_preprocessors: List[t.Any] = List(
        [
            "nbconvert.preprocessors.TagRemovePreprocessor",
            "nbconvert.preprocessors.RegexRemovePreprocessor",
            "nbconvert.preprocessors.ClearOutputPreprocessor",
            "nbconvert.preprocessors.CoalesceStreamsPreprocessor",
            "nbconvert.preprocessors.ExecutePreprocessor",
            "nbconvert.preprocessors.SVG2PDFPreprocessor",
            "nbconvert.preprocessors.LatexPreprocessor",
            "nbconvert.preprocessors.HighlightMagicsPreprocessor",
            "nbconvert.preprocessors.ExtractOutputPreprocessor",
            "nbconvert.preprocessors.ExtractAttachmentsPreprocessor",
            "nbconvert.preprocessors.ClearMetadataPreprocessor",
        ],
        help="""List of preprocessors available by default, by name, namespace,
        instance, or type.""",
    ).tag(config=True)

    def __init__(self, config=None, **kw):
        """
        Public constructor

        Parameters
        ----------
        config : ``traitlets.config.Config``
            User configuration instance.
        `**kw`
            Additional keyword arguments passed to parent __init__

        """
        with_default_config = self.default_config
        if config:
            with_default_config.merge(config)

        super().__init__(config=with_default_config, **kw)

        self._init_preprocessors()
        self._nb_metadata = {}

    @property
    def default_config(self):
        return Config()

    def from_notebook_node(
        self, nb: NotebookNode, resources: t.Any | None = None, **kw: t.Any
    ) -> tuple[NotebookNode, dict[str, t.Any]]:
        """
        Convert a notebook from a notebook node instance.

        Parameters
        ----------
        nb : :class:`~nbformat.NotebookNode`
            Notebook node (dict-like with attr-access)
        resources : dict
            Additional resources that can be accessed read/write by
            preprocessors and filters.
        `**kw`
            Ignored

        """
        nb_copy = copy.deepcopy(nb)
        resources = self._init_resources(resources)

        if "language" in nb["metadata"]:
            resources["language"] = nb["metadata"]["language"].lower()

        # Preprocess
        nb_copy, resources = self._preprocess(nb_copy, resources)
        notebook_name = ""
        if resources is not None:
            name = resources.get("metadata", {}).get("name", "")
            path = resources.get("metadata", {}).get("path", "")
            notebook_name = os.path.join(path, name)
        self._nb_metadata[notebook_name] = nb_copy.metadata
        return nb_copy, resources

    def from_filename(
        self, filename: str, resources: dict[str, t.Any] | None = None, **kw: t.Any
    ) -> tuple[NotebookNode, dict[str, t.Any]]:
        """
        Convert a notebook from a notebook file.

        Parameters
        ----------
        filename : str
            Full filename of the notebook file to open and convert.
        resources : dict
            Additional resources that can be accessed read/write by
            preprocessors and filters.
        `**kw`
            Ignored

        """
        # Pull the metadata from the filesystem.
        if resources is None:
            resources = ResourcesDict()
        if "metadata" not in resources or resources["metadata"] == "":
            resources["metadata"] = ResourcesDict()
        path, basename = os.path.split(filename)
        notebook_name = os.path.splitext(basename)[0]
        resources["metadata"]["name"] = notebook_name
        resources["metadata"]["path"] = path

        modified_date = datetime.datetime.fromtimestamp(
            os.path.getmtime(filename), tz=datetime.timezone.utc
        )
        # datetime.strftime date format for ipython
        if sys.platform == "win32":
            date_format = "%B %d, %Y"
        else:
            date_format = "%B %-d, %Y"
        resources["metadata"]["modified_date"] = modified_date.strftime(date_format)

        with open(filename, encoding="utf-8") as f:
            return self.from_file(f, resources=resources, **kw)

    def from_file(
        self, file_stream: t.Any, resources: dict[str, t.Any] | None = None, **kw: t.Any
    ) -> tuple[NotebookNode, dict[str, t.Any]]:
        """
        Convert a notebook from a notebook file.

        Parameters
        ----------
        file_stream : file-like object
            Notebook file-like object to convert.
        resources : dict
            Additional resources that can be accessed read/write by
            preprocessors and filters.
        `**kw`
            Ignored

        """
        return self.from_notebook_node(
            nbformat.read(file_stream, as_version=4), resources=resources, **kw
        )

    def register_preprocessor(self, preprocessor, enabled=False):
        """
        Register a preprocessor.
        Preprocessors are classes that act upon the notebook before it is
        passed into the Jinja templating engine. Preprocessors are also
        capable of passing additional information to the Jinja
        templating engine.

        Parameters
        ----------
        preprocessor : `nbconvert.preprocessors.Preprocessor`
            A dotted module name, a type, or an instance
        enabled : bool
            Mark the preprocessor as enabled

        """
        if preprocessor is None:
            msg = "preprocessor must not be None"
            raise TypeError(msg)
        isclass = isinstance(preprocessor, type)
        constructed = not isclass

        # Handle preprocessor's registration based on it's type
        if constructed and isinstance(
            preprocessor,
            str,
        ):
            # Preprocessor is a string, import the namespace and recursively call
            # this register_preprocessor method
            preprocessor_cls = import_item(preprocessor)
            return self.register_preprocessor(preprocessor_cls, enabled)

        if constructed and callable(preprocessor):
            # Preprocessor is a function, no need to construct it.
            # Register and return the preprocessor.
            if enabled:
                preprocessor.enabled = True
            self._preprocessors.append(preprocessor)
            return preprocessor

        if isclass and issubclass(preprocessor, HasTraits):
            # Preprocessor is configurable.  Make sure to pass in new default for
            # the enabled flag if one was specified.
            self.register_preprocessor(preprocessor(parent=self), enabled)
            return None

        if isclass:
            # Preprocessor is not configurable, construct it
            self.register_preprocessor(preprocessor(), enabled)
            return None

        # Preprocessor is an instance of something without a __call__
        # attribute.
        raise TypeError(
            "preprocessor must be callable or an importable constructor, got %r" % preprocessor
        )

    def _init_preprocessors(self):
        """
        Register all of the preprocessors needed for this exporter, disabled
        unless specified explicitly.
        """
        self._preprocessors = []

        # Load default preprocessors (not necessarily enabled by default).
        for preprocessor in self.default_preprocessors:
            self.register_preprocessor(preprocessor)

        # Load user-specified preprocessors.  Enable by default.
        for preprocessor in self.preprocessors:
            self.register_preprocessor(preprocessor, enabled=True)

    def _init_resources(self, resources):
        # Make sure the resources dict is of ResourcesDict type.
        if resources is None:
            resources = ResourcesDict()
        if not isinstance(resources, ResourcesDict):
            new_resources = ResourcesDict()
            new_resources.update(resources)
            resources = new_resources

        # Make sure the metadata extension exists in resources
        if "metadata" in resources:
            if not isinstance(resources["metadata"], ResourcesDict):
                new_metadata = ResourcesDict()
                new_metadata.update(resources["metadata"])
                resources["metadata"] = new_metadata
        else:
            resources["metadata"] = ResourcesDict()
            if not resources["metadata"]["name"]:
                resources["metadata"]["name"] = "Notebook"

        # Set the output extension
        resources["output_extension"] = self.file_extension
        return resources

    def _validate_preprocessor(self, nbc, preprocessor):
        try:
            nbformat.validate(nbc, relax_add_props=True)
        except nbformat.ValidationError:
            self.log.error("Notebook is invalid after preprocessor %s", preprocessor)
            raise

    def _preprocess(self, nb, resources):
        """
        Preprocess the notebook before passing it into the Jinja engine.
        To preprocess the notebook is to successively apply all the
        enabled preprocessors. Output from each preprocessor is passed
        along to the next one.

        Parameters
        ----------
        nb : notebook node
            notebook that is being exported.
        resources : a dict of additional resources that
            can be accessed read/write by preprocessors
        """

        # Do a copy.deepcopy first,
        # we are never safe enough with what the preprocessors could do.
        nbc = copy.deepcopy(nb)
        resc = copy.deepcopy(resources)

        if hasattr(validator, "normalize"):
            _, nbc = validator.normalize(nbc)

        # Run each preprocessor on the notebook.  Carry the output along
        # to each preprocessor
        for preprocessor in self._preprocessors:
            nbc, resc = preprocessor(nbc, resc)
            if not self.optimistic_validation:
                self._validate_preprocessor(nbc, preprocessor)

        if self.optimistic_validation:
            self._validate_preprocessor(nbc, preprocessor)

        return nbc, resc
