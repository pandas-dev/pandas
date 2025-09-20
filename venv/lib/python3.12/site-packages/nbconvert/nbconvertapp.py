#!/usr/bin/env python
"""NbConvert is a utility for conversion of .ipynb files.

Command-line interface for the NbConvert conversion utility.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import asyncio
import glob
import logging
import os
import sys
import typing as t
from textwrap import dedent, fill

from jupyter_core.application import JupyterApp, base_aliases, base_flags
from traitlets import Bool, DottedObjectName, Instance, List, Type, Unicode, default, observe
from traitlets.config import Configurable, catch_config_error
from traitlets.utils.importstring import import_item

from nbconvert import __version__, exporters, postprocessors, preprocessors, writers
from nbconvert.utils.text import indent

from .exporters.base import get_export_names, get_exporter
from .utils.base import NbConvertBase
from .utils.exceptions import ConversionException
from .utils.io import unicode_stdin_stream

# -----------------------------------------------------------------------------
# Classes and functions
# -----------------------------------------------------------------------------


class DottedOrNone(DottedObjectName):
    """A string holding a valid dotted object name in Python, such as A.b3._c
    Also allows for None type.
    """

    default_value = ""

    def validate(self, obj, value):
        """Validate an input."""
        if value is not None and len(value) > 0:
            return super().validate(obj, value)
        return value


nbconvert_aliases = {}
nbconvert_aliases.update(base_aliases)
nbconvert_aliases.update(
    {
        "to": "NbConvertApp.export_format",
        "template": "TemplateExporter.template_name",
        "template-file": "TemplateExporter.template_file",
        "theme": "HTMLExporter.theme",
        "sanitize_html": "HTMLExporter.sanitize_html",
        "writer": "NbConvertApp.writer_class",
        "post": "NbConvertApp.postprocessor_class",
        "output": "NbConvertApp.output_base",
        "output-dir": "FilesWriter.build_directory",
        "reveal-prefix": "SlidesExporter.reveal_url_prefix",
        "nbformat": "NotebookExporter.nbformat_version",
    }
)

nbconvert_flags = {}
nbconvert_flags.update(base_flags)
nbconvert_flags.update(
    {
        "execute": (
            {"ExecutePreprocessor": {"enabled": True}},
            "Execute the notebook prior to export.",
        ),
        "allow-errors": (
            {"ExecutePreprocessor": {"allow_errors": True}},
            (
                "Continue notebook execution even if one of the cells throws "
                "an error and include the error message in the cell output "
                "(the default behaviour is to abort conversion). This flag "
                "is only relevant if '--execute' was specified, too."
            ),
        ),
        "stdin": (
            {
                "NbConvertApp": {
                    "from_stdin": True,
                }
            },
            "read a single notebook file from stdin. Write the resulting notebook with default basename 'notebook.*'",
        ),
        "stdout": (
            {"NbConvertApp": {"writer_class": "StdoutWriter"}},
            "Write notebook output to stdout instead of files.",
        ),
        "inplace": (
            {
                "NbConvertApp": {
                    "use_output_suffix": False,
                    "export_format": "notebook",
                },
                "FilesWriter": {"build_directory": ""},
            },
            """Run nbconvert in place, overwriting the existing notebook (only
        relevant when converting to notebook format)""",
        ),
        "clear-output": (
            {
                "NbConvertApp": {
                    "use_output_suffix": False,
                    "export_format": "notebook",
                },
                "FilesWriter": {"build_directory": ""},
                "ClearOutputPreprocessor": {"enabled": True},
            },
            """Clear output of current file and save in place,
        overwriting the existing notebook. """,
        ),
        "coalesce-streams": (
            {
                "NbConvertApp": {"use_output_suffix": False, "export_format": "notebook"},
                "FilesWriter": {"build_directory": ""},
                "CoalesceStreamsPreprocessor": {"enabled": True},
            },
            """Coalesce consecutive stdout and stderr outputs into one stream (within each cell).""",
        ),
        "no-prompt": (
            {
                "TemplateExporter": {
                    "exclude_input_prompt": True,
                    "exclude_output_prompt": True,
                }
            },
            "Exclude input and output prompts from converted document.",
        ),
        "no-input": (
            {
                "TemplateExporter": {
                    "exclude_output_prompt": True,
                    "exclude_input": True,
                    "exclude_input_prompt": True,
                }
            },
            """Exclude input cells and output prompts from converted document.
        This mode is ideal for generating code-free reports.""",
        ),
        "allow-chromium-download": (
            {
                "WebPDFExporter": {
                    "allow_chromium_download": True,
                }
            },
            """Whether to allow downloading chromium if no suitable version is found on the system.""",
        ),
        "disable-chromium-sandbox": (
            {
                "WebPDFExporter": {
                    "disable_sandbox": True,
                }
            },
            """Disable chromium security sandbox when converting to PDF..""",
        ),
        "show-input": (
            {
                "TemplateExporter": {
                    "exclude_input": False,
                }
            },
            """Shows code input. This flag is only useful for dejavu users.""",
        ),
        "embed-images": (
            {
                "HTMLExporter": {
                    "embed_images": True,
                }
            },
            """Embed the images as base64 dataurls in the output. This flag is only useful for the HTML/WebPDF/Slides exports.""",
        ),
        "sanitize-html": (
            {
                "HTMLExporter": {
                    "sanitize_html": True,
                }
            },
            """Whether the HTML in Markdown cells and cell outputs should be sanitized..""",
        ),
    }
)


class NbConvertApp(JupyterApp):
    """Application used to convert from notebook file type (``*.ipynb``)"""

    version = __version__
    name = "jupyter-nbconvert"
    aliases = nbconvert_aliases
    flags = nbconvert_flags

    @default("log_level")
    def _log_level_default(self):
        return logging.INFO

    classes = List()  # type:ignore[assignment]

    @default("classes")
    def _classes_default(self):
        classes: list[type[t.Any]] = [NbConvertBase]
        for pkg in (exporters, preprocessors, writers, postprocessors):
            for name in dir(pkg):
                cls = getattr(pkg, name)
                if isinstance(cls, type) and issubclass(cls, Configurable):
                    classes.append(cls)

        return classes

    description = Unicode(  # type:ignore[assignment]
        """This application is used to convert notebook files (*.ipynb)
        to various other formats.

        WARNING: THE COMMANDLINE INTERFACE MAY CHANGE IN FUTURE RELEASES."""
    )

    output_base = Unicode(
        "{notebook_name}",
        help="""Overwrite base name use for output files.
            Supports pattern replacements '{notebook_name}'.
            """,
    ).tag(config=True)

    use_output_suffix = Bool(
        True,
        help="""Whether to apply a suffix prior to the extension (only relevant
            when converting to notebook format). The suffix is determined by
            the exporter, and is usually '.nbconvert'.""",
    ).tag(config=True)

    output_files_dir = Unicode(
        "{notebook_name}_files",
        help="""Directory to copy extra files (figures) to.
               '{notebook_name}' in the string will be converted to notebook
               basename.""",
    ).tag(config=True)

    examples = Unicode(
        f"""
        The simplest way to use nbconvert is

        > jupyter nbconvert mynotebook.ipynb --to html

        Options include {get_export_names()}.

        > jupyter nbconvert --to latex mynotebook.ipynb

        Both HTML and LaTeX support multiple output templates. LaTeX includes
        'base', 'article' and 'report'.  HTML includes 'basic', 'lab' and
        'classic'. You can specify the flavor of the format used.

        > jupyter nbconvert --to html --template lab mynotebook.ipynb

        You can also pipe the output to stdout, rather than a file

        > jupyter nbconvert mynotebook.ipynb --stdout

        PDF is generated via latex

        > jupyter nbconvert mynotebook.ipynb --to pdf

        You can get (and serve) a Reveal.js-powered slideshow

        > jupyter nbconvert myslides.ipynb --to slides --post serve

        Multiple notebooks can be given at the command line in a couple of
        different ways:

        > jupyter nbconvert notebook*.ipynb
        > jupyter nbconvert notebook1.ipynb notebook2.ipynb

        or you can specify the notebooks list in a config file, containing::

            c.NbConvertApp.notebooks = ["my_notebook.ipynb"]

        > jupyter nbconvert --config mycfg.py
        """
    )

    # Writer specific variables
    writer = Instance(
        "nbconvert.writers.base.WriterBase",
        help="""Instance of the writer class used to write the
                      results of the conversion.""",
        allow_none=True,
    )
    writer_class = DottedObjectName(
        "FilesWriter",
        help="""Writer class used to write the
                                    results of the conversion""",
    ).tag(config=True)
    writer_aliases = {
        "fileswriter": "nbconvert.writers.files.FilesWriter",
        "debugwriter": "nbconvert.writers.debug.DebugWriter",
        "stdoutwriter": "nbconvert.writers.stdout.StdoutWriter",
    }
    writer_factory = Type(allow_none=True)

    @observe("writer_class")
    def _writer_class_changed(self, change):
        new = change["new"]
        if new.lower() in self.writer_aliases:
            new = self.writer_aliases[new.lower()]
        self.writer_factory = import_item(new)

    # Post-processor specific variables
    postprocessor = Instance(
        "nbconvert.postprocessors.base.PostProcessorBase",
        help="""Instance of the PostProcessor class used to write the
                      results of the conversion.""",
        allow_none=True,
    )

    postprocessor_class = DottedOrNone(
        help="""PostProcessor class used to write the
                                    results of the conversion"""
    ).tag(config=True)
    postprocessor_aliases = {"serve": "nbconvert.postprocessors.serve.ServePostProcessor"}
    postprocessor_factory = Type(None, allow_none=True)

    @observe("postprocessor_class")
    def _postprocessor_class_changed(self, change):
        new = change["new"]
        if new.lower() in self.postprocessor_aliases:
            new = self.postprocessor_aliases[new.lower()]
        if new:
            self.postprocessor_factory = import_item(new)

    export_format = Unicode(  # type:ignore[call-overload]
        allow_none=False,
        help=f"""The export format to be used, either one of the built-in formats
        {get_export_names()}
        or a dotted object name that represents the import path for an
        ``Exporter`` class""",
    ).tag(config=True)

    notebooks = List(
        Unicode(),
        help="""List of notebooks to convert.
                     Wildcards are supported.
                     Filenames passed positionally will be added to the list.
                     """,
    ).tag(config=True)
    from_stdin = Bool(False, help="read a single notebook from stdin.").tag(config=True)
    recursive_glob = Bool(
        False, help="set the 'recursive' option for glob for searching wildcards."
    ).tag(config=True)

    @catch_config_error
    def initialize(self, argv=None):
        """Initialize application, notebooks, writer, and postprocessor"""
        # See https://bugs.python.org/issue37373 :(
        if sys.version_info > (3, 8) and sys.platform.startswith("win"):
            asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

        self.init_syspath()
        super().initialize(argv)
        if hasattr(self, "load_config_environ"):
            self.load_config_environ()
        self.init_notebooks()
        self.init_writer()
        self.init_postprocessor()

    def init_syspath(self):
        """Add the cwd to the sys.path ($PYTHONPATH)"""
        sys.path.insert(0, os.getcwd())

    def init_notebooks(self):
        """Construct the list of notebooks.

        If notebooks are passed on the command-line,
        they override (rather than add) notebooks specified in config files.
        Glob each notebook to replace notebook patterns with filenames.
        """

        # Specifying notebooks on the command-line overrides (rather than
        # adds) the notebook list
        patterns = self.extra_args if self.extra_args else self.notebooks

        # Use glob to replace all the notebook patterns with filenames.
        filenames = []
        for pattern in patterns:
            # Use glob to find matching filenames.  Allow the user to convert
            # notebooks without having to type the extension.
            globbed_files = glob.glob(pattern, recursive=self.recursive_glob)
            globbed_files.extend(glob.glob(pattern + ".ipynb", recursive=self.recursive_glob))
            if not globbed_files:
                self.log.warning("pattern %r matched no files", pattern)

            for filename in globbed_files:
                if filename not in filenames:
                    filenames.append(filename)
        self.notebooks = filenames

    def init_writer(self):
        """Initialize the writer (which is stateless)"""
        self._writer_class_changed({"new": self.writer_class})
        if self.writer_factory:
            self.writer = self.writer_factory(parent=self)
            if hasattr(self.writer, "build_directory") and self.writer.build_directory != "":
                self.use_output_suffix = False

    def init_postprocessor(self):
        """Initialize the postprocessor (which is stateless)"""
        self._postprocessor_class_changed({"new": self.postprocessor_class})
        if self.postprocessor_factory:
            self.postprocessor = self.postprocessor_factory(parent=self)

    def start(self):
        """Run start after initialization process has completed"""
        super().start()
        self.convert_notebooks()

    def _notebook_filename_to_name(self, notebook_filename):
        """
        Returns the notebook name from the notebook filename by
        applying `output_base` pattern and stripping extension
        """
        basename = os.path.basename(notebook_filename)
        notebook_name = basename[: basename.rfind(".")]
        notebook_name = self.output_base.format(notebook_name=notebook_name)

        return notebook_name  # noqa: RET504

    def init_single_notebook_resources(self, notebook_filename):
        """Step 1: Initialize resources

        This initializes the resources dictionary for a single notebook.

        Returns
        -------
        dict
            resources dictionary for a single notebook that MUST include the following keys:
                - config_dir: the location of the Jupyter config directory
                - unique_key: the notebook name
                - output_files_dir: a directory where output files (not
                  including the notebook itself) should be saved
        """
        notebook_name = self._notebook_filename_to_name(notebook_filename)
        self.log.debug("Notebook name is '%s'", notebook_name)

        # first initialize the resources we want to use
        resources = {}
        resources["config_dir"] = self.config_dir
        resources["unique_key"] = notebook_name

        output_files_dir = self.output_files_dir.format(notebook_name=notebook_name)

        resources["output_files_dir"] = output_files_dir

        return resources

    def export_single_notebook(self, notebook_filename, resources, input_buffer=None):
        """Step 2: Export the notebook

        Exports the notebook to a particular format according to the specified
        exporter. This function returns the output and (possibly modified)
        resources from the exporter.

        Parameters
        ----------
        notebook_filename : str
            name of notebook file.
        resources : dict
        input_buffer :
            readable file-like object returning unicode.
            if not None, notebook_filename is ignored

        Returns
        -------
        output
        dict
            resources (possibly modified)
        """
        try:
            if input_buffer is not None:
                output, resources = self.exporter.from_file(input_buffer, resources=resources)
            else:
                output, resources = self.exporter.from_filename(
                    notebook_filename, resources=resources
                )
        except ConversionException:
            self.log.error("Error while converting '%s'", notebook_filename, exc_info=True)  # noqa: G201
            self.exit(1)

        return output, resources

    def write_single_notebook(self, output, resources):
        """Step 3: Write the notebook to file

        This writes output from the exporter to file using the specified writer.
        It returns the results from the writer.

        Parameters
        ----------
        output :
        resources : dict
            resources for a single notebook including name, config directory
            and directory to save output

        Returns
        -------
        file
            results from the specified writer output of exporter
        """

        if "unique_key" not in resources:
            msg = "unique_key MUST be specified in the resources, but it is not"
            raise KeyError(msg)

        notebook_name = resources["unique_key"]
        if self.use_output_suffix and self.output_base == "{notebook_name}":
            notebook_name += resources.get("output_suffix", "")

        if not self.writer:
            msg = "No writer object defined!"
            raise ValueError(msg)
        return self.writer.write(output, resources, notebook_name=notebook_name)

    def postprocess_single_notebook(self, write_results):
        """Step 4: Post-process the written file

        Only used if a postprocessor has been specified. After the
        converted notebook is written to a file in Step 3, this post-processes
        the notebook.
        """
        # Post-process if post processor has been defined.
        if hasattr(self, "postprocessor") and self.postprocessor:
            self.postprocessor(write_results)

    def convert_single_notebook(self, notebook_filename, input_buffer=None):
        """Convert a single notebook.

        Performs the following steps:

            1. Initialize notebook resources
            2. Export the notebook to a particular format
            3. Write the exported notebook to file
            4. (Maybe) postprocess the written file

        Parameters
        ----------
        notebook_filename : str
        input_buffer :
            If input_buffer is not None, conversion is done and the buffer is
            used as source into a file basenamed by the notebook_filename
            argument.
        """
        if input_buffer is None:
            self.log.info("Converting notebook %s to %s", notebook_filename, self.export_format)
        else:
            self.log.info("Converting notebook into %s", self.export_format)

        resources = self.init_single_notebook_resources(notebook_filename)
        output, resources = self.export_single_notebook(
            notebook_filename, resources, input_buffer=input_buffer
        )
        write_results = self.write_single_notebook(output, resources)
        self.postprocess_single_notebook(write_results)

    def convert_notebooks(self):
        """Convert the notebooks in the self.notebooks traitlet"""

        # no notebooks to convert!
        if len(self.notebooks) == 0 and not self.from_stdin:
            self.print_help()
            sys.exit(-1)

        if not self.export_format:
            msg = (
                "Please specify an output format with '--to <format>'."
                f"\nThe following formats are available: {get_export_names()}"
            )
            raise ValueError(msg)

        # initialize the exporter
        cls = get_exporter(self.export_format)
        self.exporter = cls(config=self.config)

        # strip duplicate extension from output_base, to avoid Basename.ext.ext
        if getattr(self.exporter, "file_extension", False):
            base, ext = os.path.splitext(self.output_base)
            if ext == self.exporter.file_extension:
                self.output_base = base

        # convert each notebook
        if not self.from_stdin:
            for notebook_filename in self.notebooks:
                self.convert_single_notebook(notebook_filename)
        else:
            input_buffer = unicode_stdin_stream()
            # default name when conversion from stdin
            self.convert_single_notebook("notebook.ipynb", input_buffer=input_buffer)
            input_buffer.close()

    def document_flag_help(self):
        """
        Return a string containing descriptions of all the flags.
        """
        flags = "The following flags are defined:\n\n"
        for flag, (cfg, fhelp) in self.flags.items():
            flags += f"{flag}\n"
            flags += indent(fill(fhelp, 80)) + "\n\n"
            flags += indent(fill("Long Form: " + str(cfg), 80)) + "\n\n"
        return flags

    def document_alias_help(self):
        """Return a string containing all of the aliases"""

        aliases = "The following aliases are defined:\n\n"
        for alias, longname in self.aliases.items():
            aliases += f"\t**{alias}** ({longname})\n\n"
        return aliases

    def document_config_options(self):
        """
        Provides a much improves version of the configuration documentation by
        breaking the configuration options into app, exporter, writer,
        preprocessor, postprocessor, and other sections.
        """
        categories = {
            category: [c for c in self._classes_inc_parents() if category in c.__name__.lower()]
            for category in ["app", "exporter", "writer", "preprocessor", "postprocessor"]
        }
        accounted_for = {c for category in categories.values() for c in category}
        categories["other"] = [c for c in self._classes_inc_parents() if c not in accounted_for]

        header = dedent(
            """
                        {section} Options
                        -----------------------

                        """
        )
        sections = ""
        for category in categories:
            sections += header.format(section=category.title())
            if category in ["exporter", "preprocessor", "writer"]:
                sections += f".. image:: _static/{category}_inheritance.png\n\n"
            sections += "\n".join(c.class_config_rst_doc() for c in categories[category])

        return sections.replace(" : ", r" \: ")


class DejavuApp(NbConvertApp):
    """A deja vu app."""

    def initialize(self, argv=None):
        """Initialize the app."""
        self.config.TemplateExporter.exclude_input = True
        self.config.TemplateExporter.exclude_output_prompt = True
        self.config.TemplateExporter.exclude_input_prompt = True
        self.config.ExecutePreprocessor.enabled = True
        self.config.WebPDFExporter.paginate = False
        self.config.QtPDFExporter.paginate = False

        super().initialize(argv)
        if hasattr(self, "load_config_environ"):
            self.load_config_environ()

    @default("export_format")
    def _default_export_format(self):
        return "html"


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

main = launch_new_instance = NbConvertApp.launch_instance
dejavu_main = DejavuApp.launch_instance
