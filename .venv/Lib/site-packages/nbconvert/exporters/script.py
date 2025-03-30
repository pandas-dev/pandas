"""Generic script exporter class for any kernel language"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import sys

if sys.version_info < (3, 10):
    from importlib_metadata import entry_points  # type:ignore[import-not-found]
else:
    from importlib.metadata import entry_points
from traitlets import Dict, default

from .base import get_exporter
from .templateexporter import TemplateExporter


class ScriptExporter(TemplateExporter):
    """A script exporter."""

    # Caches of already looked-up and instantiated exporters for delegation:
    _exporters = Dict()
    _lang_exporters = Dict()
    export_from_notebook = "Script"

    @default("template_file")
    def _template_file_default(self):
        return "script.j2"

    @default("template_name")
    def _template_name_default(self):
        return "script"

    def _get_language_exporter(self, lang_name):
        """Find an exporter for the language name from notebook metadata.

        Uses the nbconvert.exporters.script group of entry points.
        Returns None if no exporter is found.
        """
        if lang_name not in self._lang_exporters:
            try:
                exporters = entry_points(group="nbconvert.exporters.script")
                exporter = [e for e in exporters if e.name == lang_name][0].load()  # noqa: RUF015
            except (KeyError, IndexError):
                self._lang_exporters[lang_name] = None
            else:
                # TODO: passing config is wrong, but changing this revealed more complicated issues
                self._lang_exporters[lang_name] = exporter(config=self.config, parent=self)
        return self._lang_exporters[lang_name]

    def from_notebook_node(self, nb, resources=None, **kw):
        """Convert from notebook node."""
        langinfo = nb.metadata.get("language_info", {})

        # delegate to custom exporter, if specified
        exporter_name = langinfo.get("nbconvert_exporter")
        if exporter_name and exporter_name != "script":
            self.log.debug("Loading script exporter: %s", exporter_name)
            if exporter_name not in self._exporters:
                exporter = get_exporter(exporter_name)
                # TODO: passing config is wrong, but changing this revealed more complicated issues
                self._exporters[exporter_name] = exporter(config=self.config, parent=self)
            exporter = self._exporters[exporter_name]
            return exporter.from_notebook_node(nb, resources, **kw)

        # Look up a script exporter for this notebook's language
        lang_name = langinfo.get("name")
        if lang_name:
            self.log.debug("Using script exporter for language: %s", lang_name)
            exporter = self._get_language_exporter(lang_name)
            if exporter is not None:
                return exporter.from_notebook_node(nb, resources, **kw)

        # Fall back to plain script export
        self.file_extension = langinfo.get("file_extension", ".txt")
        self.output_mimetype = langinfo.get("mimetype", "text/plain")
        return super().from_notebook_node(nb, resources, **kw)
