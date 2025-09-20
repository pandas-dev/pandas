"""This module defines TemplateExporter, a highly configurable converter
that uses Jinja2 to export notebook files into different formats.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import html
import json
import os
import typing as t
import uuid
import warnings
from pathlib import Path

from jinja2 import (
    BaseLoader,
    ChoiceLoader,
    DictLoader,
    Environment,
    FileSystemLoader,
    TemplateNotFound,
)
from jupyter_core.paths import jupyter_path
from nbformat import NotebookNode
from traitlets import Bool, Dict, HasTraits, List, Unicode, default, observe, validate
from traitlets.config import Config
from traitlets.utils.importstring import import_item

from nbconvert import filters

from .exporter import Exporter

# Jinja2 extensions to load.
JINJA_EXTENSIONS = ["jinja2.ext.loopcontrols"]

ROOT = os.path.dirname(__file__)
DEV_MODE = os.path.exists(os.path.join(ROOT, "../../.git"))


default_filters = {
    "indent": filters.indent,
    "markdown2html": filters.markdown2html,
    "markdown2asciidoc": filters.markdown2asciidoc,
    "ansi2html": filters.ansi2html,
    "filter_data_type": filters.DataTypeFilter,
    "get_lines": filters.get_lines,
    "highlight2html": filters.Highlight2HTML,
    "highlight2latex": filters.Highlight2Latex,
    "ipython2python": filters.ipython2python,
    "posix_path": filters.posix_path,
    "markdown2latex": filters.markdown2latex,
    "markdown2rst": filters.markdown2rst,
    "comment_lines": filters.comment_lines,
    "strip_ansi": filters.strip_ansi,
    "strip_dollars": filters.strip_dollars,
    "strip_files_prefix": filters.strip_files_prefix,
    "html2text": filters.html2text,
    "add_anchor": filters.add_anchor,
    "ansi2latex": filters.ansi2latex,
    "wrap_text": filters.wrap_text,
    "escape_latex": filters.escape_latex,
    "citation2latex": filters.citation2latex,
    "path2url": filters.path2url,
    "add_prompts": filters.add_prompts,
    "ascii_only": filters.ascii_only,
    "prevent_list_blocks": filters.prevent_list_blocks,
    "get_metadata": filters.get_metadata,
    "convert_pandoc": filters.convert_pandoc,
    "json_dumps": json.dumps,
    # For removing any HTML
    "escape_html": lambda s: html.escape(str(s)),
    "escape_html_keep_quotes": lambda s: html.escape(str(s), quote=False),
    "escape_html_script": lambda s: s.replace("/", "\\/"),
    # For sanitizing HTML for any XSS
    "clean_html": filters.clean_html,
    "strip_trailing_newline": filters.strip_trailing_newline,
    "text_base64": filters.text_base64,
}


# copy of https://github.com/jupyter/jupyter_server/blob/b62458a7f5ad6b5246d2f142258dedaa409de5d9/jupyter_server/config_manager.py#L19
def recursive_update(target, new):
    """Recursively update one dictionary using another.
    None values will delete their keys.
    """
    for k, v in new.items():
        if isinstance(v, dict):
            if k not in target:
                target[k] = {}
            recursive_update(target[k], v)
            if not target[k]:
                # Prune empty subdicts
                del target[k]

        elif v is None:
            target.pop(k, None)

        else:
            target[k] = v
    return target  # return for convenience


# define function at the top level to avoid pickle errors
def deprecated(msg):
    """Emit a deprecation warning."""
    warnings.warn(msg, DeprecationWarning, stacklevel=2)


class ExtensionTolerantLoader(BaseLoader):
    """A template loader which optionally adds a given extension when searching.

    Constructor takes two arguments: *loader* is another Jinja loader instance
    to wrap. *extension* is the extension, which will be added to the template
    name if finding the template without it fails. This should include the dot,
    e.g. '.tpl'.
    """

    def __init__(self, loader, extension):
        """Initialize the loader."""
        self.loader = loader
        self.extension = extension

    def get_source(self, environment, template):
        """Get the source for a template."""
        try:
            return self.loader.get_source(environment, template)
        except TemplateNotFound:
            if template.endswith(self.extension):
                raise TemplateNotFound(template) from None
            return self.loader.get_source(environment, template + self.extension)

    def list_templates(self):
        """List available templates."""
        return self.loader.list_templates()


class TemplateExporter(Exporter):
    """
    Exports notebooks into other file formats.  Uses Jinja 2 templating engine
    to output new formats.  Inherit from this class if you are creating a new
    template type along with new filters/preprocessors.  If the filters/
    preprocessors provided by default suffice, there is no need to inherit from
    this class.  Instead, override the template_file and file_extension
    traits via a config file.

    Filters available by default for templates:

    {filters}
    """

    # finish the docstring
    __doc__ = (
        __doc__.format(filters="- " + "\n    - ".join(sorted(default_filters.keys())))
        if __doc__
        else None
    )

    _template_cached = None

    def _invalidate_template_cache(self, change=None):
        self._template_cached = None

    @property
    def template(self):
        if self._template_cached is None:
            self._template_cached = self._load_template()
        return self._template_cached

    _environment_cached = None

    def _invalidate_environment_cache(self, change=None):
        self._environment_cached = None
        self._invalidate_template_cache()

    @property
    def environment(self):
        if self._environment_cached is None:
            self._environment_cached = self._create_environment()
        return self._environment_cached

    @property
    def default_config(self):
        c = Config(
            {
                "RegexRemovePreprocessor": {"enabled": True},
                "TagRemovePreprocessor": {"enabled": True},
            }
        )
        if super().default_config:
            c2 = super().default_config.copy()
            c2.merge(c)
            c = c2
        return c

    template_name = Unicode(help="Name of the template to use").tag(
        config=True, affects_template=True
    )

    template_file = Unicode(None, allow_none=True, help="Name of the template file to use").tag(
        config=True, affects_template=True
    )

    raw_template = Unicode("", help="raw template string").tag(affects_environment=True)

    enable_async = Bool(False, help="Enable Jinja async template execution").tag(
        affects_environment=True
    )

    _last_template_file = ""
    _raw_template_key = "<memory>"

    @validate("template_name")
    def _template_name_validate(self, change):
        template_name = change["value"]
        if template_name and template_name.endswith(".tpl"):
            warnings.warn(
                f"5.x style template name passed '{self.template_name}'. Use --template-name for the template directory with a index.<ext>.j2 file and/or --template-file to denote a different template.",
                DeprecationWarning,
                stacklevel=2,
            )
            directory, self.template_file = os.path.split(self.template_name)
            if directory:
                directory, template_name = os.path.split(directory)
            if directory and os.path.isabs(directory):
                self.extra_template_basedirs = [directory]
        return template_name

    @observe("template_file")
    def _template_file_changed(self, change):
        new = change["new"]
        if new == "default":
            self.template_file = self.default_template  # type:ignore[attr-defined]
            return
        # check if template_file is a file path
        # rather than a name already on template_path
        full_path = os.path.abspath(new)
        if os.path.isfile(full_path):
            directory, self.template_file = os.path.split(full_path)
            self.extra_template_paths = [directory, *self.extra_template_paths]
            # While not strictly an invalid template file name, the extension hints that there isn't a template directory involved
            if self.template_file and self.template_file.endswith(".tpl"):
                warnings.warn(
                    f"5.x style template file passed '{new}'. Use --template-name for the template directory with a index.<ext>.j2 file and/or --template-file to denote a different template.",
                    DeprecationWarning,
                    stacklevel=2,
                )

    @default("template_file")
    def _template_file_default(self):
        if self.template_extension:
            return "index" + self.template_extension
        return None

    @observe("raw_template")
    def _raw_template_changed(self, change):
        if not change["new"]:
            self.template_file = self._last_template_file
        self._invalidate_template_cache()

    template_paths = List(["."]).tag(config=True, affects_environment=True)
    extra_template_basedirs = List(Unicode()).tag(config=True, affects_environment=True)
    extra_template_paths = List(Unicode()).tag(config=True, affects_environment=True)

    @default("extra_template_basedirs")
    def _default_extra_template_basedirs(self):
        return [os.getcwd()]

    # Extension that the template files use.
    template_extension = Unicode().tag(config=True, affects_environment=True)

    template_data_paths = List(
        jupyter_path("nbconvert", "templates"), help="Path where templates can be installed too."
    ).tag(affects_environment=True)

    @default("template_extension")
    def _template_extension_default(self):
        if self.file_extension:
            return self.file_extension + ".j2"
        return self.file_extension

    exclude_input = Bool(
        False, help="This allows you to exclude code cell inputs from all templates if set to True."
    ).tag(config=True)

    exclude_input_prompt = Bool(
        False, help="This allows you to exclude input prompts from all templates if set to True."
    ).tag(config=True)

    exclude_output = Bool(
        False,
        help="This allows you to exclude code cell outputs from all templates if set to True.",
    ).tag(config=True)

    exclude_output_prompt = Bool(
        False, help="This allows you to exclude output prompts from all templates if set to True."
    ).tag(config=True)

    exclude_output_stdin = Bool(
        True,
        help="This allows you to exclude output of stdin stream from lab template if set to True.",
    ).tag(config=True)

    exclude_code_cell = Bool(
        False, help="This allows you to exclude code cells from all templates if set to True."
    ).tag(config=True)

    exclude_markdown = Bool(
        False, help="This allows you to exclude markdown cells from all templates if set to True."
    ).tag(config=True)

    exclude_raw = Bool(
        False, help="This allows you to exclude raw cells from all templates if set to True."
    ).tag(config=True)

    exclude_unknown = Bool(
        False, help="This allows you to exclude unknown cells from all templates if set to True."
    ).tag(config=True)

    extra_loaders: List[t.Any] = List(
        help="Jinja loaders to find templates. Will be tried in order "
        "before the default FileSystem ones.",
    ).tag(affects_environment=True)

    filters = Dict(
        help="""Dictionary of filters, by name and namespace, to add to the Jinja
        environment."""
    ).tag(config=True, affects_environment=True)

    raw_mimetypes = List(
        Unicode(), help="""formats of raw cells to be included in this Exporter's output."""
    ).tag(config=True)

    @default("raw_mimetypes")
    def _raw_mimetypes_default(self):
        return [self.output_mimetype, ""]

    # TODO: passing config is wrong, but changing this revealed more complicated issues
    def __init__(self, config=None, **kw):
        """
        Public constructor

        Parameters
        ----------
        config : config
            User configuration instance.
        extra_loaders : list[of Jinja Loaders]
            ordered list of Jinja loader to find templates. Will be tried in order
            before the default FileSystem ones.
        template_file : str (optional, kw arg)
            Template to use when exporting.
        """
        super().__init__(config=config, **kw)

        self.observe(
            self._invalidate_environment_cache, list(self.traits(affects_environment=True))
        )
        self.observe(self._invalidate_template_cache, list(self.traits(affects_template=True)))

    def _load_template(self):
        """Load the Jinja template object from the template file

        This is triggered by various trait changes that would change the template.
        """

        # this gives precedence to a raw_template if present
        with self.hold_trait_notifications():
            if self.template_file and (self.template_file != self._raw_template_key):
                self._last_template_file = self.template_file
            if self.raw_template:
                self.template_file = self._raw_template_key

        if not self.template_file:
            msg = "No template_file specified!"
            raise ValueError(msg)

        # First try to load the
        # template by name with extension added, then try loading the template
        # as if the name is explicitly specified.
        template_file = self.template_file
        self.log.debug("Attempting to load template %s", template_file)
        self.log.debug("    template_paths: %s", os.pathsep.join(self.template_paths))
        return self.environment.get_template(template_file)

    def from_filename(  # type:ignore[override]
        self, filename: str, resources: dict[str, t.Any] | None = None, **kw: t.Any
    ) -> tuple[str, dict[str, t.Any]]:
        """Convert a notebook from a filename."""
        return super().from_filename(filename, resources, **kw)  # type:ignore[return-value]

    def from_file(  # type:ignore[override]
        self, file_stream: t.Any, resources: dict[str, t.Any] | None = None, **kw: t.Any
    ) -> tuple[str, dict[str, t.Any]]:
        """Convert a notebook from a file."""
        return super().from_file(file_stream, resources, **kw)  # type:ignore[return-value]

    def from_notebook_node(  # type:ignore[explicit-override, override]
        self, nb: NotebookNode, resources: dict[str, t.Any] | None = None, **kw: t.Any
    ) -> tuple[str, dict[str, t.Any]]:
        """
        Convert a notebook from a notebook node instance.

        Parameters
        ----------
        nb : :class:`~nbformat.NotebookNode`
            Notebook node
        resources : dict
            Additional resources that can be accessed read/write by
            preprocessors and filters.
        """
        nb_copy, resources = super().from_notebook_node(nb, resources, **kw)
        resources.setdefault("raw_mimetypes", self.raw_mimetypes)
        resources.setdefault("output_mimetype", self.output_mimetype)
        resources["global_content_filter"] = {
            "include_code": not self.exclude_code_cell,
            "include_markdown": not self.exclude_markdown,
            "include_raw": not self.exclude_raw,
            "include_unknown": not self.exclude_unknown,
            "include_input": not self.exclude_input,
            "include_output": not self.exclude_output,
            "include_output_stdin": not self.exclude_output_stdin,
            "include_input_prompt": not self.exclude_input_prompt,
            "include_output_prompt": not self.exclude_output_prompt,
            "no_prompt": self.exclude_input_prompt and self.exclude_output_prompt,
        }

        # Top level variables are passed to the template_exporter here.
        output = self.template.render(nb=nb_copy, resources=resources)
        output = output.lstrip("\r\n")
        return output, resources

    def _register_filter(self, environ, name, jinja_filter):
        """
        Register a filter.
        A filter is a function that accepts and acts on one string.
        The filters are accessible within the Jinja templating engine.

        Parameters
        ----------
        name : str
            name to give the filter in the Jinja engine
        filter : filter
        """
        if jinja_filter is None:
            msg = "filter"
            raise TypeError(msg)
        isclass = isinstance(jinja_filter, type)
        constructed = not isclass

        # Handle filter's registration based on it's type
        if constructed and isinstance(jinja_filter, (str,)):
            # filter is a string, import the namespace and recursively call
            # this register_filter method
            filter_cls = import_item(jinja_filter)
            return self._register_filter(environ, name, filter_cls)

        if constructed and callable(jinja_filter):
            # filter is a function, no need to construct it.
            environ.filters[name] = jinja_filter
            return jinja_filter

        if isclass and issubclass(jinja_filter, HasTraits):
            # filter is configurable.  Make sure to pass in new default for
            # the enabled flag if one was specified.
            filter_instance = jinja_filter(parent=self)
            self._register_filter(environ, name, filter_instance)
            return None

        if isclass:
            # filter is not configurable, construct it
            filter_instance = jinja_filter()
            self._register_filter(environ, name, filter_instance)
            return None

        # filter is an instance of something without a __call__
        # attribute.
        msg = "filter"
        raise TypeError(msg)

    def register_filter(self, name, jinja_filter):
        """
        Register a filter.
        A filter is a function that accepts and acts on one string.
        The filters are accessible within the Jinja templating engine.

        Parameters
        ----------
        name : str
            name to give the filter in the Jinja engine
        filter : filter
        """
        return self._register_filter(self.environment, name, jinja_filter)

    def default_filters(self):
        """Override in subclasses to provide extra filters.

        This should return an iterable of 2-tuples: (name, class-or-function).
        You should call the method on the parent class and include the filters
        it provides.

        If a name is repeated, the last filter provided wins. Filters from
        user-supplied config win over filters provided by classes.
        """
        return default_filters.items()

    def _create_environment(self):
        """
        Create the Jinja templating environment.
        """
        paths = self.template_paths
        self.log.debug("Template paths:\n\t%s", "\n\t".join(paths))

        loaders = [
            *self.extra_loaders,
            ExtensionTolerantLoader(FileSystemLoader(paths), self.template_extension),
            DictLoader({self._raw_template_key: self.raw_template}),
        ]
        environment = Environment(  # noqa: S701
            loader=ChoiceLoader(loaders),
            extensions=JINJA_EXTENSIONS,
            enable_async=self.enable_async,
        )

        environment.globals["uuid4"] = uuid.uuid4

        # Add default filters to the Jinja2 environment
        for key, value in self.default_filters():
            self._register_filter(environment, key, value)

        # Load user filters.  Overwrite existing filters if need be.
        if self.filters:
            for key, user_filter in self.filters.items():
                self._register_filter(environment, key, user_filter)

        return environment

    def _init_preprocessors(self):
        super()._init_preprocessors()
        conf = self._get_conf()
        preprocessors = conf.get("preprocessors", {})
        # preprocessors is a dict for three reasons
        #  * We rely on recursive_update, which can only merge dicts, lists will be overwritten
        #  * We can use the key with numerical prefixing to guarantee ordering (/etc/*.d/XY-file style)
        #  * We can disable preprocessors by overwriting the value with None
        for _, preprocessor in sorted(preprocessors.items(), key=lambda x: x[0]):
            if preprocessor is not None:
                kwargs = preprocessor.copy()
                preprocessor_cls = kwargs.pop("type")
                preprocessor_cls = import_item(preprocessor_cls)
                if preprocessor_cls.__name__ in self.config:
                    kwargs.update(self.config[preprocessor_cls.__name__])
                preprocessor = preprocessor_cls(**kwargs)  # noqa: PLW2901
                self.register_preprocessor(preprocessor)

    def _get_conf(self):
        conf: dict[str, t.Any] = {}  # the configuration once all conf files are merged
        for path in map(Path, self.template_paths):
            conf_path = path / "conf.json"
            if conf_path.exists():
                with conf_path.open() as f:
                    conf = recursive_update(conf, json.load(f))
        return conf

    @default("template_paths")
    def _template_paths(self, prune=True, root_dirs=None):
        paths = []
        root_dirs = self.get_prefix_root_dirs()
        template_names = self.get_template_names()
        for template_name in template_names:
            for base_dir in self.extra_template_basedirs:
                path = os.path.join(base_dir, template_name)
                try:
                    if not prune or os.path.exists(path):
                        paths.append(path)
                except PermissionError:
                    pass
            for root_dir in root_dirs:
                base_dir = os.path.join(root_dir, "nbconvert", "templates")
                path = os.path.join(base_dir, template_name)
                try:
                    if not prune or os.path.exists(path):
                        paths.append(path)
                except PermissionError:
                    pass

        for root_dir in root_dirs:
            # we include root_dir for when we want to be very explicit, e.g.
            # {% extends 'nbconvert/templates/classic/base.html' %}
            paths.append(root_dir)
            # we include base_dir for when we want to be explicit, but less than root_dir, e.g.
            # {% extends 'classic/base.html' %}
            base_dir = os.path.join(root_dir, "nbconvert", "templates")
            paths.append(base_dir)

            compatibility_dir = os.path.join(root_dir, "nbconvert", "templates", "compatibility")
            paths.append(compatibility_dir)

        additional_paths = []
        for path in self.template_data_paths:
            if not prune or os.path.exists(path):
                additional_paths.append(path)

        return paths + self.extra_template_paths + additional_paths

    @classmethod
    def get_compatibility_base_template_conf(cls, name):
        """Get the base template config."""
        # Hard-coded base template confs to use for backwards compatibility for 5.x-only templates
        if name == "display_priority":
            return {"base_template": "base"}
        if name == "full":
            return {"base_template": "classic", "mimetypes": {"text/html": True}}
        return None

    def get_template_names(self):
        """Finds a list of template names where each successive template name is the base template"""
        template_names = []
        root_dirs = self.get_prefix_root_dirs()
        base_template: str | None = self.template_name
        merged_conf: dict[str, t.Any] = {}  # the configuration once all conf files are merged
        while base_template is not None:
            template_names.append(base_template)
            conf: dict[str, t.Any] = {}
            found_at_least_one = False
            for base_dir in self.extra_template_basedirs:
                template_dir = os.path.join(base_dir, base_template)
                if os.path.exists(template_dir):
                    found_at_least_one = True
                conf_file = os.path.join(template_dir, "conf.json")
                if os.path.exists(conf_file):
                    with open(conf_file) as f:
                        conf = recursive_update(json.load(f), conf)
            for root_dir in root_dirs:
                template_dir = os.path.join(root_dir, "nbconvert", "templates", base_template)
                if os.path.exists(template_dir):
                    found_at_least_one = True
                conf_file = os.path.join(template_dir, "conf.json")
                if os.path.exists(conf_file):
                    with open(conf_file) as f:
                        conf = recursive_update(json.load(f), conf)
            if not found_at_least_one:
                # Check for backwards compatibility template names
                for root_dir in root_dirs:
                    compatibility_file = base_template + ".tpl"
                    compatibility_path = os.path.join(
                        root_dir, "nbconvert", "templates", "compatibility", compatibility_file
                    )
                    if os.path.exists(compatibility_path):
                        found_at_least_one = True
                        warnings.warn(
                            f"5.x template name passed '{self.template_name}'. Use 'lab' or 'classic' for new template usage.",
                            DeprecationWarning,
                            stacklevel=2,
                        )
                        self.template_file = compatibility_file
                        conf = self.get_compatibility_base_template_conf(base_template)
                        self.template_name = t.cast(str, conf.get("base_template"))
                        break
                if not found_at_least_one:
                    paths = "\n\t".join(root_dirs)
                    msg = f"No template sub-directory with name {base_template!r} found in the following paths:\n\t{paths}"
                    raise ValueError(msg)
            merged_conf = recursive_update(dict(conf), merged_conf)
            base_template = t.cast(t.Any, conf.get("base_template"))
        conf = merged_conf
        mimetypes = [mimetype for mimetype, enabled in conf.get("mimetypes", {}).items() if enabled]
        if self.output_mimetype and self.output_mimetype not in mimetypes and mimetypes:
            supported_mimetypes = "\n\t".join(mimetypes)
            msg = f"Unsupported mimetype {self.output_mimetype!r} for template {self.template_name!r}, mimetypes supported are: \n\t{supported_mimetypes}"
            raise ValueError(msg)
        return template_names

    def get_prefix_root_dirs(self):
        """Get the prefix root dirs."""
        # We look at the usual jupyter locations, and for development purposes also
        # relative to the package directory (first entry, meaning with highest precedence)
        root_dirs = []
        if DEV_MODE:
            root_dirs.append(os.path.abspath(os.path.join(ROOT, "..", "..", "share", "jupyter")))
        root_dirs.extend(jupyter_path())
        return root_dirs

    def _init_resources(self, resources):
        resources = super()._init_resources(resources)
        resources["deprecated"] = deprecated
        return resources
