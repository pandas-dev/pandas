import os
from traitlets import (
    HasTraits, Dict, Unicode, List, Bool,
    observe, default
)
from jupyter_core.paths import jupyter_path
from jupyter_server.transutils import _i18n
from jupyter_server.utils import url_path_join


class NotebookAppTraits(HasTraits):

    ignore_minified_js = Bool(False,
                              config=True,
                              help=_i18n(
                                  'Deprecated: Use minified JS file or not, mainly use during dev to avoid JS recompilation'),
                              )

    jinja_environment_options = Dict(config=True,
                                     help=_i18n("Supply extra arguments that will be passed to Jinja environment."))

    jinja_template_vars = Dict(
        config=True,
        help=_i18n(
            "Extra variables to supply to jinja templates when rendering."),
    )

    enable_mathjax = Bool(True, config=True,
                          help="""Whether to enable MathJax for typesetting math/TeX

        MathJax is the javascript library Jupyter uses to render math/LaTeX. It is
        very large, so you may want to disable it if you have a slow internet
        connection, or for offline use of the notebook.

        When disabled, equations etc. will appear as their untransformed TeX source.
        """
                          )

    @observe('enable_mathjax')
    def _update_enable_mathjax(self, change):
        """set mathjax url to empty if mathjax is disabled"""
        if not change['new']:
            self.mathjax_url = u''

    extra_static_paths = List(Unicode(), config=True,
                              help="""Extra paths to search for serving static files.

        This allows adding javascript/css to be available from the notebook server machine,
        or overriding individual files in the IPython"""
                              )

    @property
    def static_file_path(self):
        """return extra paths + the default location"""
        return self.extra_static_paths

    static_custom_path = List(Unicode(),
                              help=_i18n(
                                  """Path to search for custom.js, css""")
                              )

    @default('static_custom_path')
    def _default_static_custom_path(self):
        return [
            os.path.join(self.config_dir, 'custom')
        ]

    extra_template_paths = List(Unicode(), config=True,
                                help=_i18n("""Extra paths to search for serving jinja templates.

        Can be used to override templates from notebook.templates.""")
                                )

    @property
    def template_file_path(self):
        """return extra paths + the default locations"""
        return self.extra_template_paths

    extra_nbextensions_path = List(Unicode(), config=True,
                                   help=_i18n(
                                       """extra paths to look for Javascript notebook extensions""")
                                   )

    @property
    def nbextensions_path(self):
        """The path to look for Javascript notebook extensions"""
        path = self.extra_nbextensions_path + jupyter_path('nbextensions')
        # FIXME: remove IPython nbextensions path after a migration period
        try:
            from IPython.paths import get_ipython_dir
        except ImportError:
            pass
        else:
            path.append(os.path.join(get_ipython_dir(), 'nbextensions'))
        return path

    mathjax_url = Unicode("", config=True,
                          help="""A custom url for MathJax.js.
        Should be in the form of a case-sensitive url to MathJax,
        for example:  /static/components/MathJax/MathJax.js
        """
                          )

    @property
    def static_url_prefix(self):
        """Get the static url prefix for serving static files."""
        return super(NotebookAppTraits, self).static_url_prefix

    @default('mathjax_url')
    def _default_mathjax_url(self):
        if not self.enable_mathjax:
            return u''
        static_url_prefix = self.static_url_prefix
        return url_path_join(static_url_prefix, 'components', 'MathJax', 'MathJax.js')

    @observe('mathjax_url')
    def _update_mathjax_url(self, change):
        new = change['new']
        if new and not self.enable_mathjax:
            # enable_mathjax=False overrides mathjax_url
            self.mathjax_url = u''
        else:
            self.log.info(_i18n("Using MathJax: %s"), new)

    mathjax_config = Unicode("TeX-AMS-MML_HTMLorMML-full,Safe", config=True,
                             help=_i18n(
                                 """The MathJax.js configuration file that is to be used.""")
                             )

    @observe('mathjax_config')
    def _update_mathjax_config(self, change):
        self.log.info(
            _i18n("Using MathJax configuration file: %s"), change['new'])

    quit_button = Bool(True, config=True,
                       help="""If True, display a button in the dashboard to quit
        (shutdown the notebook server)."""
                       )

    nbserver_extensions = Dict({}, config=True,
                               help=(_i18n("Dict of Python modules to load as notebook server extensions."
                                           "Entry values can be used to enable and disable the loading of"
                                           "the extensions. The extensions will be loaded in alphabetical "
                                           "order."))
                               )
