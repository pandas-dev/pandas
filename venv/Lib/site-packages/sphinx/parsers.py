"""A Base class for additional parsers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import docutils.parsers
import docutils.parsers.rst
from docutils.parsers.rst import states
from docutils.statemachine import StringList
from docutils.transforms.universal import SmartQuotes

from sphinx.deprecation import _deprecation_warning
from sphinx.util.rst import _append_epilogue, _prepend_prologue

if TYPE_CHECKING:
    from docutils import nodes
    from docutils.transforms import Transform

    from sphinx.application import Sphinx
    from sphinx.config import Config
    from sphinx.environment import BuildEnvironment
    from sphinx.util.typing import ExtensionMetadata


class Parser(docutils.parsers.Parser):
    """A base class for source parsers.

    Additional parsers should inherit from this class instead of
    ``docutils.parsers.Parser``.
    This class provides access to core Sphinx objects; *config* and *env*.
    """

    _config: Config
    _env: BuildEnvironment

    @property
    def config(self) -> Config:
        """The config object."""
        cls_module = self.__class__.__module__
        cls_name = self.__class__.__qualname__
        _deprecation_warning(cls_module, f'{cls_name}.config', remove=(10, 0))
        return self._config

    @property
    def env(self) -> BuildEnvironment:
        """The environment object."""
        cls_module = self.__class__.__module__
        cls_name = self.__class__.__qualname__
        _deprecation_warning(cls_module, f'{cls_name}.env', remove=(10, 0))
        return self._env

    def set_application(self, app: Sphinx) -> None:
        """set_application will be called from Sphinx to set app and other instance variables

        :param sphinx.application.Sphinx app: Sphinx application object
        """
        cls_module = self.__class__.__module__
        cls_name = self.__class__.__qualname__
        _deprecation_warning(cls_module, f'{cls_name}.set_application', remove=(10, 0))
        self._config = app.config
        self._env = app.env


class RSTParser(docutils.parsers.rst.Parser, Parser):
    """A reST parser for Sphinx."""

    def get_transforms(self) -> list[type[Transform]]:
        """Sphinx's reST parser replaces a transform class for smart-quotes by its own

        refs: sphinx.io.SphinxStandaloneReader
        """
        transforms = super(RSTParser, RSTParser()).get_transforms()
        transforms.remove(SmartQuotes)
        return transforms

    def parse(self, inputstring: str | StringList, document: nodes.document) -> None:
        """Parse text and generate a document tree."""
        self.setup_parse(inputstring, document)  # type: ignore[arg-type]
        self.statemachine = states.RSTStateMachine(
            state_classes=self.state_classes,
            initial_state=self.initial_state,
            debug=document.reporter.debug_flag,
        )

        # preprocess inputstring
        if isinstance(inputstring, str):
            lines = docutils.statemachine.string2lines(
                inputstring,
                tab_width=document.settings.tab_width,
                convert_whitespace=True,
            )

            inputlines = StringList(lines, document.current_source)
        else:
            inputlines = inputstring

        self.decorate(inputlines)
        self.statemachine.run(inputlines, document, inliner=self.inliner)
        self.finish_parse()

    def decorate(self, content: StringList) -> None:
        """Preprocess reStructuredText content before parsing."""
        _prepend_prologue(content, self._config.rst_prolog)
        _append_epilogue(content, self._config.rst_epilog)


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_source_parser(RSTParser)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
