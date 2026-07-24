"""Input/Output files"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from docutils.io import FileInput
from docutils.readers import standalone
from docutils.transforms.references import DanglingReferences
from docutils.writers import UnfilteredWriter

from sphinx.deprecation import RemovedInSphinx10Warning
from sphinx.transforms import SphinxTransformer
from sphinx.util import logging
from sphinx.util.docutils import LoggingReporter

if TYPE_CHECKING:
    from typing import Any

    from docutils import nodes
    from docutils.io import Input
    from docutils.parsers import Parser
    from docutils.transforms import Transform

    from sphinx.environment import BuildEnvironment
    from sphinx.util.docutils import _DocutilsSettings


logger = logging.getLogger(__name__)

warnings.warn('sphinx.io is deprecated', RemovedInSphinx10Warning, stacklevel=2)


class SphinxBaseReader(standalone.Reader):  # type: ignore[type-arg]
    """A base class of readers for Sphinx.

    This replaces reporter by Sphinx's on generating document.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        warnings.warn(
            'sphinx.io.SphinxBaseReader is deprecated',
            RemovedInSphinx10Warning,
            stacklevel=2,
        )

    transforms: list[type[Transform]] = []

    def get_transforms(self) -> list[type[Transform]]:
        transforms = super().get_transforms() + self.transforms

        # remove transforms which is not needed for Sphinx
        unused = [DanglingReferences]
        for transform in unused:
            if transform in transforms:
                transforms.remove(transform)

        return transforms

    def new_document(self) -> nodes.document:
        """Creates a new document object which has a special reporter object good
        for logging.
        """
        document = super().new_document()

        # substitute transformer
        document.transformer = SphinxTransformer(document)
        document.transformer.set_environment(self.settings.env)

        # substitute reporter
        reporter = document.reporter
        document.reporter = LoggingReporter.from_reporter(reporter)

        return document


class SphinxStandaloneReader(SphinxBaseReader):
    """A basic document reader for Sphinx."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        warnings.warn(
            'sphinx.io.SphinxStandaloneReader is deprecated',
            RemovedInSphinx10Warning,
            stacklevel=2,
        )

    def _setup_transforms(self, transforms: list[type[Transform]], /) -> None:
        self.transforms = self.transforms + transforms

    def read(
        self,
        source: Input,  # type: ignore[type-arg]
        parser: Parser,
        settings: _DocutilsSettings,
    ) -> nodes.document:
        self.source = source
        if not self.parser:
            self.parser = parser
        self.settings = settings
        self.input = self.read_source(settings.env)
        self.parse()
        assert self.document is not None
        return self.document

    def read_source(self, env: BuildEnvironment) -> str:
        """Read content from source and do post-process."""
        assert self.source is not None
        content = self.source.read()

        # emit "source-read" event
        arg = [content]
        env.events.emit('source-read', env.current_document.docname, arg)
        return arg[0]


class SphinxDummyWriter(UnfilteredWriter):  # type: ignore[type-arg]
    """Dummy writer module used for generating doctree."""

    def __init__(self) -> None:
        super().__init__()
        warnings.warn(
            'sphinx.io.SphinxDummyWriter is deprecated',
            RemovedInSphinx10Warning,
            stacklevel=2,
        )

    supported = ('html',)  # needed to keep "meta" nodes

    def translate(self) -> None:
        pass


def SphinxDummySourceClass(source: Any, *args: Any, **kwargs: Any) -> Any:
    """Bypass source object as is to cheat Publisher."""
    warnings.warn(
        'sphinx.io.SphinxDummySourceClass is deprecated',
        RemovedInSphinx10Warning,
        stacklevel=2,
    )
    return source


class SphinxFileInput(FileInput):
    """A basic FileInput for Sphinx."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs['error_handler'] = 'sphinx'
        super().__init__(*args, **kwargs)
        warnings.warn(
            'sphinx.io.SphinxFileInput is deprecated',
            RemovedInSphinx10Warning,
            stacklevel=2,
        )
