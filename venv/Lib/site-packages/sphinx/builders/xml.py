"""Docutils-native XML and pseudo-XML builders."""

from __future__ import annotations

from typing import TYPE_CHECKING

from docutils import nodes
from docutils.writers.docutils_xml import XMLTranslator

from sphinx.builders import Builder
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.osutil import _last_modified_time

if TYPE_CHECKING:
    from collections.abc import Iterator

    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata

logger = logging.getLogger(__name__)


class XMLBuilder(Builder):
    """Builds Docutils-native XML."""

    name = 'xml'
    format = 'xml'
    epilog = __('The XML files are in %(outdir)s.')

    out_suffix = '.xml'
    allow_parallel = True

    default_translator_class = XMLTranslator

    def init(self) -> None:
        pass

    def get_outdated_docs(self) -> Iterator[str]:
        for docname in self.env.found_docs:
            if docname not in self.env.all_docs:
                yield docname
                continue
            targetname = self.outdir / (docname + self.out_suffix)
            try:
                targetmtime = _last_modified_time(targetname)
            except Exception:
                targetmtime = 0
            try:
                srcmtime = _last_modified_time(self.env.doc2path(docname))
                if srcmtime > targetmtime:
                    yield docname
            except OSError:
                # source doesn't exist anymore
                pass

    def get_target_uri(self, docname: str, typ: str | None = None) -> str:
        return docname

    def write_doc(self, docname: str, doctree: nodes.document) -> None:
        # work around multiple string % tuple issues in docutils;
        # replace tuples in attribute values with lists
        doctree = doctree.deepcopy()
        for domain in self.env.domains.sorted():
            doctree[f'xmlns:{domain.name}'] = 'https://www.sphinx-doc.org/'
        for node in doctree.findall(nodes.Element):
            for att, value in node.attributes.items():
                if isinstance(value, tuple):
                    node.attributes[att] = list(value)
                value = node.attributes[att]
                if isinstance(value, list):
                    for i, val in enumerate(value):
                        if isinstance(val, tuple):
                            value[i] = list(val)
        output = self._translate(doctree)
        out_file_name = self.outdir / (docname + self.out_suffix)
        out_file_name.parent.mkdir(parents=True, exist_ok=True)
        try:
            out_file_name.write_text(output, encoding='utf-8')
        except OSError as err:
            logger.warning(__('error writing file %s: %s'), out_file_name, err)

    def _translate(self, doctree: nodes.document) -> str:
        doctree.settings.newlines = doctree.settings.indents = self.config.xml_pretty
        doctree.settings.xml_declaration = True
        doctree.settings.doctype_declaration = True

        # copied from docutils.writers.docutils_xml.Writer.translate()
        # so that we can override the translator class
        visitor: XMLTranslator = self.create_translator(doctree)  # type: ignore[assignment]
        doctree.walkabout(visitor)
        return ''.join(visitor.output)  # ty: ignore[unresolved-attribute]

    def finish(self) -> None:
        pass


class PseudoXMLBuilder(XMLBuilder):
    """Builds pseudo-XML for display purposes."""

    name = 'pseudoxml'
    format = 'pseudoxml'
    epilog = __('The pseudo-XML files are in %(outdir)s.')

    out_suffix = '.pseudoxml'

    def _translate(self, doctree: nodes.document) -> str:
        return doctree.pformat()


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_builder(XMLBuilder)
    app.add_builder(PseudoXMLBuilder)

    app.add_config_value('xml_pretty', True, 'env', types=frozenset({bool}))

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
