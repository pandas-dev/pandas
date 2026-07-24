from _typeshed import StrPath
from pathlib import Path
from typing import Any, Final, Generic, TypedDict, TypeVar, type_check_only
from typing_extensions import Required

from docutils import Component, nodes
from docutils.frontend import Values
from docutils.io import Output
from docutils.languages import LanguageImporter

_S = TypeVar("_S")

__docformat__: Final = "reStructuredText"

# It would probably be better to specialize writers for subclasses,
# but this gives us all possible Writer items w/o instance checks
@type_check_only
class _WriterParts(TypedDict, total=False):
    # Parts Provided by All Writers https://docutils.sourceforge.io/docs/api/publisher.html#parts-provided-by-all-writers

    # See Writer.assemble_parts
    whole: Required[str | bytes]
    encoding: Required[str]
    errors: Required[str]
    version: Required[str]

    # Parts Provided by the HTML Writers https://docutils.sourceforge.io/docs/api/publisher.html#parts-provided-by-the-html-writers

    # HTML4 Writer https://docutils.sourceforge.io/docs/api/publisher.html#html4-writer
    # + HTML5 Writer https://docutils.sourceforge.io/docs/api/publisher.html#html5-writer
    body: str
    body_prefix: str
    body_pre_docinfo: str
    body_suffix: str
    docinfo: str
    footer: str
    fragment: str
    head: str
    head_prefix: str
    header: str
    html_body: str
    html_head: str
    html_prolog: str
    html_subtitle: str
    html_title: str
    meta: str
    stylesheet: str
    subtitle: str
    title: str
    # PEP/HTML Writer https://docutils.sourceforge.io/docs/api/publisher.html#pep-html-writer
    # + S5/HTML Writer https://docutils.sourceforge.io/docs/api/publisher.html#s5-html-writer
    pepnum: str

    # Parts Provided by the (Xe)LaTeX Writers https://docutils.sourceforge.io/docs/api/publisher.html#parts-provided-by-the-xe-latex-writers

    # (commenting out those already included)
    abstract: str
    # body: str
    # body_pre_docinfo: str
    dedication: str
    # docinfo: str
    fallbacks: str
    # head_prefix: str
    latex_preamble: str
    pdfsetup: str
    requirements: str
    # stylesheet: str
    # subtitle: str
    # title: str
    titledata: str

class Writer(Component, Generic[_S]):
    parts: _WriterParts
    language: LanguageImporter | None = None
    document: nodes.document | None = None
    destination: Output | None = None
    output: _S | None = None
    def __init__(self) -> None: ...
    def write(self, document: nodes.document, destination: Output) -> str | bytes | None: ...
    def translate(self) -> None: ...
    def assemble_parts(self) -> None: ...

class UnfilteredWriter(Writer[_S]): ...

class DoctreeTranslator(nodes.NodeVisitor):
    settings: Values
    def __init__(self, document: nodes.document) -> None: ...
    def uri2path(self, uri: str, output_path: StrPath | None = None) -> Path: ...

WRITER_ALIASES: Final[dict[str, str]]

def get_writer_class(writer_name: str) -> type[Writer[Any]]: ...
