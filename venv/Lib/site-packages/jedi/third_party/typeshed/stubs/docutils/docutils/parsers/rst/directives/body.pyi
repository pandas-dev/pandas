from collections.abc import Callable
from typing import ClassVar, Final
from typing_extensions import TypeAlias

from docutils import nodes
from docutils.parsers.rst import Directive

__docformat__: Final = "reStructuredText"

_DirectiveFn: TypeAlias = Callable[[str], str | list[str]]

class BasePseudoSection(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    node_class: ClassVar[type[nodes.Node] | None]
    invalid_parents: ClassVar[
        tuple[
            type[nodes.SubStructural],
            type[nodes.Bibliographic],
            type[nodes.Decorative],
            type[nodes.Body],
            type[nodes.Part],
            type[nodes.topic],
        ]
    ]
    def run(self): ...

class Topic(BasePseudoSection):
    node_class: ClassVar[type[nodes.Node]]

class Sidebar(BasePseudoSection):
    node_class: ClassVar[type[nodes.Node]]
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class LineBlock(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class ParsedLiteral(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class CodeBlock(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class MathBlock(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class Rubric(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class BlockQuote(Directive):
    classes: ClassVar[list[str]]
    def run(self): ...

class Epigraph(BlockQuote):
    classes: ClassVar[list[str]]

class Highlights(BlockQuote):
    classes: ClassVar[list[str]]

class PullQuote(BlockQuote):
    classes: ClassVar[list[str]]

class Compound(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...

class Container(Directive):
    option_spec: ClassVar[dict[str, _DirectiveFn]]
    def run(self): ...
