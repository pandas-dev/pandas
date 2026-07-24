from typing import ClassVar

from ..lexer import RegexLexer

__all__ = ["VisualPrologLexer", "VisualPrologGrammarLexer"]

class VisualPrologBaseLexer(RegexLexer):
    minorendkw: ClassVar[tuple[str, ...]]
    minorkwexp: ClassVar[tuple[str, ...]]
    dockw: ClassVar[tuple[str, ...]]

class VisualPrologLexer(VisualPrologBaseLexer):
    majorkw: ClassVar[tuple[str, ...]]
    minorkw: ClassVar[tuple[str, ...]]
    directivekw: ClassVar[tuple[str, ...]]

class VisualPrologGrammarLexer(VisualPrologBaseLexer):
    majorkw: ClassVar[tuple[str, ...]]
    directivekw: ClassVar[tuple[str, ...]]
