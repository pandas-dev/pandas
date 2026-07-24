from _typeshed import Incomplete

from antlr4.tree.ParseTreePatternMatcher import ParseTreePatternMatcher as ParseTreePatternMatcher
from antlr4.tree.Tree import ParseTree as ParseTree
from antlr4.xpath.XPathLexer import XPathLexer as XPathLexer

class ParseTreePattern:
    __slots__ = ("matcher", "patternRuleIndex", "pattern", "patternTree")
    matcher: Incomplete
    patternRuleIndex: Incomplete
    pattern: Incomplete
    patternTree: Incomplete
    def __init__(self, matcher: ParseTreePatternMatcher, pattern: str, patternRuleIndex: int, patternTree: ParseTree) -> None: ...
    def match(self, tree: ParseTree): ...
    def matches(self, tree: ParseTree): ...
    def findAll(self, tree: ParseTree, xpath: str): ...
