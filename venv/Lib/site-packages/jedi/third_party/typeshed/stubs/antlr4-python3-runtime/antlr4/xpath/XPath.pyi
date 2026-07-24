from _typeshed import Incomplete

from antlr4 import (
    DFA as DFA,
    CommonTokenStream as CommonTokenStream,
    Lexer as Lexer,
    LexerATNSimulator as LexerATNSimulator,
    ParserRuleContext as ParserRuleContext,
    PredictionContextCache as PredictionContextCache,
    TerminalNode as TerminalNode,
)
from antlr4.atn.ATNDeserializer import ATNDeserializer as ATNDeserializer
from antlr4.error.ErrorListener import ErrorListener as ErrorListener
from antlr4.error.Errors import LexerNoViableAltException as LexerNoViableAltException
from antlr4.InputStream import InputStream as InputStream
from antlr4.Parser import Parser as Parser
from antlr4.RuleContext import RuleContext as RuleContext
from antlr4.Token import Token as Token
from antlr4.tree.Tree import ParseTree as ParseTree
from antlr4.tree.Trees import Trees as Trees
from antlr4.xpath.XPathLexer import XPathLexer as XPathLexer

class XPath:
    WILDCARD: str
    NOT: str
    parser: Incomplete
    path: Incomplete
    elements: Incomplete
    def __init__(self, parser: Parser, path: str) -> None: ...
    def split(self, path: str): ...
    def getXPathElement(self, wordToken: Token, anywhere: bool): ...
    @staticmethod
    def findAll(tree: ParseTree, xpath: str, parser: Parser): ...
    def evaluate(self, t: ParseTree): ...

class XPathElement:
    nodeName: Incomplete
    invert: bool
    def __init__(self, nodeName: str) -> None: ...

class XPathRuleAnywhereElement(XPathElement):
    ruleIndex: Incomplete
    def __init__(self, ruleName: str, ruleIndex: int) -> None: ...
    def evaluate(self, t: ParseTree): ...

class XPathRuleElement(XPathElement):
    ruleIndex: Incomplete
    def __init__(self, ruleName: str, ruleIndex: int) -> None: ...
    def evaluate(self, t: ParseTree): ...

class XPathTokenAnywhereElement(XPathElement):
    tokenType: Incomplete
    def __init__(self, ruleName: str, tokenType: int) -> None: ...
    def evaluate(self, t: ParseTree): ...

class XPathTokenElement(XPathElement):
    tokenType: Incomplete
    def __init__(self, ruleName: str, tokenType: int) -> None: ...
    def evaluate(self, t: ParseTree): ...

class XPathWildcardAnywhereElement(XPathElement):
    def __init__(self) -> None: ...
    def evaluate(self, t: ParseTree): ...

class XPathWildcardElement(XPathElement):
    def __init__(self) -> None: ...
    def evaluate(self, t: ParseTree): ...
