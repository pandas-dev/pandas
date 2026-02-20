"Provides a post-lexer for implementing Python-style indentation."

from abc import ABC, abstractmethod
from typing import List, Iterator

from .exceptions import LarkError
from .lark import PostLex
from .lexer import Token

###{standalone

class DedentError(LarkError):
    pass

class Indenter(PostLex, ABC):
    """This is a postlexer that "injects" indent/dedent tokens based on indentation.

    It keeps track of the current indentation, as well as the current level of parentheses.
    Inside parentheses, the indentation is ignored, and no indent/dedent tokens get generated.

    Note: This is an abstract class. To use it, inherit and implement all its abstract methods:
        - tab_len
        - NL_type
        - OPEN_PAREN_types, CLOSE_PAREN_types
        - INDENT_type, DEDENT_type

    See also: the ``postlex`` option in `Lark`.
    """
    paren_level: int
    indent_level: List[int]

    def __init__(self) -> None:
        self.paren_level = 0
        self.indent_level = [0]
        assert self.tab_len > 0

    def handle_NL(self, token: Token) -> Iterator[Token]:
        if self.paren_level > 0:
            return

        yield token

        indent_str = token.rsplit('\n', 1)[1] # Tabs and spaces
        indent = indent_str.count(' ') + indent_str.count('\t') * self.tab_len

        if indent > self.indent_level[-1]:
            self.indent_level.append(indent)
            yield Token.new_borrow_pos(self.INDENT_type, indent_str, token)
        else:
            while indent < self.indent_level[-1]:
                self.indent_level.pop()
                yield Token.new_borrow_pos(self.DEDENT_type, indent_str, token)

            if indent != self.indent_level[-1]:
                raise DedentError('Unexpected dedent to column %s. Expected dedent to %s' % (indent, self.indent_level[-1]))

    def _process(self, stream):
        token = None
        for token in stream:
            if token.type == self.NL_type:
                yield from self.handle_NL(token)
            else:
                yield token

            if token.type in self.OPEN_PAREN_types:
                self.paren_level += 1
            elif token.type in self.CLOSE_PAREN_types:
                self.paren_level -= 1
                assert self.paren_level >= 0

        while len(self.indent_level) > 1:
            self.indent_level.pop()
            yield Token.new_borrow_pos(self.DEDENT_type, '', token) if token else Token(self.DEDENT_type, '', 0, 0, 0, 0, 0, 0)

        assert self.indent_level == [0], self.indent_level

    def process(self, stream):
        self.paren_level = 0
        self.indent_level = [0]
        return self._process(stream)

    # XXX Hack for ContextualLexer. Maybe there's a more elegant solution?
    @property
    def always_accept(self):
        return (self.NL_type,)

    @property
    @abstractmethod
    def NL_type(self) -> str:
        "The name of the newline token"
        raise NotImplementedError()

    @property
    @abstractmethod
    def OPEN_PAREN_types(self) -> List[str]:
        "The names of the tokens that open a parenthesis"
        raise NotImplementedError()

    @property
    @abstractmethod
    def CLOSE_PAREN_types(self) -> List[str]:
        """The names of the tokens that close a parenthesis
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def INDENT_type(self) -> str:
        """The name of the token that starts an indentation in the grammar.

        See also: %declare
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def DEDENT_type(self) -> str:
        """The name of the token that end an indentation in the grammar.

        See also: %declare
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def tab_len(self) -> int:
        """How many spaces does a tab equal"""
        raise NotImplementedError()


class PythonIndenter(Indenter):
    """A postlexer that "injects" _INDENT/_DEDENT tokens based on indentation, according to the Python syntax.

    See also: the ``postlex`` option in `Lark`.
    """

    NL_type = '_NEWLINE'
    OPEN_PAREN_types = ['LPAR', 'LSQB', 'LBRACE']
    CLOSE_PAREN_types = ['RPAR', 'RSQB', 'RBRACE']
    INDENT_type = '_INDENT'
    DEDENT_type = '_DEDENT'
    tab_len = 8

###}
