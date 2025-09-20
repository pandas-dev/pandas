from .utils import logger, NO_VALUE
from typing import Mapping, Iterable, Callable, Union, TypeVar, Tuple, Any, List, Set, Optional, Collection, TYPE_CHECKING

if TYPE_CHECKING:
    from .lexer import Token
    from .parsers.lalr_interactive_parser import InteractiveParser
    from .tree import Tree

###{standalone

class LarkError(Exception):
    pass


class ConfigurationError(LarkError, ValueError):
    pass


def assert_config(value, options: Collection, msg='Got %r, expected one of %s'):
    if value not in options:
        raise ConfigurationError(msg % (value, options))


class GrammarError(LarkError):
    pass


class ParseError(LarkError):
    pass


class LexError(LarkError):
    pass

T = TypeVar('T')

class UnexpectedInput(LarkError):
    """UnexpectedInput Error.

    Used as a base class for the following exceptions:

    - ``UnexpectedCharacters``: The lexer encountered an unexpected string
    - ``UnexpectedToken``: The parser received an unexpected token
    - ``UnexpectedEOF``: The parser expected a token, but the input ended

    After catching one of these exceptions, you may call the following helper methods to create a nicer error message.
    """
    line: int
    column: int
    pos_in_stream = None
    state: Any
    _terminals_by_name = None
    interactive_parser: 'InteractiveParser'

    def get_context(self, text: str, span: int=40) -> str:
        """Returns a pretty string pinpointing the error in the text,
        with span amount of context characters around it.

        Note:
            The parser doesn't hold a copy of the text it has to parse,
            so you have to provide it again
        """
        assert self.pos_in_stream is not None, self
        pos = self.pos_in_stream
        start = max(pos - span, 0)
        end = pos + span
        if not isinstance(text, bytes):
            before = text[start:pos].rsplit('\n', 1)[-1]
            after = text[pos:end].split('\n', 1)[0]
            return before + after + '\n' + ' ' * len(before.expandtabs()) + '^\n'
        else:
            before = text[start:pos].rsplit(b'\n', 1)[-1]
            after = text[pos:end].split(b'\n', 1)[0]
            return (before + after + b'\n' + b' ' * len(before.expandtabs()) + b'^\n').decode("ascii", "backslashreplace")

    def match_examples(self, parse_fn: 'Callable[[str], Tree]',
                             examples: Union[Mapping[T, Iterable[str]], Iterable[Tuple[T, Iterable[str]]]],
                             token_type_match_fallback: bool=False,
                             use_accepts: bool=True
                         ) -> Optional[T]:
        """Allows you to detect what's wrong in the input text by matching
        against example errors.

        Given a parser instance and a dictionary mapping some label with
        some malformed syntax examples, it'll return the label for the
        example that bests matches the current error. The function will
        iterate the dictionary until it finds a matching error, and
        return the corresponding value.

        For an example usage, see `examples/error_reporting_lalr.py`

        Parameters:
            parse_fn: parse function (usually ``lark_instance.parse``)
            examples: dictionary of ``{'example_string': value}``.
            use_accepts: Recommended to keep this as ``use_accepts=True``.
        """
        assert self.state is not None, "Not supported for this exception"

        if isinstance(examples, Mapping):
            examples = examples.items()

        candidate = (None, False)
        for i, (label, example) in enumerate(examples):
            assert not isinstance(example, str), "Expecting a list"

            for j, malformed in enumerate(example):
                try:
                    parse_fn(malformed)
                except UnexpectedInput as ut:
                    if ut.state == self.state:
                        if (
                            use_accepts
                            and isinstance(self, UnexpectedToken)
                            and isinstance(ut, UnexpectedToken)
                            and ut.accepts != self.accepts
                        ):
                            logger.debug("Different accepts with same state[%d]: %s != %s at example [%s][%s]" %
                                         (self.state, self.accepts, ut.accepts, i, j))
                            continue
                        if (
                            isinstance(self, (UnexpectedToken, UnexpectedEOF))
                            and isinstance(ut, (UnexpectedToken, UnexpectedEOF))
                        ):
                            if ut.token == self.token:  # Try exact match first
                                logger.debug("Exact Match at example [%s][%s]" % (i, j))
                                return label

                            if token_type_match_fallback:
                                # Fallback to token types match
                                if (ut.token.type == self.token.type) and not candidate[-1]:
                                    logger.debug("Token Type Fallback at example [%s][%s]" % (i, j))
                                    candidate = label, True

                        if candidate[0] is None:
                            logger.debug("Same State match at example [%s][%s]" % (i, j))
                            candidate = label, False

        return candidate[0]

    def _format_expected(self, expected):
        if self._terminals_by_name:
            d = self._terminals_by_name
            expected = [d[t_name].user_repr() if t_name in d else t_name for t_name in expected]
        return "Expected one of: \n\t* %s\n" % '\n\t* '.join(expected)


class UnexpectedEOF(ParseError, UnexpectedInput):
    """An exception that is raised by the parser, when the input ends while it still expects a token.
    """
    expected: 'List[Token]'

    def __init__(self, expected, state=None, terminals_by_name=None):
        super(UnexpectedEOF, self).__init__()

        self.expected = expected
        self.state = state
        from .lexer import Token
        self.token = Token("<EOF>", "")  # , line=-1, column=-1, pos_in_stream=-1)
        self.pos_in_stream = -1
        self.line = -1
        self.column = -1
        self._terminals_by_name = terminals_by_name


    def __str__(self):
        message = "Unexpected end-of-input. "
        message += self._format_expected(self.expected)
        return message


class UnexpectedCharacters(LexError, UnexpectedInput):
    """An exception that is raised by the lexer, when it cannot match the next
    string of characters to any of its terminals.
    """

    allowed: Set[str]
    considered_tokens: Set[Any]

    def __init__(self, seq, lex_pos, line, column, allowed=None, considered_tokens=None, state=None, token_history=None,
                 terminals_by_name=None, considered_rules=None):
        super(UnexpectedCharacters, self).__init__()

        # TODO considered_tokens and allowed can be figured out using state
        self.line = line
        self.column = column
        self.pos_in_stream = lex_pos
        self.state = state
        self._terminals_by_name = terminals_by_name

        self.allowed = allowed
        self.considered_tokens = considered_tokens
        self.considered_rules = considered_rules
        self.token_history = token_history

        if isinstance(seq, bytes):
            self.char = seq[lex_pos:lex_pos + 1].decode("ascii", "backslashreplace")
        else:
            self.char = seq[lex_pos]
        self._context = self.get_context(seq)


    def __str__(self):
        message = "No terminal matches '%s' in the current parser context, at line %d col %d" % (self.char, self.line, self.column)
        message += '\n\n' + self._context
        if self.allowed:
            message += self._format_expected(self.allowed)
        if self.token_history:
            message += '\nPrevious tokens: %s\n' % ', '.join(repr(t) for t in self.token_history)
        return message


class UnexpectedToken(ParseError, UnexpectedInput):
    """An exception that is raised by the parser, when the token it received
    doesn't match any valid step forward.

    Parameters:
        token: The mismatched token
        expected: The set of expected tokens
        considered_rules: Which rules were considered, to deduce the expected tokens
        state: A value representing the parser state. Do not rely on its value or type.
        interactive_parser: An instance of ``InteractiveParser``, that is initialized to the point of failure,
                            and can be used for debugging and error handling.

    Note: These parameters are available as attributes of the instance.
    """

    expected: Set[str]
    considered_rules: Set[str]

    def __init__(self, token, expected, considered_rules=None, state=None, interactive_parser=None, terminals_by_name=None, token_history=None):
        super(UnexpectedToken, self).__init__()

        # TODO considered_rules and expected can be figured out using state
        self.line = getattr(token, 'line', '?')
        self.column = getattr(token, 'column', '?')
        self.pos_in_stream = getattr(token, 'start_pos', None)
        self.state = state

        self.token = token
        self.expected = expected  # XXX deprecate? `accepts` is better
        self._accepts = NO_VALUE
        self.considered_rules = considered_rules
        self.interactive_parser = interactive_parser
        self._terminals_by_name = terminals_by_name
        self.token_history = token_history


    @property
    def accepts(self) -> Set[str]:
        if self._accepts is NO_VALUE:
            self._accepts = self.interactive_parser and self.interactive_parser.accepts()
        return self._accepts

    def __str__(self):
        message = ("Unexpected token %r at line %s, column %s.\n%s"
                   % (self.token, self.line, self.column, self._format_expected(self.accepts or self.expected)))
        if self.token_history:
            message += "Previous tokens: %r\n" % self.token_history

        return message



class VisitError(LarkError):
    """VisitError is raised when visitors are interrupted by an exception

    It provides the following attributes for inspection:

    Parameters:
        rule: the name of the visit rule that failed
        obj: the tree-node or token that was being processed
        orig_exc: the exception that cause it to fail

    Note: These parameters are available as attributes
    """

    obj: 'Union[Tree, Token]'
    orig_exc: Exception

    def __init__(self, rule, obj, orig_exc):
        message = 'Error trying to process rule "%s":\n\n%s' % (rule, orig_exc)
        super(VisitError, self).__init__(message)

        self.rule = rule
        self.obj = obj
        self.orig_exc = orig_exc


class MissingVariableError(LarkError):
    pass

###}
