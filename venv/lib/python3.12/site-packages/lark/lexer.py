# Lexer Implementation

from abc import abstractmethod, ABC
import re
from contextlib import suppress
from typing import (
    TypeVar, Type, Dict, Iterator, Collection, Callable, Optional, FrozenSet, Any,
    ClassVar, TYPE_CHECKING, overload
)
from types import ModuleType
import warnings
try:
    import interegular
except ImportError:
    pass
if TYPE_CHECKING:
    from .common import LexerConf
    from .parsers.lalr_parser_state import ParserState

from .utils import classify, get_regexp_width, Serialize, logger
from .exceptions import UnexpectedCharacters, LexError, UnexpectedToken
from .grammar import TOKEN_DEFAULT_PRIORITY


###{standalone
from copy import copy

try:  # For the standalone parser, we need to make sure that has_interegular is False to avoid NameErrors later on
    has_interegular = bool(interegular)
except NameError:
    has_interegular = False

class Pattern(Serialize, ABC):
    "An abstraction over regular expressions."

    value: str
    flags: Collection[str]
    raw: Optional[str]
    type: ClassVar[str]

    def __init__(self, value: str, flags: Collection[str] = (), raw: Optional[str] = None) -> None:
        self.value = value
        self.flags = frozenset(flags)
        self.raw = raw

    def __repr__(self):
        return repr(self.to_regexp())

    # Pattern Hashing assumes all subclasses have a different priority!
    def __hash__(self):
        return hash((type(self), self.value, self.flags))

    def __eq__(self, other):
        return type(self) == type(other) and self.value == other.value and self.flags == other.flags

    @abstractmethod
    def to_regexp(self) -> str:
        raise NotImplementedError()

    @property
    @abstractmethod
    def min_width(self) -> int:
        raise NotImplementedError()

    @property
    @abstractmethod
    def max_width(self) -> int:
        raise NotImplementedError()

    def _get_flags(self, value):
        for f in self.flags:
            value = ('(?%s:%s)' % (f, value))
        return value


class PatternStr(Pattern):
    __serialize_fields__ = 'value', 'flags', 'raw'

    type: ClassVar[str] = "str"

    def to_regexp(self) -> str:
        return self._get_flags(re.escape(self.value))

    @property
    def min_width(self) -> int:
        return len(self.value)

    @property
    def max_width(self) -> int:
        return len(self.value)


class PatternRE(Pattern):
    __serialize_fields__ = 'value', 'flags', 'raw', '_width'

    type: ClassVar[str] = "re"

    def to_regexp(self) -> str:
        return self._get_flags(self.value)

    _width = None
    def _get_width(self):
        if self._width is None:
            self._width = get_regexp_width(self.to_regexp())
        return self._width

    @property
    def min_width(self) -> int:
        return self._get_width()[0]

    @property
    def max_width(self) -> int:
        return self._get_width()[1]


class TerminalDef(Serialize):
    "A definition of a terminal"
    __serialize_fields__ = 'name', 'pattern', 'priority'
    __serialize_namespace__ = PatternStr, PatternRE

    name: str
    pattern: Pattern
    priority: int

    def __init__(self, name: str, pattern: Pattern, priority: int = TOKEN_DEFAULT_PRIORITY) -> None:
        assert isinstance(pattern, Pattern), pattern
        self.name = name
        self.pattern = pattern
        self.priority = priority

    def __repr__(self):
        return '%s(%r, %r)' % (type(self).__name__, self.name, self.pattern)

    def user_repr(self) -> str:
        if self.name.startswith('__'):  # We represent a generated terminal
            return self.pattern.raw or self.name
        else:
            return self.name

_T = TypeVar('_T', bound="Token")

class Token(str):
    """A string with meta-information, that is produced by the lexer.

    When parsing text, the resulting chunks of the input that haven't been discarded,
    will end up in the tree as Token instances. The Token class inherits from Python's ``str``,
    so normal string comparisons and operations will work as expected.

    Attributes:
        type: Name of the token (as specified in grammar)
        value: Value of the token (redundant, as ``token.value == token`` will always be true)
        start_pos: The index of the token in the text
        line: The line of the token in the text (starting with 1)
        column: The column of the token in the text (starting with 1)
        end_line: The line where the token ends
        end_column: The next column after the end of the token. For example,
            if the token is a single character with a column value of 4,
            end_column will be 5.
        end_pos: the index where the token ends (basically ``start_pos + len(token)``)
    """
    __slots__ = ('type', 'start_pos', 'value', 'line', 'column', 'end_line', 'end_column', 'end_pos')

    __match_args__ = ('type', 'value')

    type: str
    start_pos: Optional[int]
    value: Any
    line: Optional[int]
    column: Optional[int]
    end_line: Optional[int]
    end_column: Optional[int]
    end_pos: Optional[int]


    @overload
    def __new__(
            cls,
            type: str,
            value: Any,
            start_pos: Optional[int] = None,
            line: Optional[int] = None,
            column: Optional[int] = None,
            end_line: Optional[int] = None,
            end_column: Optional[int] = None,
            end_pos: Optional[int] = None
    ) -> 'Token':
        ...

    @overload
    def __new__(
            cls,
            type_: str,
            value: Any,
            start_pos: Optional[int] = None,
            line: Optional[int] = None,
            column: Optional[int] = None,
            end_line: Optional[int] = None,
            end_column: Optional[int] = None,
            end_pos: Optional[int] = None
    ) -> 'Token':        ...

    def __new__(cls, *args, **kwargs):
        if "type_" in kwargs:
            warnings.warn("`type_` is deprecated use `type` instead", DeprecationWarning)

            if "type" in kwargs:
                raise TypeError("Error: using both 'type' and the deprecated 'type_' as arguments.")
            kwargs["type"] = kwargs.pop("type_")

        return cls._future_new(*args, **kwargs)


    @classmethod
    def _future_new(cls, type, value, start_pos=None, line=None, column=None, end_line=None, end_column=None, end_pos=None):
        inst = super(Token, cls).__new__(cls, value)

        inst.type = type
        inst.start_pos = start_pos
        inst.value = value
        inst.line = line
        inst.column = column
        inst.end_line = end_line
        inst.end_column = end_column
        inst.end_pos = end_pos
        return inst

    @overload
    def update(self, type: Optional[str] = None, value: Optional[Any] = None) -> 'Token':
        ...

    @overload
    def update(self, type_: Optional[str] = None, value: Optional[Any] = None) -> 'Token':
        ...

    def update(self, *args, **kwargs):
        if "type_" in kwargs:
            warnings.warn("`type_` is deprecated use `type` instead", DeprecationWarning)

            if "type" in kwargs:
                raise TypeError("Error: using both 'type' and the deprecated 'type_' as arguments.")
            kwargs["type"] = kwargs.pop("type_")

        return self._future_update(*args, **kwargs)

    def _future_update(self, type: Optional[str] = None, value: Optional[Any] = None) -> 'Token':
        return Token.new_borrow_pos(
            type if type is not None else self.type,
            value if value is not None else self.value,
            self
        )

    @classmethod
    def new_borrow_pos(cls: Type[_T], type_: str, value: Any, borrow_t: 'Token') -> _T:
        return cls(type_, value, borrow_t.start_pos, borrow_t.line, borrow_t.column, borrow_t.end_line, borrow_t.end_column, borrow_t.end_pos)

    def __reduce__(self):
        return (self.__class__, (self.type, self.value, self.start_pos, self.line, self.column))

    def __repr__(self):
        return 'Token(%r, %r)' % (self.type, self.value)

    def __deepcopy__(self, memo):
        return Token(self.type, self.value, self.start_pos, self.line, self.column)

    def __eq__(self, other):
        if isinstance(other, Token) and self.type != other.type:
            return False

        return str.__eq__(self, other)

    __hash__ = str.__hash__


class LineCounter:
    "A utility class for keeping track of line & column information"

    __slots__ = 'char_pos', 'line', 'column', 'line_start_pos', 'newline_char'

    def __init__(self, newline_char):
        self.newline_char = newline_char
        self.char_pos = 0
        self.line = 1
        self.column = 1
        self.line_start_pos = 0

    def __eq__(self, other):
        if not isinstance(other, LineCounter):
            return NotImplemented

        return self.char_pos == other.char_pos and self.newline_char == other.newline_char

    def feed(self, token: Token, test_newline=True):
        """Consume a token and calculate the new line & column.

        As an optional optimization, set test_newline=False if token doesn't contain a newline.
        """
        if test_newline:
            newlines = token.count(self.newline_char)
            if newlines:
                self.line += newlines
                self.line_start_pos = self.char_pos + token.rindex(self.newline_char) + 1

        self.char_pos += len(token)
        self.column = self.char_pos - self.line_start_pos + 1


class UnlessCallback:
    def __init__(self, scanner):
        self.scanner = scanner

    def __call__(self, t):
        res = self.scanner.match(t.value, 0)
        if res:
            _value, t.type = res
        return t


class CallChain:
    def __init__(self, callback1, callback2, cond):
        self.callback1 = callback1
        self.callback2 = callback2
        self.cond = cond

    def __call__(self, t):
        t2 = self.callback1(t)
        return self.callback2(t) if self.cond(t2) else t2


def _get_match(re_, regexp, s, flags):
    m = re_.match(regexp, s, flags)
    if m:
        return m.group(0)

def _create_unless(terminals, g_regex_flags, re_, use_bytes):
    tokens_by_type = classify(terminals, lambda t: type(t.pattern))
    assert len(tokens_by_type) <= 2, tokens_by_type.keys()
    embedded_strs = set()
    callback = {}
    for retok in tokens_by_type.get(PatternRE, []):
        unless = []
        for strtok in tokens_by_type.get(PatternStr, []):
            if strtok.priority != retok.priority:
                continue
            s = strtok.pattern.value
            if s == _get_match(re_, retok.pattern.to_regexp(), s, g_regex_flags):
                unless.append(strtok)
                if strtok.pattern.flags <= retok.pattern.flags:
                    embedded_strs.add(strtok)
        if unless:
            callback[retok.name] = UnlessCallback(Scanner(unless, g_regex_flags, re_, match_whole=True, use_bytes=use_bytes))

    new_terminals = [t for t in terminals if t not in embedded_strs]
    return new_terminals, callback


class Scanner:
    def __init__(self, terminals, g_regex_flags, re_, use_bytes, match_whole=False):
        self.terminals = terminals
        self.g_regex_flags = g_regex_flags
        self.re_ = re_
        self.use_bytes = use_bytes
        self.match_whole = match_whole

        self.allowed_types = {t.name for t in self.terminals}

        self._mres = self._build_mres(terminals, len(terminals))

    def _build_mres(self, terminals, max_size):
        # Python sets an unreasonable group limit (currently 100) in its re module
        # Worse, the only way to know we reached it is by catching an AssertionError!
        # This function recursively tries less and less groups until it's successful.
        postfix = '$' if self.match_whole else ''
        mres = []
        while terminals:
            pattern = u'|'.join(u'(?P<%s>%s)' % (t.name, t.pattern.to_regexp() + postfix) for t in terminals[:max_size])
            if self.use_bytes:
                pattern = pattern.encode('latin-1')
            try:
                mre = self.re_.compile(pattern, self.g_regex_flags)
            except AssertionError:  # Yes, this is what Python provides us.. :/
                return self._build_mres(terminals, max_size // 2)

            mres.append(mre)
            terminals = terminals[max_size:]
        return mres

    def match(self, text, pos):
        for mre in self._mres:
            m = mre.match(text, pos)
            if m:
                return m.group(0), m.lastgroup


def _regexp_has_newline(r: str):
    r"""Expressions that may indicate newlines in a regexp:
        - newlines (\n)
        - escaped newline (\\n)
        - anything but ([^...])
        - any-char (.) when the flag (?s) exists
        - spaces (\s)
    """
    return '\n' in r or '\\n' in r or '\\s' in r or '[^' in r or ('(?s' in r and '.' in r)


class LexerState:
    """Represents the current state of the lexer as it scans the text
    (Lexer objects are only instantiated per grammar, not per text)
    """

    __slots__ = 'text', 'line_ctr', 'last_token'

    text: str
    line_ctr: LineCounter
    last_token: Optional[Token]

    def __init__(self, text: str, line_ctr: Optional[LineCounter]=None, last_token: Optional[Token]=None):
        self.text = text
        self.line_ctr = line_ctr or LineCounter(b'\n' if isinstance(text, bytes) else '\n')
        self.last_token = last_token

    def __eq__(self, other):
        if not isinstance(other, LexerState):
            return NotImplemented

        return self.text is other.text and self.line_ctr == other.line_ctr and self.last_token == other.last_token

    def __copy__(self):
        return type(self)(self.text, copy(self.line_ctr), self.last_token)


class LexerThread:
    """A thread that ties a lexer instance and a lexer state, to be used by the parser
    """

    def __init__(self, lexer: 'Lexer', lexer_state: LexerState):
        self.lexer = lexer
        self.state = lexer_state

    @classmethod
    def from_text(cls, lexer: 'Lexer', text: str) -> 'LexerThread':
        return cls(lexer, LexerState(text))

    def lex(self, parser_state):
        return self.lexer.lex(self.state, parser_state)

    def __copy__(self):
        return type(self)(self.lexer, copy(self.state))

    _Token = Token


_Callback = Callable[[Token], Token]

class Lexer(ABC):
    """Lexer interface

    Method Signatures:
        lex(self, lexer_state, parser_state) -> Iterator[Token]
    """
    @abstractmethod
    def lex(self, lexer_state: LexerState, parser_state: Any) -> Iterator[Token]:
        return NotImplemented

    def make_lexer_state(self, text):
        "Deprecated"
        return LexerState(text)


def _check_regex_collisions(terminal_to_regexp: Dict[TerminalDef, str], comparator, strict_mode, max_collisions_to_show=8):
    if not comparator:
        comparator = interegular.Comparator.from_regexes(terminal_to_regexp)

    # When in strict mode, we only ever try to provide one example, so taking
    # a long time for that should be fine
    max_time = 2 if strict_mode else 0.2

    # We don't want to show too many collisions.
    if comparator.count_marked_pairs() >= max_collisions_to_show:
        return
    for group in classify(terminal_to_regexp, lambda t: t.priority).values():
        for a, b in comparator.check(group, skip_marked=True):
            assert a.priority == b.priority
            # Mark this pair to not repeat warnings when multiple different BasicLexers see the same collision
            comparator.mark(a, b)

            # Notify the user
            message = f"Collision between Terminals {a.name} and {b.name}. "
            try:
                example = comparator.get_example_overlap(a, b, max_time).format_multiline()
            except ValueError:
                # Couldn't find an example within max_time steps.
                example = "No example could be found fast enough. However, the collision does still exists"
            if strict_mode:
                raise LexError(f"{message}\n{example}")
            logger.warning("%s The lexer will choose between them arbitrarily.\n%s", message, example)
            if comparator.count_marked_pairs() >= max_collisions_to_show:
                logger.warning("Found 8 regex collisions, will not check for more.")
                return


class AbstractBasicLexer(Lexer):
    terminals_by_name: Dict[str, TerminalDef]

    @abstractmethod
    def __init__(self, conf: 'LexerConf', comparator=None) -> None:
        ...

    @abstractmethod
    def next_token(self, lex_state: LexerState, parser_state: Any = None) -> Token:
        ...

    def lex(self, state: LexerState, parser_state: Any) -> Iterator[Token]:
        with suppress(EOFError):
            while True:
                yield self.next_token(state, parser_state)


class BasicLexer(AbstractBasicLexer):
    terminals: Collection[TerminalDef]
    ignore_types: FrozenSet[str]
    newline_types: FrozenSet[str]
    user_callbacks: Dict[str, _Callback]
    callback: Dict[str, _Callback]
    re: ModuleType

    def __init__(self, conf: 'LexerConf', comparator=None) -> None:
        terminals = list(conf.terminals)
        assert all(isinstance(t, TerminalDef) for t in terminals), terminals

        self.re = conf.re_module

        if not conf.skip_validation:
            # Sanitization
            terminal_to_regexp = {}
            for t in terminals:
                regexp = t.pattern.to_regexp()
                try:
                    self.re.compile(regexp, conf.g_regex_flags)
                except self.re.error:
                    raise LexError("Cannot compile token %s: %s" % (t.name, t.pattern))

                if t.pattern.min_width == 0:
                    raise LexError("Lexer does not allow zero-width terminals. (%s: %s)" % (t.name, t.pattern))
                if t.pattern.type == "re":
                    terminal_to_regexp[t] = regexp

            if not (set(conf.ignore) <= {t.name for t in terminals}):
                raise LexError("Ignore terminals are not defined: %s" % (set(conf.ignore) - {t.name for t in terminals}))

            if has_interegular:
                _check_regex_collisions(terminal_to_regexp, comparator, conf.strict)
            elif conf.strict:
                raise LexError("interegular must be installed for strict mode. Use `pip install 'lark[interegular]'`.")

        # Init
        self.newline_types = frozenset(t.name for t in terminals if _regexp_has_newline(t.pattern.to_regexp()))
        self.ignore_types = frozenset(conf.ignore)

        terminals.sort(key=lambda x: (-x.priority, -x.pattern.max_width, -len(x.pattern.value), x.name))
        self.terminals = terminals
        self.user_callbacks = conf.callbacks
        self.g_regex_flags = conf.g_regex_flags
        self.use_bytes = conf.use_bytes
        self.terminals_by_name = conf.terminals_by_name

        self._scanner = None

    def _build_scanner(self):
        terminals, self.callback = _create_unless(self.terminals, self.g_regex_flags, self.re, self.use_bytes)
        assert all(self.callback.values())

        for type_, f in self.user_callbacks.items():
            if type_ in self.callback:
                # Already a callback there, probably UnlessCallback
                self.callback[type_] = CallChain(self.callback[type_], f, lambda t: t.type == type_)
            else:
                self.callback[type_] = f

        self._scanner = Scanner(terminals, self.g_regex_flags, self.re, self.use_bytes)

    @property
    def scanner(self):
        if self._scanner is None:
            self._build_scanner()
        return self._scanner

    def match(self, text, pos):
        return self.scanner.match(text, pos)

    def next_token(self, lex_state: LexerState, parser_state: Any = None) -> Token:
        line_ctr = lex_state.line_ctr
        while line_ctr.char_pos < len(lex_state.text):
            res = self.match(lex_state.text, line_ctr.char_pos)
            if not res:
                allowed = self.scanner.allowed_types - self.ignore_types
                if not allowed:
                    allowed = {"<END-OF-FILE>"}
                raise UnexpectedCharacters(lex_state.text, line_ctr.char_pos, line_ctr.line, line_ctr.column,
                                           allowed=allowed, token_history=lex_state.last_token and [lex_state.last_token],
                                           state=parser_state, terminals_by_name=self.terminals_by_name)

            value, type_ = res

            ignored = type_ in self.ignore_types
            t = None
            if not ignored or type_ in self.callback:
                t = Token(type_, value, line_ctr.char_pos, line_ctr.line, line_ctr.column)
            line_ctr.feed(value, type_ in self.newline_types)
            if t is not None:
                t.end_line = line_ctr.line
                t.end_column = line_ctr.column
                t.end_pos = line_ctr.char_pos
                if t.type in self.callback:
                    t = self.callback[t.type](t)
                if not ignored:
                    if not isinstance(t, Token):
                        raise LexError("Callbacks must return a token (returned %r)" % t)
                    lex_state.last_token = t
                    return t

        # EOF
        raise EOFError(self)


class ContextualLexer(Lexer):
    lexers: Dict[int, AbstractBasicLexer]
    root_lexer: AbstractBasicLexer

    BasicLexer: Type[AbstractBasicLexer] = BasicLexer

    def __init__(self, conf: 'LexerConf', states: Dict[int, Collection[str]], always_accept: Collection[str]=()) -> None:
        terminals = list(conf.terminals)
        terminals_by_name = conf.terminals_by_name

        trad_conf = copy(conf)
        trad_conf.terminals = terminals

        if has_interegular and not conf.skip_validation:
            comparator = interegular.Comparator.from_regexes({t: t.pattern.to_regexp() for t in terminals})
        else:
            comparator = None
        lexer_by_tokens: Dict[FrozenSet[str], AbstractBasicLexer] = {}
        self.lexers = {}
        for state, accepts in states.items():
            key = frozenset(accepts)
            try:
                lexer = lexer_by_tokens[key]
            except KeyError:
                accepts = set(accepts) | set(conf.ignore) | set(always_accept)
                lexer_conf = copy(trad_conf)
                lexer_conf.terminals = [terminals_by_name[n] for n in accepts if n in terminals_by_name]
                lexer = self.BasicLexer(lexer_conf, comparator)
                lexer_by_tokens[key] = lexer

            self.lexers[state] = lexer

        assert trad_conf.terminals is terminals
        trad_conf.skip_validation = True  # We don't need to verify all terminals again
        self.root_lexer = self.BasicLexer(trad_conf, comparator)

    def lex(self, lexer_state: LexerState, parser_state: 'ParserState') -> Iterator[Token]:
        try:
            while True:
                lexer = self.lexers[parser_state.position]
                yield lexer.next_token(lexer_state, parser_state)
        except EOFError:
            pass
        except UnexpectedCharacters as e:
            # In the contextual lexer, UnexpectedCharacters can mean that the terminal is defined, but not in the current context.
            # This tests the input against the global context, to provide a nicer error.
            try:
                last_token = lexer_state.last_token  # Save last_token. Calling root_lexer.next_token will change this to the wrong token
                token = self.root_lexer.next_token(lexer_state, parser_state)
                raise UnexpectedToken(token, e.allowed, state=parser_state, token_history=[last_token], terminals_by_name=self.root_lexer.terminals_by_name)
            except UnexpectedCharacters:
                raise e  # Raise the original UnexpectedCharacters. The root lexer raises it with the wrong expected set.

###}
