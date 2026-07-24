from abc import ABC, abstractmethod
import getpass
import sys, os, pickle
import tempfile
import types
import re
from typing import (
    TypeVar, Type, List, Dict, Iterator, Callable, Union, Optional, Sequence,
    Tuple, Iterable, IO, Any, TYPE_CHECKING, Collection
)
if TYPE_CHECKING:
    from .parsers.lalr_interactive_parser import InteractiveParser
    from .tree import ParseTree
    from .visitors import Transformer
    from typing import Literal
    from .parser_frontends import ParsingFrontend

from .exceptions import ConfigurationError, assert_config, UnexpectedInput
from .utils import Serialize, SerializeMemoizer, FS, logger, TextOrSlice, LarkInput
from .load_grammar import load_grammar, FromPackageLoader, Grammar, verify_used_files, PackageResource, sha256_digest
from .tree import Tree
from .common import LexerConf, ParserConf, _ParserArgType, _LexerArgType

from .lexer import Lexer, BasicLexer, TerminalDef, LexerThread, Token
from .parse_tree_builder import ParseTreeBuilder
from .parser_frontends import _validate_frontend_args, _get_lexer_callbacks, _deserialize_parsing_frontend, _construct_parsing_frontend
from .grammar import Rule


try:
    import regex
    _has_regex = True
except ImportError:
    _has_regex = False


###{standalone


class PostLex(ABC):
    @abstractmethod
    def process(self, stream: Iterator[Token]) -> Iterator[Token]:
        return stream

    always_accept: Iterable[str] = ()

class LarkOptions(Serialize):
    """Specifies the options for Lark

    """

    start: List[str]
    debug: bool
    strict: bool
    transformer: 'Optional[Transformer]'
    propagate_positions: Union[bool, str]
    maybe_placeholders: bool
    cache: Union[bool, str]
    cache_grammar: bool
    regex: bool
    g_regex_flags: int
    keep_all_tokens: bool
    tree_class: Optional[Callable[[str, List], Any]]
    parser: _ParserArgType
    lexer: _LexerArgType
    ambiguity: 'Literal["auto", "resolve", "explicit", "forest"]'
    postlex: Optional[PostLex]
    priority: 'Optional[Literal["auto", "normal", "invert"]]'
    lexer_callbacks: Dict[str, Callable[[Token], Token]]
    use_bytes: bool
    ordered_sets: bool
    edit_terminals: Optional[Callable[[TerminalDef], TerminalDef]]
    import_paths: 'List[Union[str, Callable[[Union[None, str, PackageResource], str], Tuple[str, str]]]]'
    source_path: Optional[str]

    OPTIONS_DOC = r"""
    **===  General Options  ===**

    start
            The start symbol. Either a string, or a list of strings for multiple possible starts (Default: "start")
    debug
            Display debug information and extra warnings. Use only when debugging (Default: ``False``)
            When used with Earley, it generates a forest graph as "sppf.png", if 'dot' is installed.
    strict
            Throw an exception on any potential ambiguity, including shift/reduce conflicts, and regex collisions.
    transformer
            Applies the transformer to every parse tree (equivalent to applying it after the parse, but faster)
    propagate_positions
            Propagates positional attributes into the 'meta' attribute of all tree branches.
            Sets attributes: (line, column, end_line, end_column, start_pos, end_pos,
                              container_line, container_column, container_end_line, container_end_column)
            Accepts ``False``, ``True``, or a callable, which will filter which nodes to ignore when propagating.
    maybe_placeholders
            When ``True``, the ``[]`` operator returns ``None`` when not matched.
            When ``False``,  ``[]`` behaves like the ``?`` operator, and returns no value at all.
            (default= ``True``)
    cache
            Cache the results of the Lark grammar analysis, for x2 to x3 faster loading. LALR only for now.

            - When ``False``, does nothing (default)
            - When ``True``, caches to a temporary file in the local directory
            - When given a string, caches to the path pointed by the string
    cache_grammar
            For use with ``cache`` option. When ``True``, the unanalyzed grammar is also included in the cache.
            Useful for classes that require the ``Lark.grammar`` to be present (e.g. Reconstructor).
            (default= ``False``)
    regex
            When True, uses the ``regex`` module instead of the stdlib ``re``.
    g_regex_flags
            Flags that are applied to all terminals (both regex and strings)
    keep_all_tokens
            Prevent the tree builder from automagically removing "punctuation" tokens (Default: ``False``)
    tree_class
            Lark will produce trees comprised of instances of this class instead of the default ``lark.Tree``.

    **=== Algorithm Options ===**

    parser
            Decides which parser engine to use. Accepts "earley" or "lalr". (Default: "earley").
            (there is also a "cyk" option for legacy)
    lexer
            Decides whether or not to use a lexer stage

            - "auto" (default): Choose for me based on the parser
            - "basic": Use a basic lexer
            - "contextual": Stronger lexer (only works with parser="lalr")
            - "dynamic": Flexible and powerful (only with parser="earley")
            - "dynamic_complete": Same as dynamic, but tries *every* variation of tokenizing possible.
    ambiguity
            Decides how to handle ambiguity in the parse. Only relevant if parser="earley"

            - "resolve": The parser will automatically choose the simplest derivation
              (it chooses consistently: greedy for tokens, non-greedy for rules)
            - "explicit": The parser will return all derivations wrapped in "_ambig" tree nodes (i.e. a forest).
            - "forest": The parser will return the root of the shared packed parse forest.

    **=== Misc. / Domain Specific Options ===**

    postlex
            Lexer post-processing (Default: ``None``) Only works with the basic and contextual lexers.
    priority
            How priorities should be evaluated - "auto", ``None``, "normal", "invert" (Default: "auto")
    lexer_callbacks
            Dictionary of callbacks for the lexer. May alter tokens during lexing. Use with caution.
    use_bytes
            Accept an input of type ``bytes`` instead of ``str``.
    ordered_sets
            Should Earley use ordered-sets to achieve stable output (~10% slower than regular sets. Default: True)
    edit_terminals
            A callback for editing the terminals before parse.
    import_paths
            A List of either paths or loader functions to specify from where grammars are imported
    source_path
            Override the source of from where the grammar was loaded. Useful for relative imports and unconventional grammar loading
    **=== End of Options ===**
    """
    if __doc__:
        __doc__ += OPTIONS_DOC


    # Adding a new option needs to be done in multiple places:
    # - In the dictionary below. This is the primary truth of which options `Lark.__init__` accepts
    # - In the docstring above. It is used both for the docstring of `LarkOptions` and `Lark`, and in readthedocs
    # - As an attribute of `LarkOptions` above
    # - Potentially in `_LOAD_ALLOWED_OPTIONS` below this class, when the option doesn't change how the grammar is loaded
    # - Potentially in `lark.tools.__init__`, if it makes sense, and it can easily be passed as a cmd argument
    _defaults: Dict[str, Any] = {
        'debug': False,
        'strict': False,
        'keep_all_tokens': False,
        'tree_class': None,
        'cache': False,
        'cache_grammar': False,
        'postlex': None,
        'parser': 'earley',
        'lexer': 'auto',
        'transformer': None,
        'start': 'start',
        'priority': 'auto',
        'ambiguity': 'auto',
        'regex': False,
        'propagate_positions': False,
        'lexer_callbacks': {},
        'maybe_placeholders': True,
        'edit_terminals': None,
        'g_regex_flags': 0,
        'use_bytes': False,
        'ordered_sets': True,
        'import_paths': [],
        'source_path': None,
        '_plugins': {},
    }

    def __init__(self, options_dict: Dict[str, Any]) -> None:
        o = dict(options_dict)

        options = {}
        for name, default in self._defaults.items():
            if name in o:
                value = o.pop(name)
                if isinstance(default, bool) and name not in ('cache', 'use_bytes', 'propagate_positions'):
                    value = bool(value)
            else:
                value = default

            options[name] = value

        if isinstance(options['start'], str):
            options['start'] = [options['start']]

        self.__dict__['options'] = options


        assert_config(self.parser, ('earley', 'lalr', 'cyk', None))

        if self.parser == 'earley' and self.transformer:
            raise ConfigurationError('Cannot specify an embedded transformer when using the Earley algorithm. '
                             'Please use your transformer on the resulting parse tree, or use a different algorithm (i.e. LALR)')

        if self.cache_grammar and not self.cache:
            raise ConfigurationError('cache_grammar cannot be set when cache is disabled')

        if o:
            raise ConfigurationError("Unknown options: %s" % o.keys())

    def __getattr__(self, name: str) -> Any:
        try:
            return self.__dict__['options'][name]
        except KeyError as e:
            raise AttributeError(e)

    def __setattr__(self, name: str, value: str) -> None:
        assert_config(name, self.options.keys(), "%r isn't a valid option. Expected one of: %s")
        self.options[name] = value

    def serialize(self, memo = None) -> Dict[str, Any]:
        return self.options

    @classmethod
    def deserialize(cls, data: Dict[str, Any], memo: Dict[int, Union[TerminalDef, Rule]]) -> "LarkOptions":
        return cls(data)


# Options that can be passed to the Lark parser, even when it was loaded from cache/standalone.
# These options are only used outside of `load_grammar`.
_LOAD_ALLOWED_OPTIONS = {'postlex', 'transformer', 'lexer_callbacks', 'use_bytes', 'debug', 'g_regex_flags', 'regex', 'propagate_positions', 'tree_class', '_plugins'}

_VALID_PRIORITY_OPTIONS = ('auto', 'normal', 'invert', None)
_VALID_AMBIGUITY_OPTIONS = ('auto', 'resolve', 'explicit', 'forest')


_T = TypeVar('_T', bound="Lark")

class Lark(Serialize):
    """Main interface for the library.

    It's mostly a thin wrapper for the many different parsers, and for the tree constructor.

    Parameters:
        grammar: a string or file-object containing the grammar spec (using Lark's ebnf syntax)
        options: a dictionary controlling various aspects of Lark.

    Example:
        >>> Lark(r'''start: "foo" ''')
        Lark(...)
    """

    source_path: str
    source_grammar: str
    grammar: 'Grammar'
    options: LarkOptions
    lexer: Lexer
    parser: 'ParsingFrontend'
    terminals: Collection[TerminalDef]

    __serialize_fields__ = ['parser', 'rules', 'options']

    def __init__(self, grammar: 'Union[Grammar, str, IO[str]]', **options) -> None:
        self.options = LarkOptions(options)
        re_module: types.ModuleType

        # Update which fields are serialized
        if self.options.cache_grammar:
            self.__serialize_fields__ = self.__serialize_fields__ + ['grammar']

        # Set regex or re module
        use_regex = self.options.regex
        if use_regex:
            if _has_regex:
                re_module = regex
            else:
                raise ImportError('`regex` module must be installed if calling `Lark(regex=True)`.')
        else:
            re_module = re

        # Some, but not all file-like objects have a 'name' attribute
        if self.options.source_path is None:
            try:
                self.source_path = grammar.name  # type: ignore[union-attr]
            except AttributeError:
                self.source_path = '<string>'
        else:
            self.source_path = self.options.source_path

        # Drain file-like objects to get their contents
        try:
            read = grammar.read  # type: ignore[union-attr]
        except AttributeError:
            pass
        else:
            grammar = read()

        cache_fn = None
        cache_sha256 = None
        if isinstance(grammar, str):
            self.source_grammar = grammar
            if self.options.use_bytes:
                if not grammar.isascii():
                    raise ConfigurationError("Grammar must be ascii only, when use_bytes=True")

            if self.options.cache:
                if self.options.parser != 'lalr':
                    raise ConfigurationError("cache only works with parser='lalr' for now")

                unhashable = ('transformer', 'postlex', 'lexer_callbacks', 'edit_terminals', '_plugins')
                options_str = ''.join(k+str(v) for k, v in options.items() if k not in unhashable)
                from . import __version__
                s = grammar + options_str + __version__ + str(sys.version_info[:2])
                cache_sha256 = sha256_digest(s)

                if isinstance(self.options.cache, str):
                    cache_fn = self.options.cache
                else:
                    if self.options.cache is not True:
                        raise ConfigurationError("cache argument must be bool or str")

                    try:
                        username = getpass.getuser()
                    except Exception:
                        # The exception raised may be ImportError or OSError in
                        # the future.  For the cache, we don't care about the
                        # specific reason - we just want a username.
                        username = "unknown"


                    cache_fn = tempfile.gettempdir() + "/.lark_%s_%s_%s_%s_%s.tmp" % (
                        "cache_grammar" if self.options.cache_grammar else "cache", username, cache_sha256, *sys.version_info[:2])

                old_options = self.options
                try:
                    with FS.open(cache_fn, 'rb') as f:
                        logger.debug('Loading grammar from cache: %s', cache_fn)
                        # Remove options that aren't relevant for loading from cache
                        for name in (set(options) - _LOAD_ALLOWED_OPTIONS):
                            del options[name]
                        file_sha256 = f.readline().rstrip(b'\n')
                        cached_used_files = pickle.load(f)
                        if file_sha256 == cache_sha256.encode('utf8') and verify_used_files(cached_used_files):
                            cached_parser_data = pickle.load(f)
                            self._load(cached_parser_data, **options)
                            return
                except FileNotFoundError:
                    # The cache file doesn't exist; parse and compose the grammar as normal
                    pass
                except Exception: # We should probably narrow done which errors we catch here.
                    logger.exception("Failed to load Lark from cache: %r. We will try to carry on.", cache_fn)

                    # In theory, the Lark instance might have been messed up by the call to `_load`.
                    # In practice the only relevant thing that might have been overwritten should be `options`
                    self.options = old_options


            # Parse the grammar file and compose the grammars
            self.grammar, used_files = load_grammar(grammar, self.source_path, self.options.import_paths, self.options.keep_all_tokens)
        else:
            assert isinstance(grammar, Grammar)
            self.grammar = grammar


        if self.options.lexer == 'auto':
            if self.options.parser == 'lalr':
                self.options.lexer = 'contextual'
            elif self.options.parser == 'earley':
                if self.options.postlex is not None:
                    logger.info("postlex can't be used with the dynamic lexer, so we use 'basic' instead. "
                                "Consider using lalr with contextual instead of earley")
                    self.options.lexer = 'basic'
                else:
                    self.options.lexer = 'dynamic'
            elif self.options.parser == 'cyk':
                self.options.lexer = 'basic'
            else:
                assert False, self.options.parser
        lexer = self.options.lexer
        if isinstance(lexer, type):
            assert issubclass(lexer, Lexer)     # XXX Is this really important? Maybe just ensure interface compliance
        else:
            assert_config(lexer, ('basic', 'contextual', 'dynamic', 'dynamic_complete'))
            if self.options.postlex is not None and 'dynamic' in lexer:
                raise ConfigurationError("Can't use postlex with a dynamic lexer. Use basic or contextual instead")

        if self.options.ambiguity == 'auto':
            if self.options.parser == 'earley':
                self.options.ambiguity = 'resolve'
        else:
            assert_config(self.options.parser, ('earley', 'cyk'), "%r doesn't support disambiguation. Use one of these parsers instead: %s")

        if self.options.priority == 'auto':
            self.options.priority = 'normal'

        if self.options.priority not in _VALID_PRIORITY_OPTIONS:
            raise ConfigurationError("invalid priority option: %r. Must be one of %r" % (self.options.priority, _VALID_PRIORITY_OPTIONS))
        if self.options.ambiguity not in _VALID_AMBIGUITY_OPTIONS:
            raise ConfigurationError("invalid ambiguity option: %r. Must be one of %r" % (self.options.ambiguity, _VALID_AMBIGUITY_OPTIONS))

        if self.options.parser is None:
            terminals_to_keep = '*'     # For lexer-only mode, keep all terminals
        elif self.options.postlex is not None:
            terminals_to_keep = set(self.options.postlex.always_accept)
        else:
            terminals_to_keep = set()

        # Compile the EBNF grammar into BNF
        self.terminals, self.rules, self.ignore_tokens = self.grammar.compile(self.options.start, terminals_to_keep)

        if self.options.edit_terminals:
            for t in self.terminals:
                self.options.edit_terminals(t)

        self._terminals_dict = {t.name: t for t in self.terminals}

        # If the user asked to invert the priorities, negate them all here.
        if self.options.priority == 'invert':
            for rule in self.rules:
                if rule.options.priority is not None:
                    rule.options.priority = -rule.options.priority
            for term in self.terminals:
                term.priority = -term.priority
        # Else, if the user asked to disable priorities, strip them from the
        # rules and terminals. This allows the Earley parsers to skip an extra forest walk
        # for improved performance, if you don't need them (or didn't specify any).
        elif self.options.priority is None:
            for rule in self.rules:
                if rule.options.priority is not None:
                    rule.options.priority = None
            for term in self.terminals:
                term.priority = 0

        # TODO Deprecate lexer_callbacks?
        self.lexer_conf = LexerConf(
                self.terminals, re_module, self.ignore_tokens, self.options.postlex,
                self.options.lexer_callbacks, self.options.g_regex_flags, use_bytes=self.options.use_bytes, strict=self.options.strict
            )

        if self.options.parser:
            self.parser = self._build_parser()
        elif lexer:
            self.lexer = self._build_lexer()

        if cache_fn:
            logger.debug('Saving grammar to cache: %s', cache_fn)
            try:
                with FS.open(cache_fn, 'wb') as f:
                    assert cache_sha256 is not None
                    f.write(cache_sha256.encode('utf8') + b'\n')
                    pickle.dump(used_files, f)
                    self.save(f, _LOAD_ALLOWED_OPTIONS)
            except IOError as e:
                logger.exception("Failed to save Lark to cache: %r.", cache_fn, e)

    if __doc__:
        __doc__ += "\n\n" + LarkOptions.OPTIONS_DOC

    def _build_lexer(self, dont_ignore: bool=False) -> BasicLexer:
        lexer_conf = self.lexer_conf
        if dont_ignore:
            from copy import copy
            lexer_conf = copy(lexer_conf)
            lexer_conf.ignore = ()
        return BasicLexer(lexer_conf)

    def _prepare_callbacks(self) -> None:
        self._callbacks = {}
        # we don't need these callbacks if we aren't building a tree
        if self.options.ambiguity != 'forest':
            self._parse_tree_builder = ParseTreeBuilder(
                    self.rules,
                    self.options.tree_class or Tree,
                    self.options.propagate_positions,
                    self.options.parser != 'lalr' and self.options.ambiguity == 'explicit',
                    self.options.maybe_placeholders
                )
            self._callbacks = self._parse_tree_builder.create_callback(self.options.transformer)
        self._callbacks.update(_get_lexer_callbacks(self.options.transformer, self.terminals))

    def _build_parser(self) -> "ParsingFrontend":
        self._prepare_callbacks()
        _validate_frontend_args(self.options.parser, self.options.lexer)
        parser_conf = ParserConf(self.rules, self._callbacks, self.options.start)
        return _construct_parsing_frontend(
            self.options.parser,
            self.options.lexer,
            self.lexer_conf,
            parser_conf,
            options=self.options
        )

    def save(self, f, exclude_options: Collection[str] = ()) -> None:
        """Saves the instance into the given file object

        Useful for caching and multiprocessing.
        """
        if self.options.parser != 'lalr':
            raise NotImplementedError("Lark.save() is only implemented for the LALR(1) parser.")
        data, m = self.memo_serialize([TerminalDef, Rule])
        if exclude_options:
            data["options"] = {n: v for n, v in data["options"].items() if n not in exclude_options}
        pickle.dump({'data': data, 'memo': m}, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls: Type[_T], f) -> _T:
        """Loads an instance from the given file object

        Useful for caching and multiprocessing.
        """
        inst = cls.__new__(cls)
        return inst._load(f)

    def _deserialize_lexer_conf(self, data: Dict[str, Any], memo: Dict[int, Union[TerminalDef, Rule]], options: LarkOptions) -> LexerConf:
        lexer_conf = LexerConf.deserialize(data['lexer_conf'], memo)
        lexer_conf.callbacks = options.lexer_callbacks or {}
        lexer_conf.re_module = regex if options.regex else re
        lexer_conf.use_bytes = options.use_bytes
        lexer_conf.g_regex_flags = options.g_regex_flags
        lexer_conf.skip_validation = True
        lexer_conf.postlex = options.postlex
        return lexer_conf

    def _load(self: _T, f: Any, **kwargs) -> _T:
        if isinstance(f, dict):
            d = f
        else:
            d = pickle.load(f)
        memo_json = d['memo']
        data = d['data']

        assert memo_json
        memo = SerializeMemoizer.deserialize(memo_json, {'Rule': Rule, 'TerminalDef': TerminalDef}, {})
        if 'grammar' in data:
            self.grammar = Grammar.deserialize(data['grammar'], memo)
        options = dict(data['options'])
        if (set(kwargs) - _LOAD_ALLOWED_OPTIONS) & set(LarkOptions._defaults):
            raise ConfigurationError("Some options are not allowed when loading a Parser: {}"
                             .format(set(kwargs) - _LOAD_ALLOWED_OPTIONS))
        options.update(kwargs)
        self.options = LarkOptions.deserialize(options, memo)
        self.rules = [Rule.deserialize(r, memo) for r in data['rules']]
        self.source_path = '<deserialized>'
        _validate_frontend_args(self.options.parser, self.options.lexer)
        self.lexer_conf = self._deserialize_lexer_conf(data['parser'], memo, self.options)
        self.terminals = self.lexer_conf.terminals
        self._prepare_callbacks()
        self._terminals_dict = {t.name: t for t in self.terminals}
        self.parser = _deserialize_parsing_frontend(
            data['parser'],
            memo,
            self.lexer_conf,
            self._callbacks,
            self.options,  # Not all, but multiple attributes are used
        )
        return self

    @classmethod
    def _load_from_dict(cls, data, memo, **kwargs):
        inst = cls.__new__(cls)
        return inst._load({'data': data, 'memo': memo}, **kwargs)

    @classmethod
    def open(cls: Type[_T], grammar_filename: str, rel_to: Optional[str]=None, **options) -> _T:
        """Create an instance of Lark with the grammar given by its filename

        If ``rel_to`` is provided, the function will find the grammar filename in relation to it.

        Example:

            >>> Lark.open("grammar_file.lark", rel_to=__file__, parser="lalr")
            Lark(...)

        """
        if rel_to:
            basepath = os.path.dirname(rel_to)
            grammar_filename = os.path.join(basepath, grammar_filename)
        with open(grammar_filename, encoding='utf8') as f:
            return cls(f, **options)

    @classmethod
    def open_from_package(cls: Type[_T], package: str, grammar_path: str, search_paths: 'Sequence[str]'=[""], **options) -> _T:
        """Create an instance of Lark with the grammar loaded from within the package `package`.
        This allows grammar loading from zipapps.

        Imports in the grammar will use the `package` and `search_paths` provided, through `FromPackageLoader`

        Example:

            Lark.open_from_package(__name__, "example.lark", ("grammars",), parser=...)
        """
        package_loader = FromPackageLoader(package, search_paths)
        full_path, text = package_loader(None, grammar_path)
        options.setdefault('source_path', full_path)
        options.setdefault('import_paths', [])
        options['import_paths'].append(package_loader)
        return cls(text, **options)

    def __repr__(self):
        return 'Lark(open(%r), parser=%r, lexer=%r, ...)' % (self.source_path, self.options.parser, self.options.lexer)


    def lex(self, text: TextOrSlice, dont_ignore: bool=False) -> Iterator[Token]:
        """Only lex (and postlex) the text, without parsing it. Only relevant when lexer='basic'

        When dont_ignore=True, the lexer will return all tokens, even those marked for %ignore.

        :raises UnexpectedCharacters: In case the lexer cannot find a suitable match.
        """
        lexer: Lexer
        if not hasattr(self, 'lexer') or dont_ignore:
            lexer = self._build_lexer(dont_ignore)
        else:
            lexer = self.lexer
        lexer_thread = LexerThread.from_text(lexer, text)
        stream = lexer_thread.lex(None)
        if self.options.postlex:
            return self.options.postlex.process(stream)
        return stream

    def get_terminal(self, name: str) -> TerminalDef:
        """Get information about a terminal"""
        return self._terminals_dict[name]

    def parse_interactive(self, text: Optional[LarkInput]=None, start: Optional[str]=None) -> 'InteractiveParser':
        """Start an interactive parsing session. Only works when parser='lalr'.

        Parameters:
            text (LarkInput, optional): Text to be parsed. Required for ``resume_parse()``.
            start (str, optional): Start symbol

        Returns:
            A new InteractiveParser instance.

        See Also: ``Lark.parse()``
        """
        return self.parser.parse_interactive(text, start=start)

    def parse(self, text: LarkInput, start: Optional[str]=None, on_error: 'Optional[Callable[[UnexpectedInput], bool]]'=None) -> 'ParseTree':
        """Parse the given text, according to the options provided.

        Parameters:
            text (LarkInput): Text to be parsed, as `str` or `bytes`.
                TextSlice may also be used, but only when lexer='basic' or 'contextual'.
                If Lark was created with a custom lexer, this may be an object of any type.
            start (str, optional): Required if Lark was given multiple possible start symbols (using the start option).
            on_error (function, optional): if provided, will be called on UnexpectedInput error,
                with the exception as its argument. Return true to resume parsing, or false to raise the exception.
                LALR only. See examples/advanced/error_handling.py for an example of how to use on_error.

        Returns:
            If a transformer is supplied to ``__init__``, returns whatever is the
            result of the transformation. Otherwise, returns a Tree instance.

        :raises UnexpectedInput: On a parse error, one of these sub-exceptions will rise:
                ``UnexpectedCharacters``, ``UnexpectedToken``, or ``UnexpectedEOF``.
                For convenience, these sub-exceptions also inherit from ``ParserError`` and ``LexerError``.

        """
        if on_error is not None and self.options.parser != 'lalr':
            raise NotImplementedError("The on_error option is only implemented for the LALR(1) parser.")
        return self.parser.parse(text, start=start, on_error=on_error)


###}
