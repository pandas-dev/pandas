"""Parses and compiles Lark grammars into an internal representation.
"""

import hashlib
import os.path
import sys
from collections import namedtuple
from copy import copy, deepcopy
import pkgutil
from ast import literal_eval
from contextlib import suppress
from typing import List, Tuple, Union, Callable, Dict, Optional, Sequence, Generator

from .utils import bfs, logger, classify_bool, is_id_continue, is_id_start, bfs_all_unique, small_factors, OrderedSet, Serialize
from .lexer import Token, TerminalDef, PatternStr, PatternRE, Pattern

from .parse_tree_builder import ParseTreeBuilder
from .parser_frontends import ParsingFrontend
from .common import LexerConf, ParserConf
from .grammar import RuleOptions, Rule, Terminal, NonTerminal, Symbol, TOKEN_DEFAULT_PRIORITY
from .utils import classify, dedup_list
from .exceptions import GrammarError, UnexpectedCharacters, UnexpectedToken, ParseError, UnexpectedInput

from .tree import Tree, SlottedTree as ST
from .visitors import Transformer, Visitor, v_args, Transformer_InPlace, Transformer_NonRecursive
inline_args = v_args(inline=True)

IMPORT_PATHS = ['grammars']

EXT = '.lark'

_RE_FLAGS = 'imslux'

_EMPTY = Symbol('__empty__')

_TERMINAL_NAMES = {
    '.' : 'DOT',
    ',' : 'COMMA',
    ':' : 'COLON',
    ';' : 'SEMICOLON',
    '+' : 'PLUS',
    '-' : 'MINUS',
    '*' : 'STAR',
    '/' : 'SLASH',
    '\\' : 'BACKSLASH',
    '|' : 'VBAR',
    '?' : 'QMARK',
    '!' : 'BANG',
    '@' : 'AT',
    '#' : 'HASH',
    '$' : 'DOLLAR',
    '%' : 'PERCENT',
    '^' : 'CIRCUMFLEX',
    '&' : 'AMPERSAND',
    '_' : 'UNDERSCORE',
    '<' : 'LESSTHAN',
    '>' : 'MORETHAN',
    '=' : 'EQUAL',
    '"' : 'DBLQUOTE',
    '\'' : 'QUOTE',
    '`' : 'BACKQUOTE',
    '~' : 'TILDE',
    '(' : 'LPAR',
    ')' : 'RPAR',
    '{' : 'LBRACE',
    '}' : 'RBRACE',
    '[' : 'LSQB',
    ']' : 'RSQB',
    '\n' : 'NEWLINE',
    '\r\n' : 'CRLF',
    '\t' : 'TAB',
    ' ' : 'SPACE',
}

# Grammar Parser
TERMINALS = {
    '_LPAR': r'\(',
    '_RPAR': r'\)',
    '_LBRA': r'\[',
    '_RBRA': r'\]',
    '_LBRACE': r'\{',
    '_RBRACE': r'\}',
    'OP': '[+*]|[?](?![a-z_])',
    '_COLON': ':',
    '_COMMA': ',',
    '_OR': r'\|',
    '_DOT': r'\.(?!\.)',
    '_DOTDOT': r'\.\.',
    'TILDE': '~',
    'RULE_MODIFIERS': '(!|![?]?|[?]!?)(?=[_a-z])',
    'RULE': '_?[a-z][_a-z0-9]*',
    'TERMINAL': '_?[A-Z][_A-Z0-9]*',
    'STRING': r'"(\\"|\\\\|[^"\n])*?"i?',
    'REGEXP': r'/(?!/)(\\/|\\\\|[^/])*?/[%s]*' % _RE_FLAGS,
    '_NL': r'(\r?\n)+\s*',
    '_NL_OR': r'(\r?\n)+\s*\|',
    'WS': r'[ \t]+',
    'COMMENT': r'\s*//[^\n]*|\s*#[^\n]*',
    'BACKSLASH': r'\\[ ]*\n',
    '_TO': '->',
    '_IGNORE': r'%ignore',
    '_OVERRIDE': r'%override',
    '_DECLARE': r'%declare',
    '_EXTEND': r'%extend',
    '_IMPORT': r'%import',
    'NUMBER': r'[+-]?\d+',
}

RULES = {
    'start': ['_list'],
    '_list':  ['_item', '_list _item'],
    '_item':  ['rule', 'term', 'ignore', 'import', 'declare', 'override', 'extend', '_NL'],

    'rule': ['rule_modifiers RULE template_params priority _COLON expansions _NL'],
    'rule_modifiers': ['RULE_MODIFIERS',
                       ''],
    'priority': ['_DOT NUMBER',
                 ''],
    'template_params': ['_LBRACE _template_params _RBRACE',
                        ''],
    '_template_params': ['RULE',
                         '_template_params _COMMA RULE'],
    'expansions': ['_expansions'],
    '_expansions': ['alias',
                    '_expansions _OR alias',
                    '_expansions _NL_OR alias'],

    '?alias':     ['expansion _TO nonterminal', 'expansion'],
    'expansion': ['_expansion'],

    '_expansion': ['', '_expansion expr'],

    '?expr': ['atom',
              'atom OP',
              'atom TILDE NUMBER',
              'atom TILDE NUMBER _DOTDOT NUMBER',
              ],

    '?atom': ['_LPAR expansions _RPAR',
              'maybe',
              'value'],

    'value': ['terminal',
              'nonterminal',
              'literal',
              'range',
              'template_usage'],

    'terminal': ['TERMINAL'],
    'nonterminal': ['RULE'],

    '?name': ['RULE', 'TERMINAL'],
    '?symbol': ['terminal', 'nonterminal'],

    'maybe': ['_LBRA expansions _RBRA'],
    'range': ['STRING _DOTDOT STRING'],

    'template_usage': ['nonterminal _LBRACE _template_args _RBRACE'],
    '_template_args': ['value',
                       '_template_args _COMMA value'],

    'term': ['TERMINAL _COLON expansions _NL',
             'TERMINAL _DOT NUMBER _COLON expansions _NL'],
    'override': ['_OVERRIDE rule',
                 '_OVERRIDE term'],
    'extend': ['_EXTEND rule',
               '_EXTEND term'],
    'ignore': ['_IGNORE expansions _NL'],
    'declare': ['_DECLARE _declare_args _NL'],
    'import': ['_IMPORT _import_path _NL',
               '_IMPORT _import_path _LPAR name_list _RPAR _NL',
               '_IMPORT _import_path _TO name _NL'],

    '_import_path': ['import_lib', 'import_rel'],
    'import_lib': ['_import_args'],
    'import_rel': ['_DOT _import_args'],
    '_import_args': ['name', '_import_args _DOT name'],

    'name_list': ['_name_list'],
    '_name_list': ['name', '_name_list _COMMA name'],

    '_declare_args': ['symbol', '_declare_args symbol'],
    'literal': ['REGEXP', 'STRING'],
}


# Value 5 keeps the number of states in the lalr parser somewhat minimal
# It isn't optimal, but close to it. See PR #949
SMALL_FACTOR_THRESHOLD = 5
# The Threshold whether repeat via ~ are split up into different rules
# 50 is chosen since it keeps the number of states low and therefore lalr analysis time low,
# while not being to overaggressive and unnecessarily creating rules that might create shift/reduce conflicts.
# (See PR #949)
REPEAT_BREAK_THRESHOLD = 50


class FindRuleSize(Transformer):
    def __init__(self, keep_all_tokens: bool):
        self.keep_all_tokens = keep_all_tokens

    def _will_not_get_removed(self, sym: Symbol) -> bool:
        if isinstance(sym, NonTerminal):
            return not sym.name.startswith('_')
        if isinstance(sym, Terminal):
            return self.keep_all_tokens or not sym.filter_out
        if sym is _EMPTY:
            return False
        assert False, sym

    def _args_as_int(self, args: List[Union[int, Symbol]]) -> Generator[int, None, None]:
        for a in args:
            if isinstance(a, int):
                yield a
            elif isinstance(a, Symbol):
                yield 1 if self._will_not_get_removed(a) else 0
            else:
                assert False

    def expansion(self, args) -> int:
        return sum(self._args_as_int(args))

    def expansions(self, args) -> int:
        return max(self._args_as_int(args))


@inline_args
class EBNF_to_BNF(Transformer_InPlace):
    def __init__(self):
        self.new_rules = []
        self.rules_cache = {}
        self.prefix = 'anon'
        self.i = 0
        self.rule_options = None

    def _name_rule(self, inner: str):
        new_name = '__%s_%s_%d' % (self.prefix, inner, self.i)
        self.i += 1
        return new_name

    def _add_rule(self, key, name, expansions):
        t = NonTerminal(name)
        self.new_rules.append((name, expansions, self.rule_options))
        self.rules_cache[key] = t
        return t

    def _add_recurse_rule(self, type_: str, expr: Tree):
        try:
            return self.rules_cache[expr]
        except KeyError:
            new_name = self._name_rule(type_)
            t = NonTerminal(new_name)
            tree = ST('expansions', [
                ST('expansion', [expr]),
                ST('expansion', [t, expr])
            ])
            return self._add_rule(expr, new_name, tree)

    def _add_repeat_rule(self, a, b, target, atom):
        """Generate a rule that repeats target ``a`` times, and repeats atom ``b`` times.

        When called recursively (into target), it repeats atom for x(n) times, where:
            x(0) = 1
            x(n) = a(n) * x(n-1) + b

        Example rule when a=3, b=4:

            new_rule: target target target atom atom atom atom

        """
        key = (a, b, target, atom)
        try:
            return self.rules_cache[key]
        except KeyError:
            new_name = self._name_rule('repeat_a%d_b%d' % (a, b))
            tree = ST('expansions', [ST('expansion', [target] * a + [atom] * b)])
            return self._add_rule(key, new_name, tree)

    def _add_repeat_opt_rule(self, a, b, target, target_opt, atom):
        """Creates a rule that matches atom 0 to (a*n+b)-1 times.

        When target matches n times atom, and target_opt 0 to n-1 times target_opt,

        First we generate target * i followed by target_opt, for i from 0 to a-1
        These match 0 to n*a - 1 times atom

        Then we generate target * a followed by atom * i, for i from 0 to b-1
        These match n*a to n*a + b-1 times atom

        The created rule will not have any shift/reduce conflicts so that it can be used with lalr

        Example rule when a=3, b=4:

            new_rule: target_opt
                    | target target_opt
                    | target target target_opt

                    | target target target
                    | target target target atom
                    | target target target atom atom
                    | target target target atom atom atom

        """
        key = (a, b, target, atom, "opt")
        try:
            return self.rules_cache[key]
        except KeyError:
            new_name = self._name_rule('repeat_a%d_b%d_opt' % (a, b))
            tree = ST('expansions', [
                ST('expansion', [target]*i + [target_opt]) for i in range(a)
            ] + [
                ST('expansion', [target]*a + [atom]*i) for i in range(b)
            ])
            return self._add_rule(key, new_name, tree)

    def _generate_repeats(self, rule: Tree, mn: int, mx: int):
        """Generates a rule tree that repeats ``rule`` exactly between ``mn`` to ``mx`` times.
        """
        # For a small number of repeats, we can take the naive approach
        if mx < REPEAT_BREAK_THRESHOLD:
            return ST('expansions', [ST('expansion', [rule] * n) for n in range(mn, mx + 1)])

        # For large repeat values, we break the repetition into sub-rules.
        # We treat ``rule~mn..mx`` as ``rule~mn rule~0..(diff=mx-mn)``.
        # We then use small_factors to split up mn and diff up into values [(a, b), ...]
        # This values are used with the help of _add_repeat_rule and _add_repeat_rule_opt
        # to generate a complete rule/expression that matches the corresponding number of repeats
        mn_target = rule
        for a, b in small_factors(mn, SMALL_FACTOR_THRESHOLD):
            mn_target = self._add_repeat_rule(a, b, mn_target, rule)
        if mx == mn:
            return mn_target

        diff = mx - mn + 1  # We add one because _add_repeat_opt_rule generates rules that match one less
        diff_factors = small_factors(diff, SMALL_FACTOR_THRESHOLD)
        diff_target = rule  # Match rule 1 times
        diff_opt_target = ST('expansion', [])  # match rule 0 times (e.g. up to 1 -1 times)
        for a, b in diff_factors[:-1]:
            diff_opt_target = self._add_repeat_opt_rule(a, b, diff_target, diff_opt_target, rule)
            diff_target = self._add_repeat_rule(a, b, diff_target, rule)

        a, b = diff_factors[-1]
        diff_opt_target = self._add_repeat_opt_rule(a, b, diff_target, diff_opt_target, rule)

        return ST('expansions', [ST('expansion', [mn_target] + [diff_opt_target])])

    def expr(self, rule: Tree, op: Token, *args):
        if op.value == '?':
            empty = ST('expansion', [])
            return ST('expansions', [rule, empty])
        elif op.value == '+':
            # a : b c+ d
            #   -->
            # a : b _c d
            # _c : _c c | c;
            return self._add_recurse_rule('plus', rule)
        elif op.value == '*':
            # a : b c* d
            #   -->
            # a : b _c? d
            # _c : _c c | c;
            new_name = self._add_recurse_rule('star', rule)
            return ST('expansions', [new_name, ST('expansion', [])])
        elif op.value == '~':
            if len(args) == 1:
                mn = mx = int(args[0])
            else:
                mn, mx = map(int, args)
                if mx < mn or mn < 0:
                    raise GrammarError("Bad Range for %s (%d..%d isn't allowed)" % (rule, mn, mx))

            return self._generate_repeats(rule, mn, mx)

        assert False, op

    def maybe(self, rule: Tree):
        keep_all_tokens = self.rule_options and self.rule_options.keep_all_tokens
        rule_size = FindRuleSize(keep_all_tokens).transform(rule)
        empty = ST('expansion', [_EMPTY] * rule_size)
        return ST('expansions', [rule, empty])


class SimplifyRule_Visitor(Visitor):

    @staticmethod
    def _flatten(tree: Tree):
        while tree.expand_kids_by_data(tree.data):
            pass

    def expansion(self, tree: Tree):
        # rules_list unpacking
        # a : b (c|d) e
        #  -->
        # a : b c e | b d e
        #
        # In AST terms:
        # expansion(b, expansions(c, d), e)
        #   -->
        # expansions( expansion(b, c, e), expansion(b, d, e) )

        self._flatten(tree)

        for i, child in enumerate(tree.children):
            if isinstance(child, Tree) and child.data == 'expansions':
                tree.data = 'expansions'
                tree.children = [self.visit(ST('expansion', [option if i == j else other
                                                             for j, other in enumerate(tree.children)]))
                                 for option in dedup_list(child.children)]
                self._flatten(tree)
                break

    def alias(self, tree):
        rule, alias_name = tree.children
        if rule.data == 'expansions':
            aliases = []
            for child in tree.children[0].children:
                aliases.append(ST('alias', [child, alias_name]))
            tree.data = 'expansions'
            tree.children = aliases

    def expansions(self, tree: Tree):
        self._flatten(tree)
        # Ensure all children are unique
        if len(set(tree.children)) != len(tree.children):
            tree.children = dedup_list(tree.children)   # dedup is expensive, so try to minimize its use


class RuleTreeToText(Transformer):
    def expansions(self, x):
        return x

    def expansion(self, symbols):
        return symbols, None

    def alias(self, x):
        (expansion, _alias), alias = x
        assert _alias is None, (alias, expansion, '-', _alias)  # Double alias not allowed
        return expansion, alias.name


class PrepareAnonTerminals(Transformer_InPlace):
    """Create a unique list of anonymous terminals. Attempt to give meaningful names to them when we add them"""

    def __init__(self, terminals):
        self.terminals = terminals
        self.term_set = {td.name for td in self.terminals}
        self.term_reverse = {td.pattern: td for td in terminals}
        self.i = 0
        self.rule_options = None

    @inline_args
    def pattern(self, p):
        value = p.value
        if p in self.term_reverse and p.flags != self.term_reverse[p].pattern.flags:
            raise GrammarError(u'Conflicting flags for the same terminal: %s' % p)

        term_name = None

        if isinstance(p, PatternStr):
            try:
                # If already defined, use the user-defined terminal name
                term_name = self.term_reverse[p].name
            except KeyError:
                # Try to assign an indicative anon-terminal name
                try:
                    term_name = _TERMINAL_NAMES[value]
                except KeyError:
                    if value and is_id_continue(value) and is_id_start(value[0]) and value.upper() not in self.term_set:
                        term_name = value.upper()

                if term_name in self.term_set:
                    term_name = None

        elif isinstance(p, PatternRE):
            if p in self.term_reverse:  # Kind of a weird placement.name
                term_name = self.term_reverse[p].name
        else:
            assert False, p

        if term_name is None:
            term_name = '__ANON_%d' % self.i
            self.i += 1

        if term_name not in self.term_set:
            assert p not in self.term_reverse
            self.term_set.add(term_name)
            termdef = TerminalDef(term_name, p)
            self.term_reverse[p] = termdef
            self.terminals.append(termdef)

        filter_out = False if self.rule_options and self.rule_options.keep_all_tokens else isinstance(p, PatternStr)

        return Terminal(term_name, filter_out=filter_out)


class _ReplaceSymbols(Transformer_InPlace):
    """Helper for ApplyTemplates"""

    def __init__(self):
        self.names = {}

    def value(self, c):
        if len(c) == 1 and isinstance(c[0], Symbol) and c[0].name in self.names:
            return self.names[c[0].name]
        return self.__default__('value', c, None)

    def template_usage(self, c):
        name = c[0].name
        if name in self.names:
            return self.__default__('template_usage', [self.names[name]] + c[1:], None)
        return self.__default__('template_usage', c, None)


class ApplyTemplates(Transformer_InPlace):
    """Apply the templates, creating new rules that represent the used templates"""

    def __init__(self, rule_defs):
        self.rule_defs = rule_defs
        self.replacer = _ReplaceSymbols()
        self.created_templates = set()

    def template_usage(self, c):
        name = c[0].name
        args = c[1:]
        result_name = "%s{%s}" % (name, ",".join(a.name for a in args))
        if result_name not in self.created_templates:
            self.created_templates.add(result_name)
            (_n, params, tree, options) ,= (t for t in self.rule_defs if t[0] == name)
            assert len(params) == len(args), args
            result_tree = deepcopy(tree)
            self.replacer.names = dict(zip(params, args))
            self.replacer.transform(result_tree)
            self.rule_defs.append((result_name, [], result_tree, deepcopy(options)))
        return NonTerminal(result_name)


def _rfind(s, choices):
    return max(s.rfind(c) for c in choices)


def eval_escaping(s):
    w = ''
    i = iter(s)
    for n in i:
        w += n
        if n == '\\':
            try:
                n2 = next(i)
            except StopIteration:
                raise GrammarError("Literal ended unexpectedly (bad escaping): `%r`" % s)
            if n2 == '\\':
                w += '\\\\'
            elif n2 not in 'Uuxnftr':
                w += '\\'
            w += n2
    w = w.replace('\\"', '"').replace("'", "\\'")

    to_eval = "u'''%s'''" % w
    try:
        s = literal_eval(to_eval)
    except SyntaxError as e:
        raise GrammarError(s, e)

    return s


def _literal_to_pattern(literal):
    assert isinstance(literal, Token)
    v = literal.value
    flag_start = _rfind(v, '/"')+1
    assert flag_start > 0
    flags = v[flag_start:]
    assert all(f in _RE_FLAGS for f in flags), flags

    if literal.type == 'STRING' and '\n' in v:
        raise GrammarError('You cannot put newlines in string literals')

    if literal.type == 'REGEXP' and '\n' in v and 'x' not in flags:
        raise GrammarError('You can only use newlines in regular expressions '
                           'with the `x` (verbose) flag')

    v = v[:flag_start]
    assert v[0] == v[-1] and v[0] in '"/'
    x = v[1:-1]

    s = eval_escaping(x)

    if s == "":
        raise GrammarError("Empty terminals are not allowed (%s)" % literal)

    if literal.type == 'STRING':
        s = s.replace('\\\\', '\\')
        return PatternStr(s, flags, raw=literal.value)
    elif literal.type == 'REGEXP':
        return PatternRE(s, flags, raw=literal.value)
    else:
        assert False, 'Invariant failed: literal.type not in ["STRING", "REGEXP"]'


@inline_args
class PrepareLiterals(Transformer_InPlace):
    def literal(self, literal):
        return ST('pattern', [_literal_to_pattern(literal)])

    def range(self, start, end):
        assert start.type == end.type == 'STRING'
        start = start.value[1:-1]
        end = end.value[1:-1]
        assert len(eval_escaping(start)) == len(eval_escaping(end)) == 1
        regexp = '[%s-%s]' % (start, end)
        return ST('pattern', [PatternRE(regexp)])


def _make_joined_pattern(regexp, flags_set) -> PatternRE:
    return PatternRE(regexp, ())

class TerminalTreeToPattern(Transformer_NonRecursive):
    def pattern(self, ps):
        p ,= ps
        return p

    def expansion(self, items: List[Pattern]) -> Pattern:
        if not items:
            return PatternStr('')

        if len(items) == 1:
            return items[0]

        pattern = ''.join(i.to_regexp() for i in items)
        return _make_joined_pattern(pattern, {i.flags for i in items})

    def expansions(self, exps: List[Pattern]) -> Pattern:
        if len(exps) == 1:
            return exps[0]

        # Do a bit of sorting to make sure that the longest option is returned
        # (Python's re module otherwise prefers just 'l' when given (l|ll) and both could match)
        exps.sort(key=lambda x: (-x.max_width, -x.min_width, -len(x.value)))

        pattern = '(?:%s)' % ('|'.join(i.to_regexp() for i in exps))
        return _make_joined_pattern(pattern, {i.flags for i in exps})

    def expr(self, args) -> Pattern:
        inner: Pattern
        inner, op = args[:2]
        if op == '~':
            if len(args) == 3:
                op = "{%d}" % int(args[2])
            else:
                mn, mx = map(int, args[2:])
                if mx < mn:
                    raise GrammarError("Bad Range for %s (%d..%d isn't allowed)" % (inner, mn, mx))
                op = "{%d,%d}" % (mn, mx)
        else:
            assert len(args) == 2
        return PatternRE('(?:%s)%s' % (inner.to_regexp(), op), inner.flags)

    def maybe(self, expr):
        return self.expr(expr + ['?'])

    def alias(self, t):
        raise GrammarError("Aliasing not allowed in terminals (You used -> in the wrong place)")

    def value(self, v):
        return v[0]


class ValidateSymbols(Transformer_InPlace):
    def value(self, v):
        v ,= v
        assert isinstance(v, (Tree, Symbol))
        return v


def nr_deepcopy_tree(t):
    """Deepcopy tree `t` without recursion"""
    return Transformer_NonRecursive(False).transform(t)


class Grammar(Serialize):

    term_defs: List[Tuple[str, Tuple[Tree, int]]]
    rule_defs: List[Tuple[str, Tuple[str, ...], Tree, RuleOptions]]
    ignore: List[str]

    def __init__(self, rule_defs: List[Tuple[str, Tuple[str, ...], Tree, RuleOptions]], term_defs: List[Tuple[str, Tuple[Tree, int]]], ignore: List[str]) -> None:
        self.term_defs = term_defs
        self.rule_defs = rule_defs
        self.ignore = ignore

    __serialize_fields__ = 'term_defs', 'rule_defs', 'ignore'

    def compile(self, start, terminals_to_keep) -> Tuple[List[TerminalDef], List[Rule], List[str]]:
        # We change the trees in-place (to support huge grammars)
        # So deepcopy allows calling compile more than once.
        term_defs = [(n, (nr_deepcopy_tree(t), p)) for n, (t, p) in self.term_defs]
        rule_defs = [(n, p, nr_deepcopy_tree(t), o) for n, p, t, o in self.rule_defs]

        # ===================
        #  Compile Terminals
        # ===================

        # Convert terminal-trees to strings/regexps

        for name, (term_tree, priority) in term_defs:
            if term_tree is None:  # Terminal added through %declare
                continue
            expansions = list(term_tree.find_data('expansion'))
            if len(expansions) == 1 and not expansions[0].children:
                raise GrammarError("Terminals cannot be empty (%s)" % name)

        transformer = PrepareLiterals() * TerminalTreeToPattern()
        terminals = [TerminalDef(name, transformer.transform(term_tree), priority)
                     for name, (term_tree, priority) in term_defs if term_tree]

        # =================
        #  Compile Rules
        # =================

        # 1. Pre-process terminals
        anon_tokens_transf = PrepareAnonTerminals(terminals)
        transformer = PrepareLiterals() * ValidateSymbols() * anon_tokens_transf  # Adds to terminals

        # 2. Inline Templates

        transformer *= ApplyTemplates(rule_defs)

        # 3. Convert EBNF to BNF (and apply step 1 & 2)
        ebnf_to_bnf = EBNF_to_BNF()
        rules = []
        i = 0
        while i < len(rule_defs):  # We have to do it like this because rule_defs might grow due to templates
            name, params, rule_tree, options = rule_defs[i]
            i += 1
            if len(params) != 0:  # Dont transform templates
                continue
            rule_options = RuleOptions(keep_all_tokens=True) if options and options.keep_all_tokens else None
            ebnf_to_bnf.rule_options = rule_options
            ebnf_to_bnf.prefix = name
            anon_tokens_transf.rule_options = rule_options
            tree = transformer.transform(rule_tree)
            res: Tree = ebnf_to_bnf.transform(tree)
            rules.append((name, res, options))
        rules += ebnf_to_bnf.new_rules

        assert len(rules) == len({name for name, _t, _o in rules}), "Whoops, name collision"

        # 4. Compile tree to Rule objects
        rule_tree_to_text = RuleTreeToText()

        simplify_rule = SimplifyRule_Visitor()
        compiled_rules: List[Rule] = []
        for rule_content in rules:
            name, tree, options = rule_content
            simplify_rule.visit(tree)
            expansions = rule_tree_to_text.transform(tree)

            for i, (expansion, alias) in enumerate(expansions):
                if alias and name.startswith('_'):
                    raise GrammarError("Rule %s is marked for expansion (it starts with an underscore) and isn't allowed to have aliases (alias=%s)"% (name, alias))

                empty_indices = tuple(x==_EMPTY for x in expansion)
                if any(empty_indices):
                    exp_options = copy(options) or RuleOptions()
                    exp_options.empty_indices = empty_indices
                    expansion = [x for x in expansion if x!=_EMPTY]
                else:
                    exp_options = options

                for sym in expansion:
                    assert isinstance(sym, Symbol)
                    if sym.is_term and exp_options and exp_options.keep_all_tokens:
                        assert isinstance(sym, Terminal)
                        sym.filter_out = False
                rule = Rule(NonTerminal(name), expansion, i, alias, exp_options)
                compiled_rules.append(rule)

        # Remove duplicates of empty rules, throw error for non-empty duplicates
        if len(set(compiled_rules)) != len(compiled_rules):
            duplicates = classify(compiled_rules, lambda x: x)
            for dups in duplicates.values():
                if len(dups) > 1:
                    if dups[0].expansion:
                        raise GrammarError("Rules defined twice: %s\n\n(Might happen due to colliding expansion of optionals: [] or ?)"
                                           % ''.join('\n  * %s' % i for i in dups))

                    # Empty rule; assert all other attributes are equal
                    assert len({(r.alias, r.order, r.options) for r in dups}) == len(dups)

            # Remove duplicates
            compiled_rules = list(OrderedSet(compiled_rules))

        # Filter out unused rules
        while True:
            c = len(compiled_rules)
            used_rules = {s for r in compiled_rules
                            for s in r.expansion
                            if isinstance(s, NonTerminal)
                            and s != r.origin}
            used_rules |= {NonTerminal(s) for s in start}
            compiled_rules, unused = classify_bool(compiled_rules, lambda r: r.origin in used_rules)
            for r in unused:
                logger.debug("Unused rule: %s", r)
            if len(compiled_rules) == c:
                break

        # Filter out unused terminals
        if terminals_to_keep != '*':
            used_terms = {t.name for r in compiled_rules
                                 for t in r.expansion
                                 if isinstance(t, Terminal)}
            terminals, unused = classify_bool(terminals, lambda t: t.name in used_terms or t.name in self.ignore or t.name in terminals_to_keep)
            if unused:
                logger.debug("Unused terminals: %s", [t.name for t in unused])

        return terminals, compiled_rules, self.ignore


PackageResource = namedtuple('PackageResource', 'pkg_name path')


class FromPackageLoader:
    """
    Provides a simple way of creating custom import loaders that load from packages via ``pkgutil.get_data`` instead of using `open`.
    This allows them to be compatible even from within zip files.

    Relative imports are handled, so you can just freely use them.

    pkg_name: The name of the package. You can probably provide `__name__` most of the time
    search_paths: All the path that will be search on absolute imports.
    """

    pkg_name: str
    search_paths: Sequence[str]

    def __init__(self, pkg_name: str, search_paths: Sequence[str]=("", )) -> None:
        self.pkg_name = pkg_name
        self.search_paths = search_paths

    def __repr__(self):
        return "%s(%r, %r)" % (type(self).__name__, self.pkg_name, self.search_paths)

    def __call__(self, base_path: Union[None, str, PackageResource], grammar_path: str) -> Tuple[PackageResource, str]:
        if base_path is None:
            to_try = self.search_paths
        else:
            # Check whether or not the importing grammar was loaded by this module.
            if not isinstance(base_path, PackageResource) or base_path.pkg_name != self.pkg_name:
                # Technically false, but FileNotFound doesn't exist in python2.7, and this message should never reach the end user anyway
                raise IOError()
            to_try = [base_path.path]

        err = None
        for path in to_try:
            full_path = os.path.join(path, grammar_path)
            try:
                text: Optional[bytes] = pkgutil.get_data(self.pkg_name, full_path)
            except IOError as e:
                err = e
                continue
            else:
                return PackageResource(self.pkg_name, full_path), (text.decode() if text else '')

        raise IOError('Cannot find grammar in given paths') from err


stdlib_loader = FromPackageLoader('lark', IMPORT_PATHS)



def resolve_term_references(term_dict):
    # TODO Solve with transitive closure (maybe)

    while True:
        changed = False
        for name, token_tree in term_dict.items():
            if token_tree is None:  # Terminal added through %declare
                continue
            for exp in token_tree.find_data('value'):
                item ,= exp.children
                if isinstance(item, NonTerminal):
                    raise GrammarError("Rules aren't allowed inside terminals (%s in %s)" % (item, name))
                elif isinstance(item, Terminal):
                    try:
                        term_value = term_dict[item.name]
                    except KeyError:
                        raise GrammarError("Terminal used but not defined: %s" % item.name)
                    assert term_value is not None
                    exp.children[0] = term_value
                    changed = True
                else:
                    assert isinstance(item, Tree)
        if not changed:
            break

    for name, term in term_dict.items():
        if term:    # Not just declared
            for child in term.children:
                ids = [id(x) for x in child.iter_subtrees()]
                if id(term) in ids:
                    raise GrammarError("Recursion in terminal '%s' (recursion is only allowed in rules, not terminals)" % name)



def symbol_from_strcase(s):
    assert isinstance(s, str)
    return Terminal(s, filter_out=s.startswith('_')) if s.isupper() else NonTerminal(s)

@inline_args
class PrepareGrammar(Transformer_InPlace):
    def terminal(self, name):
        return Terminal(str(name), filter_out=name.startswith('_'))

    def nonterminal(self, name):
        return NonTerminal(name.value)


def _find_used_symbols(tree):
    assert tree.data == 'expansions'
    return {t.name for x in tree.find_data('expansion')
            for t in x.scan_values(lambda t: isinstance(t, Symbol))}


def _get_parser():
    try:
        return _get_parser.cache
    except AttributeError:
        terminals = [TerminalDef(name, PatternRE(value)) for name, value in TERMINALS.items()]

        rules = [(name.lstrip('?'), x, RuleOptions(expand1=name.startswith('?')))
                for name, x in RULES.items()]
        rules = [Rule(NonTerminal(r), [symbol_from_strcase(s) for s in x.split()], i, None, o)
                 for r, xs, o in rules for i, x in enumerate(xs)]

        callback = ParseTreeBuilder(rules, ST).create_callback()
        import re
        lexer_conf = LexerConf(terminals, re, ['WS', 'COMMENT', 'BACKSLASH'])
        parser_conf = ParserConf(rules, callback, ['start'])
        lexer_conf.lexer_type = 'basic'
        parser_conf.parser_type = 'lalr'
        _get_parser.cache = ParsingFrontend(lexer_conf, parser_conf, None)
        return _get_parser.cache

GRAMMAR_ERRORS = [
        ('Incorrect type of value', ['a: 1\n']),
        ('Unclosed parenthesis', ['a: (\n']),
        ('Unmatched closing parenthesis', ['a: )\n', 'a: [)\n', 'a: (]\n']),
        ('Expecting rule or terminal definition (missing colon)', ['a\n', 'A\n', 'a->\n', 'A->\n', 'a A\n']),
        ('Illegal name for rules or terminals', ['Aa:\n']),
        ('Alias expects lowercase name', ['a: -> "a"\n']),
        ('Unexpected colon', ['a::\n', 'a: b:\n', 'a: B:\n', 'a: "a":\n']),
        ('Misplaced operator', ['a: b??', 'a: b(?)', 'a:+\n', 'a:?\n', 'a:*\n', 'a:|*\n']),
        ('Expecting option ("|") or a new rule or terminal definition', ['a:a\n()\n']),
        ('Terminal names cannot contain dots', ['A.B\n']),
        ('Expecting rule or terminal definition', ['"a"\n']),
        ('%import expects a name', ['%import "a"\n']),
        ('%ignore expects a value', ['%ignore %import\n']),
    ]

def _translate_parser_exception(parse, e):
        error = e.match_examples(parse, GRAMMAR_ERRORS, use_accepts=True)
        if error:
            return error
        elif 'STRING' in e.expected:
            return "Expecting a value"

def _parse_grammar(text, name, start='start'):
    try:
        tree = _get_parser().parse(text + '\n', start)
    except UnexpectedCharacters as e:
        context = e.get_context(text)
        raise GrammarError("Unexpected input at line %d column %d in %s: \n\n%s" %
                           (e.line, e.column, name, context))
    except UnexpectedToken as e:
        context = e.get_context(text)
        error = _translate_parser_exception(_get_parser().parse, e)
        if error:
            raise GrammarError("%s, at line %s column %s\n\n%s" % (error, e.line, e.column, context))
        raise

    return PrepareGrammar().transform(tree)

def _error_repr(error):
    if isinstance(error, UnexpectedToken):
        error2 = _translate_parser_exception(_get_parser().parse, error)
        if error2:
            return error2
        expected = ', '.join(error.accepts or error.expected)
        return "Unexpected token %r. Expected one of: {%s}" % (str(error.token), expected)
    else:
        return str(error)

def _search_interactive_parser(interactive_parser, predicate):
    def expand(node):
        path, p = node
        for choice in p.choices():
            t = Token(choice, '')
            try:
                new_p = p.feed_token(t)
            except ParseError:    # Illegal
                pass
            else:
                yield path + (choice,), new_p

    for path, p in bfs_all_unique([((), interactive_parser)], expand):
        if predicate(p):
            return path, p

def find_grammar_errors(text: str, start: str='start') -> List[Tuple[UnexpectedInput, str]]:
    errors = []
    def on_error(e):
        errors.append((e, _error_repr(e)))

        # recover to a new line
        token_path, _ = _search_interactive_parser(e.interactive_parser.as_immutable(), lambda p: '_NL' in p.choices())
        for token_type in token_path:
            e.interactive_parser.feed_token(Token(token_type, ''))
        e.interactive_parser.feed_token(Token('_NL', '\n'))
        return True

    _tree = _get_parser().parse(text + '\n', start, on_error=on_error)

    errors_by_line = classify(errors, lambda e: e[0].line)
    errors = [el[0] for el in errors_by_line.values()]      # already sorted

    for e in errors:
        e[0].interactive_parser = None
    return errors


def _get_mangle(prefix, aliases, base_mangle=None):
    def mangle(s):
        if s in aliases:
            s = aliases[s]
        else:
            if s[0] == '_':
                s = '_%s__%s' % (prefix, s[1:])
            else:
                s = '%s__%s' % (prefix, s)
        if base_mangle is not None:
            s = base_mangle(s)
        return s
    return mangle

def _mangle_definition_tree(exp, mangle):
    if mangle is None:
        return exp
    exp = deepcopy(exp) # TODO: is this needed?
    for t in exp.iter_subtrees():
        for i, c in enumerate(t.children):
            if isinstance(c, Symbol):
                t.children[i] = c.renamed(mangle)

    return exp

def _make_rule_tuple(modifiers_tree, name, params, priority_tree, expansions):
    if modifiers_tree.children:
        m ,= modifiers_tree.children
        expand1 = '?' in m
        if expand1 and name.startswith('_'):
            raise GrammarError("Inlined rules (_rule) cannot use the ?rule modifier.")
        keep_all_tokens = '!' in m
    else:
        keep_all_tokens = False
        expand1 = False

    if priority_tree.children:
        p ,= priority_tree.children
        priority = int(p)
    else:
        priority = None

    if params is not None:
        params = [t.value for t in params.children]  # For the grammar parser

    return name, params, expansions, RuleOptions(keep_all_tokens, expand1, priority=priority,
                                                 template_source=(name if params else None))


class Definition:
    def __init__(self, is_term, tree, params=(), options=None):
        self.is_term = is_term
        self.tree = tree
        self.params = tuple(params)
        self.options = options

class GrammarBuilder:

    global_keep_all_tokens: bool
    import_paths: List[Union[str, Callable]]
    used_files: Dict[str, str]

    _definitions: Dict[str, Definition]
    _ignore_names: List[str]

    def __init__(self, global_keep_all_tokens: bool=False, import_paths: Optional[List[Union[str, Callable]]]=None, used_files: Optional[Dict[str, str]]=None) -> None:
        self.global_keep_all_tokens = global_keep_all_tokens
        self.import_paths = import_paths or []
        self.used_files = used_files or {}

        self._definitions: Dict[str, Definition] = {}
        self._ignore_names: List[str] = []

    def _grammar_error(self, is_term, msg, *names):
        args = {}
        for i, name in enumerate(names, start=1):
            postfix = '' if i == 1 else str(i)
            args['name' + postfix] = name
            args['type' + postfix] = lowercase_type = ("rule", "terminal")[is_term]
            args['Type' + postfix] = lowercase_type.title()
        raise GrammarError(msg.format(**args))

    def _check_options(self, is_term, options):
        if is_term:
            if options is None:
                options = 1
            elif not isinstance(options, int):
                raise GrammarError("Terminal require a single int as 'options' (e.g. priority), got %s" % (type(options),))
        else:
            if options is None:
                options = RuleOptions()
            elif not isinstance(options, RuleOptions):
                raise GrammarError("Rules require a RuleOptions instance as 'options'")
            if self.global_keep_all_tokens:
                options.keep_all_tokens = True
        return options


    def _define(self, name, is_term, exp, params=(), options=None, *, override=False):
        if name in self._definitions:
            if not override:
                self._grammar_error(is_term, "{Type} '{name}' defined more than once", name)
        elif override:
            self._grammar_error(is_term, "Cannot override a nonexisting {type} {name}", name)

        if name.startswith('__'):
            self._grammar_error(is_term, 'Names starting with double-underscore are reserved (Error at {name})', name)

        self._definitions[name] = Definition(is_term, exp, params, self._check_options(is_term, options))

    def _extend(self, name, is_term, exp, params=(), options=None):
        if name not in self._definitions:
            self._grammar_error(is_term, "Can't extend {type} {name} as it wasn't defined before", name)

        d = self._definitions[name]

        if is_term != d.is_term:
            self._grammar_error(is_term, "Cannot extend {type} {name} - one is a terminal, while the other is not.", name)
        if tuple(params) != d.params:
            self._grammar_error(is_term, "Cannot extend {type} with different parameters: {name}", name)

        if d.tree is None:
            self._grammar_error(is_term, "Can't extend {type} {name} - it is abstract.", name)

        # TODO: think about what to do with 'options'
        base = d.tree

        assert isinstance(base, Tree) and base.data == 'expansions'
        base.children.insert(0, exp)

    def _ignore(self, exp_or_name):
        if isinstance(exp_or_name, str):
            self._ignore_names.append(exp_or_name)
        else:
            assert isinstance(exp_or_name, Tree)
            t = exp_or_name
            if t.data == 'expansions' and len(t.children) == 1:
                t2 ,= t.children
                if t2.data=='expansion' and len(t2.children) == 1:
                    item ,= t2.children
                    if item.data == 'value':
                        item ,= item.children
                        if isinstance(item, Terminal):
                            # Keep terminal name, no need to create a new definition
                            self._ignore_names.append(item.name)
                            return

            name = '__IGNORE_%d'% len(self._ignore_names)
            self._ignore_names.append(name)
            self._definitions[name] = Definition(True, t, options=TOKEN_DEFAULT_PRIORITY)

    def _unpack_import(self, stmt, grammar_name):
        if len(stmt.children) > 1:
            path_node, arg1 = stmt.children
        else:
            path_node, = stmt.children
            arg1 = None

        if isinstance(arg1, Tree):  # Multi import
            dotted_path = tuple(path_node.children)
            names = arg1.children
            aliases = dict(zip(names, names))  # Can't have aliased multi import, so all aliases will be the same as names
        else:  # Single import
            dotted_path = tuple(path_node.children[:-1])
            if not dotted_path:
                name ,= path_node.children
                raise GrammarError("Nothing was imported from grammar `%s`" % name)
            name = path_node.children[-1]  # Get name from dotted path
            aliases = {name.value: (arg1 or name).value}  # Aliases if exist

        if path_node.data == 'import_lib':  # Import from library
            base_path = None
        else:  # Relative import
            if grammar_name == '<string>':  # Import relative to script file path if grammar is coded in script
                try:
                    base_file = os.path.abspath(sys.modules['__main__'].__file__)
                except AttributeError:
                    base_file = None
            else:
                base_file = grammar_name  # Import relative to grammar file path if external grammar file
            if base_file:
                if isinstance(base_file, PackageResource):
                    base_path = PackageResource(base_file.pkg_name, os.path.split(base_file.path)[0])
                else:
                    base_path = os.path.split(base_file)[0]
            else:
                base_path = os.path.abspath(os.path.curdir)

        return dotted_path, base_path, aliases

    def _unpack_definition(self, tree, mangle):

        if tree.data == 'rule':
            name, params, exp, opts = _make_rule_tuple(*tree.children)
            is_term = False
        else:
            name = tree.children[0].value
            params = ()     # TODO terminal templates
            opts = int(tree.children[1]) if len(tree.children) == 3 else TOKEN_DEFAULT_PRIORITY # priority
            exp = tree.children[-1]
            is_term = True

        if mangle is not None:
            params = tuple(mangle(p) for p in params)
            name = mangle(name)

        exp = _mangle_definition_tree(exp, mangle)
        return name, is_term, exp, params, opts


    def load_grammar(self, grammar_text: str, grammar_name: str="<?>", mangle: Optional[Callable[[str], str]]=None) -> None:
        tree = _parse_grammar(grammar_text, grammar_name)

        imports: Dict[Tuple[str, ...], Tuple[Optional[str], Dict[str, str]]] = {}

        for stmt in tree.children:
            if stmt.data == 'import':
                dotted_path, base_path, aliases = self._unpack_import(stmt, grammar_name)
                try:
                    import_base_path, import_aliases = imports[dotted_path]
                    assert base_path == import_base_path, 'Inconsistent base_path for %s.' % '.'.join(dotted_path)
                    import_aliases.update(aliases)
                except KeyError:
                    imports[dotted_path] = base_path, aliases

        for dotted_path, (base_path, aliases) in imports.items():
            self.do_import(dotted_path, base_path, aliases, mangle)

        for stmt in tree.children:
            if stmt.data in ('term', 'rule'):
                self._define(*self._unpack_definition(stmt, mangle))
            elif stmt.data == 'override':
                r ,= stmt.children
                self._define(*self._unpack_definition(r, mangle), override=True)
            elif stmt.data == 'extend':
                r ,= stmt.children
                self._extend(*self._unpack_definition(r, mangle))
            elif stmt.data == 'ignore':
                # if mangle is not None, we shouldn't apply ignore, since we aren't in a toplevel grammar
                if mangle is None:
                    self._ignore(*stmt.children)
            elif stmt.data == 'declare':
                for symbol in stmt.children:
                    assert isinstance(symbol, Symbol), symbol
                    is_term = isinstance(symbol, Terminal)
                    if mangle is None:
                        name = symbol.name
                    else:
                        name = mangle(symbol.name)
                    self._define(name, is_term, None)
            elif stmt.data == 'import':
                pass
            else:
                assert False, stmt


        term_defs = { name: d.tree
            for name, d in self._definitions.items()
            if d.is_term
        }
        resolve_term_references(term_defs)


    def _remove_unused(self, used):
        def rule_dependencies(symbol):
            try:
                d = self._definitions[symbol]
            except KeyError:
                return []
            if d.is_term:
                return []
            return _find_used_symbols(d.tree) - set(d.params)

        _used = set(bfs(used, rule_dependencies))
        self._definitions = {k: v for k, v in self._definitions.items() if k in _used}


    def do_import(self, dotted_path: Tuple[str, ...], base_path: Optional[str], aliases: Dict[str, str], base_mangle: Optional[Callable[[str], str]]=None) -> None:
        assert dotted_path
        mangle = _get_mangle('__'.join(dotted_path), aliases, base_mangle)
        grammar_path = os.path.join(*dotted_path) + EXT
        to_try = self.import_paths + ([base_path] if base_path is not None else []) + [stdlib_loader]
        for source in to_try:
            try:
                if callable(source):
                    joined_path, text = source(base_path, grammar_path)
                else:
                    joined_path = os.path.join(source, grammar_path)
                    with open(joined_path, encoding='utf8') as f:
                        text = f.read()
            except IOError:
                continue
            else:
                h = sha256_digest(text)
                if self.used_files.get(joined_path, h) != h:
                    raise RuntimeError("Grammar file was changed during importing")
                self.used_files[joined_path] = h

                gb = GrammarBuilder(self.global_keep_all_tokens, self.import_paths, self.used_files)
                gb.load_grammar(text, joined_path, mangle)
                gb._remove_unused(map(mangle, aliases))
                for name in gb._definitions:
                    if name in self._definitions:
                        raise GrammarError("Cannot import '%s' from '%s': Symbol already defined." % (name, grammar_path))

                self._definitions.update(**gb._definitions)
                break
        else:
            # Search failed. Make Python throw a nice error.
            open(grammar_path, encoding='utf8')
            assert False, "Couldn't import grammar %s, but a corresponding file was found at a place where lark doesn't search for it" % (dotted_path,)


    def validate(self) -> None:
        for name, d in self._definitions.items():
            params = d.params
            exp = d.tree

            for i, p in enumerate(params):
                if p in self._definitions:
                    raise GrammarError("Template Parameter conflicts with rule %s (in template %s)" % (p, name))
                if p in params[:i]:
                    raise GrammarError("Duplicate Template Parameter %s (in template %s)" % (p, name))

            if exp is None: # Remaining checks don't apply to abstract rules/terminals (created with %declare)
                continue

            for temp in exp.find_data('template_usage'):
                sym = temp.children[0].name
                args = temp.children[1:]
                if sym not in params:
                    if sym not in self._definitions:
                        self._grammar_error(d.is_term, "Template '%s' used but not defined (in {type} {name})" % sym, name)
                    if len(args) != len(self._definitions[sym].params):
                        expected, actual = len(self._definitions[sym].params), len(args)
                        self._grammar_error(d.is_term, "Wrong number of template arguments used for {name} "
                                            "(expected %s, got %s) (in {type2} {name2})" % (expected, actual), sym, name)

            for sym in _find_used_symbols(exp):
                if sym not in self._definitions and sym not in params:
                    self._grammar_error(d.is_term, "{Type} '{name}' used but not defined (in {type2} {name2})", sym, name)

        if not set(self._definitions).issuperset(self._ignore_names):
            raise GrammarError("Terminals %s were marked to ignore but were not defined!" % (set(self._ignore_names) - set(self._definitions)))

    def build(self) -> Grammar:
        self.validate()
        rule_defs = []
        term_defs = []
        for name, d in self._definitions.items():
            (params, exp, options) = d.params, d.tree, d.options
            if d.is_term:
                assert len(params) == 0
                term_defs.append((name, (exp, options)))
            else:
                rule_defs.append((name, params, exp, options))
        # resolve_term_references(term_defs)
        return Grammar(rule_defs, term_defs, self._ignore_names)


def verify_used_files(file_hashes):
    for path, old in file_hashes.items():
        text = None
        if isinstance(path, str) and os.path.exists(path):
            with open(path, encoding='utf8') as f:
                text = f.read()
        elif isinstance(path, PackageResource):
            with suppress(IOError):
                text = pkgutil.get_data(*path).decode('utf-8')
        if text is None: # We don't know how to load the path. ignore it.
            continue

        current = sha256_digest(text)
        if old != current:
            logger.info("File %r changed, rebuilding Parser" % path)
            return False
    return True

def list_grammar_imports(grammar, import_paths=[]):
    "Returns a list of paths to the lark grammars imported by the given grammar (recursively)"
    builder = GrammarBuilder(False, import_paths)
    builder.load_grammar(grammar, '<string>')
    return list(builder.used_files.keys())

def load_grammar(grammar, source, import_paths, global_keep_all_tokens):
    builder = GrammarBuilder(global_keep_all_tokens, import_paths)
    builder.load_grammar(grammar, source)
    return builder.build(), builder.used_files


def sha256_digest(s: str) -> str:
    """Get the sha256 digest of a string

    Supports the `usedforsecurity` argument for Python 3.9+ to allow running on
    a FIPS-enabled system.
    """
    if sys.version_info >= (3, 9):
        return hashlib.sha256(s.encode('utf8'), usedforsecurity=False).hexdigest()
    else:
        return hashlib.sha256(s.encode('utf8')).hexdigest()
