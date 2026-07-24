"""This module implements an Earley parser.

The core Earley algorithm used here is based on Elizabeth Scott's implementation, here:
    https://www.sciencedirect.com/science/article/pii/S1571066108001497

That is probably the best reference for understanding the algorithm here.

The Earley parser outputs an SPPF-tree as per that document. The SPPF tree format
is explained here: https://lark-parser.readthedocs.io/en/latest/_static/sppf/sppf.html
"""

from typing import TYPE_CHECKING, Callable, Optional, List, Any
from collections import deque

from ..lexer import Token
from ..tree import Tree
from ..exceptions import UnexpectedEOF, UnexpectedToken
from ..utils import logger, OrderedSet, dedup_list
from .grammar_analysis import GrammarAnalyzer
from ..grammar import NonTerminal
from .earley_common import Item
from .earley_forest import ForestSumVisitor, SymbolNode, StableSymbolNode, TokenNode, ForestToParseTree

if TYPE_CHECKING:
    from ..common import LexerConf, ParserConf

class Parser:
    lexer_conf: 'LexerConf'
    parser_conf: 'ParserConf'
    debug: bool

    def __init__(self, lexer_conf: 'LexerConf', parser_conf: 'ParserConf', term_matcher: Callable,
                 resolve_ambiguity: bool=True, debug: bool=False,
                 tree_class: Optional[Callable[[str, List], Any]]=Tree, ordered_sets: bool=True):
        analysis = GrammarAnalyzer(parser_conf)
        self.lexer_conf = lexer_conf
        self.parser_conf = parser_conf
        self.resolve_ambiguity = resolve_ambiguity
        self.debug = debug
        self.Tree = tree_class
        self.Set = OrderedSet if ordered_sets else set
        self.SymbolNode = StableSymbolNode if ordered_sets else SymbolNode

        self.FIRST = analysis.FIRST
        self.NULLABLE = analysis.NULLABLE
        self.callbacks = parser_conf.callbacks
        # TODO add typing info
        self.predictions = {}   # type: ignore[var-annotated]

        ## These could be moved to the grammar analyzer. Pre-computing these is *much* faster than
        #  the slow 'isupper' in is_terminal.
        self.TERMINALS = { sym for r in parser_conf.rules for sym in r.expansion if sym.is_term }
        self.NON_TERMINALS = { sym for r in parser_conf.rules for sym in r.expansion if not sym.is_term }

        self.forest_sum_visitor = None
        for rule in parser_conf.rules:
            if rule.origin not in self.predictions:
                self.predictions[rule.origin] = [x.rule for x in analysis.expand_rule(rule.origin)]

            ## Detect if any rules/terminals have priorities set. If the user specified priority = None, then
            #  the priorities will be stripped from all rules/terminals before they reach us, allowing us to
            #  skip the extra tree walk. We'll also skip this if the user just didn't specify priorities
            #  on any rules/terminals.
            if self.forest_sum_visitor is None and rule.options.priority is not None:
                self.forest_sum_visitor = ForestSumVisitor

        # Check terminals for priorities
        # Ignore terminal priorities if the basic lexer is used
        if self.lexer_conf.lexer_type != 'basic' and self.forest_sum_visitor is None:
            for term in self.lexer_conf.terminals:
                if term.priority:
                    self.forest_sum_visitor = ForestSumVisitor
                    break

        self.term_matcher = term_matcher


    def predict_and_complete(self, i, to_scan, columns, transitives, node_cache):
        """The core Earley Predictor and Completer.

        At each stage of the input, we handling any completed items (things
        that matched on the last cycle) and use those to predict what should
        come next in the input stream. The completions and any predicted
        non-terminals are recursively processed until we reach a set of,
        which can be added to the scan list for the next scanner cycle."""
        # Held Completions (H in E.Scotts paper).
        held_completions = {}

        column = columns[i]
        # R (items) = Ei (column.items)
        items = deque(column)
        while items:
            item = items.pop()    # remove an element, A say, from R

            ### The Earley completer
            if item.is_complete:   ### (item.s == string)
                if item.node is None:
                    label = (item.s, item.start, i)
                    item.node = node_cache[label] if label in node_cache else node_cache.setdefault(label, self.SymbolNode(*label))
                    item.node.add_family(item.s, item.rule, item.start, None, None)

                # create_leo_transitives(item.rule.origin, item.start)

                ###R Joop Leo right recursion Completer
                if item.rule.origin in transitives[item.start]:
                    transitive = transitives[item.start][item.s]
                    if transitive.previous in transitives[transitive.column]:
                        root_transitive = transitives[transitive.column][transitive.previous]
                    else:
                        root_transitive = transitive

                    new_item = Item(transitive.rule, transitive.ptr, transitive.start)
                    label = (root_transitive.s, root_transitive.start, i)
                    new_item.node = node_cache[label] if label in node_cache else node_cache.setdefault(label, self.SymbolNode(*label))
                    new_item.node.add_path(root_transitive, item.node)
                    if new_item.expect in self.TERMINALS:
                        # Add (B :: aC.B, h, y) to Q
                        to_scan.add(new_item)
                    elif new_item not in column:
                        # Add (B :: aC.B, h, y) to Ei and R
                        column.add(new_item)
                        items.append(new_item)
                ###R Regular Earley completer
                else:
                    # Empty has 0 length. If we complete an empty symbol in a particular
                    # parse step, we need to be able to use that same empty symbol to complete
                    # any predictions that result, that themselves require empty. Avoids
                    # infinite recursion on empty symbols.
                    # held_completions is 'H' in E.Scott's paper.
                    is_empty_item = item.start == i
                    if is_empty_item:
                        held_completions[item.rule.origin] = item.node

                    originators = [originator for originator in columns[item.start] if originator.expect is not None and originator.expect == item.s]
                    for originator in originators:
                        new_item = originator.advance()
                        label = (new_item.s, originator.start, i)
                        new_item.node = node_cache[label] if label in node_cache else node_cache.setdefault(label, self.SymbolNode(*label))
                        new_item.node.add_family(new_item.s, new_item.rule, i, originator.node, item.node)
                        if new_item.expect in self.TERMINALS:
                            # Add (B :: aC.B, h, y) to Q
                            to_scan.add(new_item)
                        elif new_item not in column:
                            # Add (B :: aC.B, h, y) to Ei and R
                            column.add(new_item)
                            items.append(new_item)

            ### The Earley predictor
            elif item.expect in self.NON_TERMINALS: ### (item.s == lr0)
                new_items = []
                for rule in self.predictions[item.expect]:
                    new_item = Item(rule, 0, i)
                    new_items.append(new_item)

                # Process any held completions (H).
                if item.expect in held_completions:
                    new_item = item.advance()
                    label = (new_item.s, item.start, i)
                    new_item.node = node_cache[label] if label in node_cache else node_cache.setdefault(label, self.SymbolNode(*label))
                    new_item.node.add_family(new_item.s, new_item.rule, new_item.start, item.node, held_completions[item.expect])
                    new_items.append(new_item)

                for new_item in new_items:
                    if new_item.expect in self.TERMINALS:
                        to_scan.add(new_item)
                    elif new_item not in column:
                        column.add(new_item)
                        items.append(new_item)

    def _parse(self, lexer, columns, to_scan, start_symbol=None):

        def is_quasi_complete(item):
            if item.is_complete:
                return True

            quasi = item.advance()
            while not quasi.is_complete:
                if quasi.expect not in self.NULLABLE:
                    return False
                if quasi.rule.origin == start_symbol and quasi.expect == start_symbol:
                    return False
                quasi = quasi.advance()
            return True

        # def create_leo_transitives(origin, start):
        #   ...   # removed at commit 4c1cfb2faf24e8f8bff7112627a00b94d261b420

        def scan(i, token, to_scan):
            """The core Earley Scanner.

            This is a custom implementation of the scanner that uses the
            Lark lexer to match tokens. The scan list is built by the
            Earley predictor, based on the previously completed tokens.
            This ensures that at each phase of the parse we have a custom
            lexer context, allowing for more complex ambiguities."""
            next_to_scan = self.Set()
            next_set = self.Set()
            columns.append(next_set)
            transitives.append({})
            node_cache = {}

            for item in self.Set(to_scan):
                if match(item.expect, token):
                    new_item = item.advance()
                    label = (new_item.s, new_item.start, i + 1)
                    # 'terminals' may not contain token.type when using %declare
                    # Additionally, token is not always a Token
                    # For example, it can be a Tree when using TreeMatcher
                    term = terminals.get(token.type) if isinstance(token, Token) else None
                    # Set the priority of the token node to 0 so that the
                    # terminal priorities do not affect the Tree chosen by
                    # ForestSumVisitor after the basic lexer has already
                    # "used up" the terminal priorities
                    token_node = TokenNode(token, term, priority=0)
                    new_item.node = node_cache[label] if label in node_cache else node_cache.setdefault(label, self.SymbolNode(*label))
                    new_item.node.add_family(new_item.s, item.rule, new_item.start, item.node, token_node)

                    if new_item.expect in self.TERMINALS:
                        # add (B ::= Aai+1.B, h, y) to Q'
                        next_to_scan.add(new_item)
                    else:
                        # add (B ::= Aa+1.B, h, y) to Ei+1
                        next_set.add(new_item)

            if not next_set and not next_to_scan:
                expect = {i.expect.name for i in to_scan}
                raise UnexpectedToken(token, expect, considered_rules=set(to_scan), state=frozenset(i.s for i in to_scan))

            return next_to_scan, node_cache


        # Define parser functions
        match = self.term_matcher

        terminals = self.lexer_conf.terminals_by_name

        # Cache for nodes & tokens created in a particular parse step.
        transitives = [{}]

        ## The main Earley loop.
        # Run the Prediction/Completion cycle for any Items in the current Earley set.
        # Completions will be added to the SPPF tree, and predictions will be recursively
        # processed down to terminals/empty nodes to be added to the scanner for the next
        # step.
        expects = {i.expect for i in to_scan}
        i = 0
        node_cache = {}
        for token in lexer.lex(expects):
            self.predict_and_complete(i, to_scan, columns, transitives, node_cache)

            to_scan, node_cache = scan(i, token, to_scan)
            i += 1

            expects.clear()
            expects |= {i.expect for i in to_scan}

        self.predict_and_complete(i, to_scan, columns, transitives, node_cache)

        ## Column is now the final column in the parse.
        assert i == len(columns)-1
        return to_scan

    def parse(self, lexer, start):
        assert start, start
        start_symbol = NonTerminal(start)

        columns = [self.Set()]
        to_scan = self.Set()     # The scan buffer. 'Q' in E.Scott's paper.

        ## Predict for the start_symbol.
        # Add predicted items to the first Earley set (for the predictor) if they
        # result in a non-terminal, or the scanner if they result in a terminal.
        for rule in self.predictions[start_symbol]:
            item = Item(rule, 0, 0)
            if item.expect in self.TERMINALS:
                to_scan.add(item)
            else:
                columns[0].add(item)

        to_scan = self._parse(lexer, columns, to_scan, start_symbol)

        # If the parse was successful, the start
        # symbol should have been completed in the last step of the Earley cycle, and will be in
        # this column. Find the item for the start_symbol, which is the root of the SPPF tree.
        solutions = dedup_list(n.node for n in columns[-1] if n.is_complete and n.node is not None and n.s == start_symbol and n.start == 0)
        if not solutions:
            expected_terminals = [t.expect.name for t in to_scan]
            raise UnexpectedEOF(expected_terminals, state=frozenset(i.s for i in to_scan))
        if len(solutions) > 1:
            raise RuntimeError('Earley should not generate multiple start symbol items! Please report this bug.')
        solution ,= solutions

        if self.debug:
            from .earley_forest import ForestToPyDotVisitor
            try:
                debug_walker = ForestToPyDotVisitor()
            except ImportError:
                logger.warning("Cannot find dependency 'pydot', will not generate sppf debug image")
            else:
                debug_walker.visit(solution, "sppf.png")


        if self.Tree is not None:
            # Perform our SPPF -> AST conversion
            # Disable the ForestToParseTree cache when ambiguity='resolve'
            # to prevent a tree construction bug. See issue #1283
            use_cache = not self.resolve_ambiguity
            transformer = ForestToParseTree(self.Tree, self.callbacks, self.forest_sum_visitor and self.forest_sum_visitor(), self.resolve_ambiguity, use_cache)
            return transformer.transform(solution)

        # return the root of the SPPF
        return solution
