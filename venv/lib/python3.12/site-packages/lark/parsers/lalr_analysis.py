"""This module builds a LALR(1) transition-table for lalr_parser.py

For now, shift/reduce conflicts are automatically resolved as shifts.
"""

# Author: Erez Shinan (2017)
# Email : erezshin@gmail.com

from typing import Dict, Set, Iterator, Tuple, List, TypeVar, Generic
from collections import defaultdict

from ..utils import classify, classify_bool, bfs, fzset, Enumerator, logger
from ..exceptions import GrammarError

from .grammar_analysis import GrammarAnalyzer, Terminal, LR0ItemSet, RulePtr, State
from ..grammar import Rule, Symbol
from ..common import ParserConf

###{standalone

class Action:
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return self.name
    def __repr__(self):
        return str(self)

Shift = Action('Shift')
Reduce = Action('Reduce')

StateT = TypeVar("StateT")

class ParseTableBase(Generic[StateT]):
    states: Dict[StateT, Dict[str, Tuple]]
    start_states: Dict[str, StateT]
    end_states: Dict[str, StateT]

    def __init__(self, states, start_states, end_states):
        self.states = states
        self.start_states = start_states
        self.end_states = end_states

    def serialize(self, memo):
        tokens = Enumerator()

        states = {
            state: {tokens.get(token): ((1, arg.serialize(memo)) if action is Reduce else (0, arg))
                    for token, (action, arg) in actions.items()}
            for state, actions in self.states.items()
        }

        return {
            'tokens': tokens.reversed(),
            'states': states,
            'start_states': self.start_states,
            'end_states': self.end_states,
        }

    @classmethod
    def deserialize(cls, data, memo):
        tokens = data['tokens']
        states = {
            state: {tokens[token]: ((Reduce, Rule.deserialize(arg, memo)) if action==1 else (Shift, arg))
                    for token, (action, arg) in actions.items()}
            for state, actions in data['states'].items()
        }
        return cls(states, data['start_states'], data['end_states'])

class ParseTable(ParseTableBase['State']):
    """Parse-table whose key is State, i.e. set[RulePtr]

    Slower than IntParseTable, but useful for debugging
    """
    pass


class IntParseTable(ParseTableBase[int]):
    """Parse-table whose key is int. Best for performance."""

    @classmethod
    def from_ParseTable(cls, parse_table: ParseTable):
        enum = list(parse_table.states)
        state_to_idx: Dict['State', int] = {s:i for i,s in enumerate(enum)}
        int_states = {}

        for s, la in parse_table.states.items():
            la = {k:(v[0], state_to_idx[v[1]]) if v[0] is Shift else v
                  for k,v in la.items()}
            int_states[ state_to_idx[s] ] = la


        start_states = {start:state_to_idx[s] for start, s in parse_table.start_states.items()}
        end_states = {start:state_to_idx[s] for start, s in parse_table.end_states.items()}
        return cls(int_states, start_states, end_states)

###}


# digraph and traverse, see The Theory and Practice of Compiler Writing

# computes F(x) = G(x) union (union { G(y) | x R y })
# X: nodes
# R: relation (function mapping node -> list of nodes that satisfy the relation)
# G: set valued function
def digraph(X, R, G):
    F = {}
    S = []
    N = dict.fromkeys(X, 0)
    for x in X:
        # this is always true for the first iteration, but N[x] may be updated in traverse below
        if N[x] == 0:
            traverse(x, S, N, X, R, G, F)
    return F

# x: single node
# S: stack
# N: weights
# X: nodes
# R: relation (see above)
# G: set valued function
# F: set valued function we are computing (map of input -> output)
def traverse(x, S, N, X, R, G, F):
    S.append(x)
    d = len(S)
    N[x] = d
    F[x] = G[x]
    for y in R[x]:
        if N[y] == 0:
            traverse(y, S, N, X, R, G, F)
        n_x = N[x]
        assert(n_x > 0)
        n_y = N[y]
        assert(n_y != 0)
        if (n_y > 0) and (n_y < n_x):
            N[x] = n_y
        F[x].update(F[y])
    if N[x] == d:
        f_x = F[x]
        while True:
            z = S.pop()
            N[z] = -1
            F[z] = f_x
            if z == x:
                break


class LALR_Analyzer(GrammarAnalyzer):
    lr0_itemsets: Set[LR0ItemSet]
    nonterminal_transitions: List[Tuple[LR0ItemSet, Symbol]]
    lookback: Dict[Tuple[LR0ItemSet, Symbol], Set[Tuple[LR0ItemSet, Rule]]]
    includes: Dict[Tuple[LR0ItemSet, Symbol], Set[Tuple[LR0ItemSet, Symbol]]]
    reads: Dict[Tuple[LR0ItemSet, Symbol], Set[Tuple[LR0ItemSet, Symbol]]]
    directly_reads: Dict[Tuple[LR0ItemSet, Symbol], Set[Symbol]]


    def __init__(self, parser_conf: ParserConf, debug: bool=False, strict: bool=False):
        GrammarAnalyzer.__init__(self, parser_conf, debug, strict)
        self.nonterminal_transitions = []
        self.directly_reads = defaultdict(set)
        self.reads = defaultdict(set)
        self.includes = defaultdict(set)
        self.lookback = defaultdict(set)


    def compute_lr0_states(self) -> None:
        self.lr0_itemsets = set()
        # map of kernels to LR0ItemSets
        cache: Dict['State', LR0ItemSet] = {}

        def step(state: LR0ItemSet) -> Iterator[LR0ItemSet]:
            _, unsat = classify_bool(state.closure, lambda rp: rp.is_satisfied)

            d = classify(unsat, lambda rp: rp.next)
            for sym, rps in d.items():
                kernel = fzset({rp.advance(sym) for rp in rps})
                new_state = cache.get(kernel, None)
                if new_state is None:
                    closure = set(kernel)
                    for rp in kernel:
                        if not rp.is_satisfied and not rp.next.is_term:
                            closure |= self.expand_rule(rp.next, self.lr0_rules_by_origin)
                    new_state = LR0ItemSet(kernel, closure)
                    cache[kernel] = new_state

                state.transitions[sym] = new_state
                yield new_state

            self.lr0_itemsets.add(state)

        for _ in bfs(self.lr0_start_states.values(), step):
            pass

    def compute_reads_relations(self):
        # handle start state
        for root in self.lr0_start_states.values():
            assert(len(root.kernel) == 1)
            for rp in root.kernel:
                assert(rp.index == 0)
                self.directly_reads[(root, rp.next)] = set([ Terminal('$END') ])

        for state in self.lr0_itemsets:
            seen = set()
            for rp in state.closure:
                if rp.is_satisfied:
                    continue
                s = rp.next
                # if s is a not a nonterminal
                if s not in self.lr0_rules_by_origin:
                    continue
                if s in seen:
                    continue
                seen.add(s)
                nt = (state, s)
                self.nonterminal_transitions.append(nt)
                dr = self.directly_reads[nt]
                r = self.reads[nt]
                next_state = state.transitions[s]
                for rp2 in next_state.closure:
                    if rp2.is_satisfied:
                        continue
                    s2 = rp2.next
                    # if s2 is a terminal
                    if s2 not in self.lr0_rules_by_origin:
                        dr.add(s2)
                    if s2 in self.NULLABLE:
                        r.add((next_state, s2))

    def compute_includes_lookback(self):
        for nt in self.nonterminal_transitions:
            state, nonterminal = nt
            includes = []
            lookback = self.lookback[nt]
            for rp in state.closure:
                if rp.rule.origin != nonterminal:
                    continue
                # traverse the states for rp(.rule)
                state2 = state
                for i in range(rp.index, len(rp.rule.expansion)):
                    s = rp.rule.expansion[i]
                    nt2 = (state2, s)
                    state2 = state2.transitions[s]
                    if nt2 not in self.reads:
                        continue
                    for j in range(i + 1, len(rp.rule.expansion)):
                        if rp.rule.expansion[j] not in self.NULLABLE:
                            break
                    else:
                        includes.append(nt2)
                # state2 is at the final state for rp.rule
                if rp.index == 0:
                    for rp2 in state2.closure:
                        if (rp2.rule == rp.rule) and rp2.is_satisfied:
                            lookback.add((state2, rp2.rule))
            for nt2 in includes:
                self.includes[nt2].add(nt)

    def compute_lookaheads(self):
        read_sets = digraph(self.nonterminal_transitions, self.reads, self.directly_reads)
        follow_sets = digraph(self.nonterminal_transitions, self.includes, read_sets)

        for nt, lookbacks in self.lookback.items():
            for state, rule in lookbacks:
                for s in follow_sets[nt]:
                    state.lookaheads[s].add(rule)

    def compute_lalr1_states(self) -> None:
        m: Dict[LR0ItemSet, Dict[str, Tuple]] = {}
        reduce_reduce = []
        for itemset in self.lr0_itemsets:
            actions: Dict[Symbol, Tuple] = {la: (Shift, next_state.closure)
                                                      for la, next_state in itemset.transitions.items()}
            for la, rules in itemset.lookaheads.items():
                if len(rules) > 1:
                    # Try to resolve conflict based on priority
                    p = [(r.options.priority or 0, r) for r in rules]
                    p.sort(key=lambda r: r[0], reverse=True)
                    best, second_best = p[:2]
                    if best[0] > second_best[0]:
                        rules = {best[1]}
                    else:
                        reduce_reduce.append((itemset, la, rules))
                        continue

                rule ,= rules
                if la in actions:
                    if self.strict:
                        raise GrammarError(f"Shift/Reduce conflict for terminal {la.name}. [strict-mode]\n ")
                    elif self.debug:
                        logger.warning('Shift/Reduce conflict for terminal %s: (resolving as shift)', la.name)
                        logger.warning(' * %s', rule)
                    else:
                        logger.debug('Shift/Reduce conflict for terminal %s: (resolving as shift)', la.name)
                        logger.debug(' * %s', rule)
                else:
                    actions[la] = (Reduce, rule)
            m[itemset] = { k.name: v for k, v in actions.items() }

        if reduce_reduce:
            msgs = []
            for itemset, la, rules in reduce_reduce:
                msg = 'Reduce/Reduce collision in %s between the following rules: %s' % (la, ''.join([ '\n\t- ' + str(r) for r in rules ]))
                if self.debug:
                    msg += '\n    collision occurred in state: {%s\n    }' % ''.join(['\n\t' + str(x) for x in itemset.closure])
                msgs.append(msg)
            raise GrammarError('\n\n'.join(msgs))

        states = { k.closure: v for k, v in m.items() }

        # compute end states
        end_states: Dict[str, 'State'] = {}
        for state in states:
            for rp in state:
                for start in self.lr0_start_states:
                    if rp.rule.origin.name == ('$root_' + start) and rp.is_satisfied:
                        assert start not in end_states
                        end_states[start] = state

        start_states = { start: state.closure for start, state in self.lr0_start_states.items() }
        _parse_table = ParseTable(states, start_states, end_states)

        if self.debug:
            self.parse_table = _parse_table
        else:
            self.parse_table = IntParseTable.from_ParseTable(_parse_table)

    def compute_lalr(self):
        self.compute_lr0_states()
        self.compute_reads_relations()
        self.compute_includes_lookback()
        self.compute_lookaheads()
        self.compute_lalr1_states()
