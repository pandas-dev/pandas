"""
Contains tests and a prototype implementation for the fanout algorithm in
the LLVM refprune pass.
"""

try:
    from graphviz import Digraph
except ImportError:
    pass
from collections import defaultdict

# The entry block. It's always the same.
ENTRY = "A"


# The following caseNN() functions returns a 3-tuple of
# (nodes, edges, expected).
# `nodes` maps BB nodes to incref/decref inside the block.
# `edges` maps BB nodes to their successor BB.
# `expected` maps BB-node with incref to a set of BB-nodes with the decrefs, or
#            the value can be None, indicating invalid prune.

def case1():
    edges = {
        "A": ["B"],
        "B": ["C", "D"],
        "C": [],
        "D": ["E", "F"],
        "E": ["G"],
        "F": [],
        "G": ["H", "I"],
        "I": ["G", "F"],
        "H": ["J", "K"],
        "J": ["L", "M"],
        "K": [],
        "L": ["Z"],
        "M": ["Z", "O", "P"],
        "O": ["Z"],
        "P": ["Z"],
        "Z": [],
    }
    nodes = defaultdict(list)
    nodes["D"] = ["incref"]
    nodes["H"] = ["decref"]
    nodes["F"] = ["decref", "decref"]
    expected = {"D": {"H", "F"}}
    return nodes, edges, expected


def case2():
    edges = {
        "A": ["B", "C"],
        "B": ["C"],
        "C": [],
    }
    nodes = defaultdict(list)
    nodes["A"] = ["incref"]
    nodes["B"] = ["decref"]
    nodes["C"] = ["decref"]
    expected = {"A": None}
    return nodes, edges, expected


def case3():
    nodes, edges, _ = case1()
    # adds an invalid edge
    edges["H"].append("F")
    expected = {"D": None}
    return nodes, edges, expected


def case4():
    nodes, edges, _ = case1()
    # adds an invalid edge
    edges["H"].append("E")
    expected = {"D": None}
    return nodes, edges, expected


def case5():
    nodes, edges, _ = case1()
    # adds backedge to go before incref
    edges["B"].append("I")
    expected = {"D": None}
    return nodes, edges, expected


def case6():
    nodes, edges, _ = case1()
    # adds backedge to go before incref
    edges["I"].append("B")
    expected = {"D": None}
    return nodes, edges, expected


def case7():
    nodes, edges, _ = case1()
    # adds forward jump outside
    edges["I"].append("M")
    expected = {"D": None}
    return nodes, edges, expected


def case8():
    edges = {
        "A": ["B", "C"],
        "B": ["C"],
        "C": [],
    }
    nodes = defaultdict(list)
    nodes["A"] = ["incref"]
    nodes["C"] = ["decref"]
    expected = {"A": {"C"}}
    return nodes, edges, expected


def case9():
    nodes, edges, _ = case8()
    # adds back edge
    edges["C"].append("B")
    expected = {"A": None}
    return nodes, edges, expected


def case10():
    nodes, edges, _ = case8()
    # adds back edge to A
    edges["C"].append("A")
    expected = {"A": {"C"}}
    return nodes, edges, expected


def case11():
    nodes, edges, _ = case8()
    edges["C"].append("D")
    edges["D"] = []
    expected = {"A": {"C"}}
    return nodes, edges, expected


def case12():
    nodes, edges, _ = case8()
    edges["C"].append("D")
    edges["D"] = ["A"]
    expected = {"A": {"C"}}
    return nodes, edges, expected


def case13():
    nodes, edges, _ = case8()
    edges["C"].append("D")
    edges["D"] = ["B"]
    expected = {"A": None}
    return nodes, edges, expected


def make_predecessor_map(edges):
    d = defaultdict(set)
    for src, outgoings in edges.items():
        for dst in outgoings:
            d[dst].add(src)
    return d


class FanoutAlgorithm:
    def __init__(self, nodes, edges, verbose=False):
        self.nodes = nodes
        self.edges = edges
        self.rev_edges = make_predecessor_map(edges)
        self.print = print if verbose else self._null_print

    def run(self):
        return self.find_fanout_in_function()

    def _null_print(self, *args, **kwargs):
        pass

    def find_fanout_in_function(self):
        got = {}
        for cur_node in self.edges:
            for incref in (x for x in self.nodes[cur_node] if x == "incref"):
                decref_blocks = self.find_fanout(cur_node)
                self.print(">>", cur_node, "===", decref_blocks)
                got[cur_node] = decref_blocks
        return got

    def find_fanout(self, head_node):
        decref_blocks = self.find_decref_candidates(head_node)
        self.print("candidates", decref_blocks)
        if not decref_blocks:
            return None
        if not self.verify_non_overlapping(
            head_node, decref_blocks, entry=ENTRY
        ):
            return None
        return set(decref_blocks)

    def verify_non_overlapping(self, head_node, decref_blocks, entry):
        self.print("verify_non_overlapping".center(80, "-"))
        # reverse walk for each decref_blocks
        # they should end at head_node
        todo = list(decref_blocks)
        while todo:
            cur_node = todo.pop()
            visited = set()

            workstack = [cur_node]
            del cur_node
            while workstack:
                cur_node = workstack.pop()
                self.print("cur_node", cur_node, "|", workstack)
                if cur_node in visited:
                    continue  # skip
                if cur_node == entry:
                    # Entry node
                    self.print(
                        "!! failed because we arrived at entry", cur_node
                    )
                    return False
                visited.add(cur_node)
                # check all predecessors
                self.print(
                    f"   {cur_node} preds {self.get_predecessors(cur_node)}"
                )
                for pred in self.get_predecessors(cur_node):
                    if pred in decref_blocks:
                        # reject because there's a predecessor in decref_blocks
                        self.print(
                            "!! reject because predecessor in decref_blocks"
                        )
                        return False
                    if pred != head_node:

                        workstack.append(pred)

        return True

    def get_successors(self, node):
        return tuple(self.edges[node])

    def get_predecessors(self, node):
        return tuple(self.rev_edges[node])

    def has_decref(self, node):
        return "decref" in self.nodes[node]

    def walk_child_for_decref(
        self, cur_node, path_stack, decref_blocks, depth=10
    ):
        indent = " " * len(path_stack)
        self.print(indent, "walk", path_stack, cur_node)
        if depth <= 0:
            return False  # missing
        if cur_node in path_stack:
            if cur_node == path_stack[0]:
                return False  # reject interior node backedge
            return True  # skip
        if self.has_decref(cur_node):
            decref_blocks.add(cur_node)
            self.print(indent, "found decref")
            return True

        depth -= 1
        path_stack += (cur_node,)
        found = False
        for child in self.get_successors(cur_node):
            if not self.walk_child_for_decref(
                child, path_stack, decref_blocks
            ):
                found = False
                break
            else:
                found = True

        self.print(indent, f"ret {found}")
        return found

    def find_decref_candidates(self, cur_node):
        # Forward pass
        self.print("find_decref_candidates".center(80, "-"))
        path_stack = (cur_node,)
        found = False
        decref_blocks = set()
        for child in self.get_successors(cur_node):
            if not self.walk_child_for_decref(
                child, path_stack, decref_blocks
            ):
                found = False
                break
            else:
                found = True
        if not found:
            return set()
        else:
            return decref_blocks


def check_once():
    nodes, edges, expected = case13()

    # Render graph
    G = Digraph()
    for node in edges:
        G.node(node, shape="rect", label=f"{node}\n" + r"\l".join(nodes[node]))
    for node, children in edges.items():
        for child in children:
            G.edge(node, child)

    G.view()

    algo = FanoutAlgorithm(nodes, edges, verbose=True)
    got = algo.run()
    assert expected == got


def check_all():
    for k, fn in list(globals().items()):
        if k.startswith("case"):
            print(f"{fn}".center(80, "-"))
            nodes, edges, expected = fn()
            algo = FanoutAlgorithm(nodes, edges)
            got = algo.run()
            assert expected == got
    print("ALL PASSED")


if __name__ == "__main__":
    # check_once()
    check_all()
