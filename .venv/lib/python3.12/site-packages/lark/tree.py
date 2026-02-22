import sys
from copy import deepcopy

from typing import List, Callable, Iterator, Union, Optional, Generic, TypeVar, TYPE_CHECKING

from .lexer import Token

if TYPE_CHECKING:
    from .lexer import TerminalDef
    try:
        import rich
    except ImportError:
        pass
    from typing import Literal

###{standalone

class Meta:

    empty: bool
    line: int
    column: int
    start_pos: int
    end_line: int
    end_column: int
    end_pos: int
    orig_expansion: 'List[TerminalDef]'
    match_tree: bool

    def __init__(self):
        self.empty = True


_Leaf_T = TypeVar("_Leaf_T")
Branch = Union[_Leaf_T, 'Tree[_Leaf_T]']


class Tree(Generic[_Leaf_T]):
    """The main tree class.

    Creates a new tree, and stores "data" and "children" in attributes of the same name.
    Trees can be hashed and compared.

    Parameters:
        data: The name of the rule or alias
        children: List of matched sub-rules and terminals
        meta: Line & Column numbers (if ``propagate_positions`` is enabled).
            meta attributes: (line, column, end_line, end_column, start_pos, end_pos,
                              container_line, container_column, container_end_line, container_end_column)
            container_* attributes consider all symbols, including those that have been inlined in the tree.
            For example, in the rule 'a: _A B _C', the regular attributes will mark the start and end of B,
            but the container_* attributes will also include _A and _C in the range. However, rules that
            contain 'a' will consider it in full, including _A and _C for all attributes.
    """

    data: str
    children: 'List[Branch[_Leaf_T]]'

    def __init__(self, data: str, children: 'List[Branch[_Leaf_T]]', meta: Optional[Meta]=None) -> None:
        self.data = data
        self.children = children
        self._meta = meta

    @property
    def meta(self) -> Meta:
        if self._meta is None:
            self._meta = Meta()
        return self._meta

    def __repr__(self):
        return 'Tree(%r, %r)' % (self.data, self.children)

    __match_args__ = ("data", "children")

    def _pretty_label(self):
        return self.data

    def _pretty(self, level, indent_str):
        yield f'{indent_str*level}{self._pretty_label()}'
        if len(self.children) == 1 and not isinstance(self.children[0], Tree):
            yield f'\t{self.children[0]}\n'
        else:
            yield '\n'
            for n in self.children:
                if isinstance(n, Tree):
                    yield from n._pretty(level+1, indent_str)
                else:
                    yield f'{indent_str*(level+1)}{n}\n'

    def pretty(self, indent_str: str='  ') -> str:
        """Returns an indented string representation of the tree.

        Great for debugging.
        """
        return ''.join(self._pretty(0, indent_str))

    def __rich__(self, parent:Optional['rich.tree.Tree']=None) -> 'rich.tree.Tree':
        """Returns a tree widget for the 'rich' library.

        Example:
            ::
                from rich import print
                from lark import Tree

                tree = Tree('root', ['node1', 'node2'])
                print(tree)
        """
        return self._rich(parent)

    def _rich(self, parent):
        if parent:
            tree = parent.add(f'[bold]{self.data}[/bold]')
        else:
            import rich.tree
            tree = rich.tree.Tree(self.data)

        for c in self.children:
            if isinstance(c, Tree):
                c._rich(tree)
            else:
                tree.add(f'[green]{c}[/green]')

        return tree

    def __eq__(self, other):
        try:
            return self.data == other.data and self.children == other.children
        except AttributeError:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self) -> int:
        return hash((self.data, tuple(self.children)))

    def iter_subtrees(self) -> 'Iterator[Tree[_Leaf_T]]':
        """Depth-first iteration.

        Iterates over all the subtrees, never returning to the same node twice (Lark's parse-tree is actually a DAG).
        """
        queue = [self]
        subtrees = dict()
        for subtree in queue:
            subtrees[id(subtree)] = subtree
            queue += [c for c in reversed(subtree.children)
                      if isinstance(c, Tree) and id(c) not in subtrees]

        del queue
        return reversed(list(subtrees.values()))

    def iter_subtrees_topdown(self):
        """Breadth-first iteration.

        Iterates over all the subtrees, return nodes in order like pretty() does.
        """
        stack = [self]
        stack_append = stack.append
        stack_pop = stack.pop
        while stack:
            node = stack_pop()
            if not isinstance(node, Tree):
                continue
            yield node
            for child in reversed(node.children):
                stack_append(child)

    def find_pred(self, pred: 'Callable[[Tree[_Leaf_T]], bool]') -> 'Iterator[Tree[_Leaf_T]]':
        """Returns all nodes of the tree that evaluate pred(node) as true."""
        return filter(pred, self.iter_subtrees())

    def find_data(self, data: str) -> 'Iterator[Tree[_Leaf_T]]':
        """Returns all nodes of the tree whose data equals the given data."""
        return self.find_pred(lambda t: t.data == data)

###}

    def find_token(self, token_type: str) -> Iterator[_Leaf_T]:
        """Returns all tokens whose type equals the given token_type.

        This is a recursive function that will find tokens in all the subtrees.

        Example:
            >>> term_tokens = tree.find_token('TERM')
        """
        return self.scan_values(lambda v: isinstance(v, Token) and v.type == token_type)

    def expand_kids_by_data(self, *data_values):
        """Expand (inline) children with any of the given data values. Returns True if anything changed"""
        changed = False
        for i in range(len(self.children)-1, -1, -1):
            child = self.children[i]
            if isinstance(child, Tree) and child.data in data_values:
                self.children[i:i+1] = child.children
                changed = True
        return changed


    def scan_values(self, pred: 'Callable[[Branch[_Leaf_T]], bool]') -> Iterator[_Leaf_T]:
        """Return all values in the tree that evaluate pred(value) as true.

        This can be used to find all the tokens in the tree.

        Example:
            >>> all_tokens = tree.scan_values(lambda v: isinstance(v, Token))
        """
        for c in self.children:
            if isinstance(c, Tree):
                for t in c.scan_values(pred):
                    yield t
            else:
                if pred(c):
                    yield c

    def __deepcopy__(self, memo):
        return type(self)(self.data, deepcopy(self.children, memo), meta=self._meta)

    def copy(self) -> 'Tree[_Leaf_T]':
        return type(self)(self.data, self.children)

    def set(self, data: str, children: 'List[Branch[_Leaf_T]]') -> None:
        self.data = data
        self.children = children


ParseTree = Tree['Token']


class SlottedTree(Tree):
    __slots__ = 'data', 'children', 'rule', '_meta'


def pydot__tree_to_png(tree: Tree, filename: str, rankdir: 'Literal["TB", "LR", "BT", "RL"]'="LR", **kwargs) -> None:
    graph = pydot__tree_to_graph(tree, rankdir, **kwargs)
    graph.write_png(filename)


def pydot__tree_to_dot(tree: Tree, filename, rankdir="LR", **kwargs):
    graph = pydot__tree_to_graph(tree, rankdir, **kwargs)
    graph.write(filename)


def pydot__tree_to_graph(tree: Tree, rankdir="LR", **kwargs):
    """Creates a colorful image that represents the tree (data+children, without meta)

    Possible values for `rankdir` are "TB", "LR", "BT", "RL", corresponding to
    directed graphs drawn from top to bottom, from left to right, from bottom to
    top, and from right to left, respectively.

    `kwargs` can be any graph attribute (e. g. `dpi=200`). For a list of
    possible attributes, see https://www.graphviz.org/doc/info/attrs.html.
    """

    import pydot  # type: ignore[import-not-found]
    graph = pydot.Dot(graph_type='digraph', rankdir=rankdir, **kwargs)

    i = [0]

    def new_leaf(leaf):
        node = pydot.Node(i[0], label=repr(leaf))
        i[0] += 1
        graph.add_node(node)
        return node

    def _to_pydot(subtree):
        color = hash(subtree.data) & 0xffffff
        color |= 0x808080

        subnodes = [_to_pydot(child) if isinstance(child, Tree) else new_leaf(child)
                    for child in subtree.children]
        node = pydot.Node(i[0], style="filled", fillcolor="#%x" % color, label=subtree.data)
        i[0] += 1
        graph.add_node(node)

        for subnode in subnodes:
            graph.add_edge(pydot.Edge(node, subnode))

        return node

    _to_pydot(tree)
    return graph
