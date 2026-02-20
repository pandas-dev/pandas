"""This utility prints the contents of an HDF5 file as a tree.

Pass the flag -h to this for help on usage.

"""

import sys
import argparse
import warnings
from pathlib import Path
from collections import defaultdict, deque

import numpy as np

import tables as tb


def _get_parser():
    parser = argparse.ArgumentParser(
        description=(
            "`pttree` is designed to give a quick overview of the contents "
            "of a PyTables HDF5 file by printing a depth-indented list of "
            "nodes, similar to the output of the Unix `tree` function. "
            "It can also display the size, shape and compression states of "
            "individual nodes, as well as summary information for the whole "
            "file. "
            "For a more verbose output (including metadata), see `ptdump`. "
        ),
    )

    parser.add_argument(
        "-L",
        "--max-level",
        type=int,
        dest="max_depth",
        help="maximum branch depth of tree to display (-1 == no limit)",
    )
    parser.add_argument(
        "-S",
        "--sort-by",
        type=str,
        dest="sort_by",
        help='artificially order nodes, can be either "size", "name" or "none"',
    )
    parser.add_argument(
        "--print-size",
        action="store_true",
        dest="print_size",
        help="print size of each node/branch",
    )
    parser.add_argument(
        "--no-print-size",
        action="store_false",
        dest="print_size",
    )
    parser.add_argument(
        "--print-shape",
        action="store_true",
        dest="print_shape",
        help="print shape of each node",
    )
    parser.add_argument(
        "--no-print-shape",
        action="store_false",
        dest="print_shape",
    )
    parser.add_argument(
        "--print-compression",
        action="store_true",
        dest="print_compression",
        help="print compression library(level) for each compressed node",
    )
    parser.add_argument(
        "--no-print-compression",
        action="store_false",
        dest="print_compression",
    )
    parser.add_argument(
        "--print-percent",
        action="store_true",
        dest="print_percent",
        help="print size of each node as a %% of the total tree size on disk",
    )
    parser.add_argument(
        "--no-print-percent",
        action="store_false",
        dest="print_percent",
    )
    parser.add_argument(
        "--use-si-units",
        action="store_true",
        dest="use_si_units",
        help="report sizes in SI units (1 MB == 10^6 B)",
    )
    parser.add_argument(
        "--use-binary-units",
        action="store_false",
        dest="use_si_units",
        help="report sizes in binary units (1 MiB == 2^20 B)",
    )

    parser.add_argument(
        "src",
        metavar="filename[:nodepath]",
        help="path to the root of the tree structure",
    )

    parser.set_defaults(
        max_depth=1,
        sort_by="size",
        print_size=True,
        print_percent=True,
        print_shape=False,
        print_compression=False,
        use_si_units=False,
    )

    return parser


def main():
    """Implement the main CLI interface."""
    parser = _get_parser()
    args = parser.parse_args()

    # Catch the files passed as the last arguments
    src = args.__dict__.pop("src").rsplit(":", 1)
    if len(src) == 1:
        filename, nodename = src[0], "/"
    else:
        filename, nodename = src
        if nodename == "":
            # case where filename == "filename:" instead of "filename:/"
            nodename = "/"

    with tb.open_file(filename, "r") as f:
        tree_str = get_tree_str(f, nodename, **args.__dict__)
        print(tree_str)

    pass


def get_tree_str(
    f,
    where="/",
    max_depth=-1,
    print_class=True,
    print_size=True,
    print_percent=True,
    print_shape=False,
    print_compression=False,
    print_total=True,
    sort_by=None,
    use_si_units=False,
):
    """Return a string representing the tree structure, and the summary info."""
    root = f.get_node(where)
    root._g_check_open()
    start_depth = root._v_depth
    if max_depth < 0:
        max_depth = sys.maxsize

    b2h = bytes2human(use_si_units)

    # we will pass over each node in the tree twice

    # on the first pass we'll start at the root node and recurse down the
    # branches, finding all of the leaf nodes and calculating the total size
    # over all tables and arrays
    total_in_mem = 0
    total_on_disk = 0
    total_items = 0

    # defaultdicts for holding the cumulative branch sizes at each node
    in_mem = defaultdict(int)
    on_disk = defaultdict(int)
    leaf_count = defaultdict(int)

    # keep track of node addresses within the HDF5 file so that we don't count
    # nodes with multiple references (i.e. hardlinks) multiple times
    ref_count = defaultdict(int)
    ref_idx = defaultdict(int)
    hl_addresses = defaultdict(lambda: None)
    hl_targets = defaultdict(str)

    stack = deque(root)
    leaves = deque()

    while stack:

        node = stack.pop()

        if isinstance(node, tb.link.Link):
            # we treat links like leaves, except we don't dereference them to
            # get their sizes or addresses
            leaves.append(node)
            continue

        path = node._v_pathname
        addr, rc = node._get_obj_info()
        ref_count[addr] += 1
        ref_idx[path] = ref_count[addr]
        hl_addresses[path] = addr

        if isinstance(node, tb.UnImplemented):
            leaves.append(node)

        elif isinstance(node, tb.Leaf):

            # only count the size of a hardlinked leaf the first time it is
            # visited
            if ref_count[addr] == 1:

                try:
                    m = node.size_in_memory
                    d = node.size_on_disk

                    # size of this node
                    in_mem[path] += m
                    on_disk[path] += d
                    leaf_count[path] += 1

                    # total over all nodes
                    total_in_mem += m
                    total_on_disk += d
                    total_items += 1

                    # arbitrarily treat this node as the 'target' for all other
                    # hardlinks that point to the same address
                    hl_targets[addr] = path

                except NotImplementedError as e:
                    # size_on_disk is not implemented for VLArrays
                    warnings.warn(str(e))

            # push leaf nodes onto the stack for the next pass
            leaves.append(node)

        elif isinstance(node, tb.Group):

            # don't recurse down the same hardlinked branch multiple times!
            if ref_count[addr] == 1:
                stack.extend(list(node._v_children.values()))
                hl_targets[addr] = path

            # if we've already visited this group's address, treat it as a leaf
            # instead
            else:
                leaves.append(node)

    # on the second pass we start at each leaf and work upwards towards the
    # root node, computing the cumulative size of each branch at each node, and
    # instantiating a PrettyTree object for each node to create an ASCII
    # representation of the tree structure

    # this will store the PrettyTree objects for every node we're printing
    pretty = {}

    stack = leaves

    while stack:

        node = stack.pop()
        path = node._v_pathname

        parent = node._v_parent
        parent_path = parent._v_pathname

        # cumulative size at parent node
        in_mem[parent_path] += in_mem[path]
        on_disk[parent_path] += on_disk[path]
        leaf_count[parent_path] += leaf_count[path]

        depth = node._v_depth - start_depth

        # if we're deeper than the max recursion depth, we print nothing
        if not depth > max_depth:

            # create a PrettyTree representation of this node
            name = node._v_name
            if print_class:
                name += " (%s)" % node.__class__.__name__

            labels = []
            ratio = on_disk[path] / total_on_disk

            # if the address of this object has a ref_count > 1, it has
            # multiple hardlinks
            if ref_count[hl_addresses[path]] > 1:
                name += ", addr=%i, ref=%i/%i" % (
                    hl_addresses[path],
                    ref_idx[path],
                    ref_count[hl_addresses[path]],
                )

            if isinstance(node, tb.link.Link):
                labels.append("softlink --> %s" % node.target)

            elif ref_idx[path] > 1:
                labels.append(
                    "hardlink --> %s" % hl_targets[hl_addresses[path]]
                )

            elif isinstance(node, (tb.Array, tb.Table)):

                if print_size:
                    sizestr = "mem={}, disk={}".format(
                        b2h(in_mem[path]), b2h(on_disk[path])
                    )
                    if print_percent:
                        sizestr += f" [{ratio:5.1%}]"
                    labels.append(sizestr)

                if print_shape:
                    labels.append("shape=%s" % repr(node.shape))

                if print_compression:
                    lib = node.filters.complib
                    level = node.filters.complevel
                    if level:
                        compstr = "%s(%i)" % (lib, level)
                    else:
                        compstr = "None"
                    labels.append("compression=%s" % compstr)

            # if we're at our max recursion depth, we'll print summary
            # information for this branch
            elif depth == max_depth:
                itemstr = "... %i leaves" % leaf_count[path]
                if print_size:
                    itemstr += ", mem={}, disk={}".format(
                        b2h(in_mem[path]), b2h(on_disk[path])
                    )
                if print_percent:
                    itemstr += f" [{ratio:5.1%}]"
                labels.append(itemstr)

            # create a PrettyTree for this node, if one doesn't exist already
            if path not in pretty:
                pretty.update({path: PrettyTree()})
            pretty[path].name = name
            pretty[path].labels = labels
            if sort_by == "size":
                # descending size order
                pretty[path].sort_by = -ratio
            elif sort_by == "name":
                pretty[path].sort_by = node._v_name
            else:
                # natural order
                if path == "/":
                    # root is not in root._v_children
                    pretty[path].sort_by = 0
                else:
                    pretty[path].sort_by = list(
                        parent._v_children.values()
                    ).index(node)

            # exclude root node or we'll get infinite recursions (since '/' is
            # the parent of '/')
            if path != "/":

                # create a PrettyTree for the parent of this node, if one
                # doesn't exist already
                if parent_path not in pretty:
                    pretty.update({parent_path: PrettyTree()})

                # make this PrettyTree a child of the parent PrettyTree
                pretty[parent_path].add_child(pretty[path])

        if node is not root and parent not in stack:
            # we append to the 'bottom' of the stack, so that we exhaust all of
            # the nodes at this level before going up a level in the tree
            stack.appendleft(parent)

    out_str = "\n" + "-" * 60 + "\n" * 2
    out_str += str(pretty[root._v_pathname]) + "\n" * 2

    if print_total:
        avg_ratio = total_on_disk / total_in_mem
        fsize = Path(f.filename).stat().st_size

        out_str += "-" * 60 + "\n"
        out_str += "Total branch leaves:    %i\n" % total_items
        out_str += "Total branch size:      {} in memory, {} on disk\n".format(
            b2h(total_in_mem), b2h(total_on_disk)
        )
        out_str += "Mean compression ratio: %.2f\n" % avg_ratio
        out_str += "HDF5 file size:         %s\n" % b2h(fsize)
        out_str += "-" * 60 + "\n"

    return out_str


class PrettyTree:
    """Pretty ASCII representation of a recursive tree structure.

    Each node can have multiple labels, given as a list of strings.

    Example:
    --------
        A = PrettyTree('A', labels=['wow'])
        B = PrettyTree('B', labels=['such tree'])
        C = PrettyTree('C', children=[A, B])
        D = PrettyTree('D', labels=['so recursive'])
        root = PrettyTree('root', labels=['many nodes'], children=[C, D])
        print root

    Credit to Andrew Cooke's blog:
    <http://www.acooke.org/cute/ASCIIDispl0.html>

    """

    def __init__(self, name=None, children=None, labels=None, sort_by=None):

        # NB: do NOT assign default list/dict arguments in the function
        # declaration itself - these objects are shared between ALL instances
        # of PrettyTree, and by assigning to them it's easy to get into
        # infinite recursions, e.g. when 'self in self.children == True'
        if children is None:
            children = []
        if labels is None:
            labels = []

        self.name = name
        self.children = children
        self.labels = labels
        self.sort_by = sort_by

    def add_child(self, child):
        """Add a child to the tree."""
        # some basic checks to help to avoid infinite recursion
        assert child is not self
        assert self not in child.children
        if child not in self.children:
            self.children.append(child)

    def tree_lines(self):
        """Generate lines of teh string representation of a tree."""
        yield self.name
        for label in self.labels:
            yield "   " + label
        children = sorted(self.children, key=(lambda c: c.sort_by))
        last = children[-1] if children else None
        for child in children:
            prefix = "`--" if child is last else "+--"
            for line in child.tree_lines():
                yield prefix + line
                prefix = "   " if child is last else "|  "

    def __str__(self):
        return "\n".join(self.tree_lines())

    def __repr__(self):
        return f"<{self.__class__.__name__} at 0x{id(self):x}>"


def bytes2human(use_si_units=False):
    """Return the string representation of the number of bytes with units."""
    if use_si_units:
        prefixes = "TB", "GB", "MB", "kB", "B"
        values = 1e12, 1e9, 1e6, 1e3, 1
    else:
        prefixes = "TiB", "GiB", "MiB", "KiB", "B"
        values = 2**40, 2**30, 2**20, 2**10, 1

    def b2h(nbytes):

        for prefix, value in zip(prefixes, values):
            scaled = nbytes / value
            if scaled >= 1:
                break

        return f"{scaled:.1f}{prefix}"

    return b2h


def make_test_file(prefix="/tmp"):
    """Create a test file."""
    f = tb.open_file(str(Path(prefix) / "test_pttree.hdf5"), "w")

    g1 = f.create_group("/", "group1")
    g1a = f.create_group(g1, "group1a")
    g1b = f.create_group(g1, "group1b")

    filters = tb.Filters(complevel=5, complib="bzip2")

    for gg in g1a, g1b:
        f.create_carray(
            gg,
            "zeros128b",
            obj=np.zeros(32, dtype=np.float64),
            filters=filters,
        )
        f.create_carray(
            gg, "random128b", obj=np.random.rand(32), filters=filters
        )

    g2 = f.create_group("/", "group2")

    f.create_soft_link(g2, "softlink_g1_z128", "/group1/group1a/zeros128b")
    f.create_hard_link(g2, "hardlink_g1a_z128", "/group1/group1a/zeros128b")
    f.create_hard_link(g2, "hardlink_g1a", "/group1/group1a")

    return f
