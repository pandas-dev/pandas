"""This utility lets you look into the data and metadata of your data files.

Pass the flag -h to this for help on usage.

"""

import argparse
import operator

import tables as tb

# default options
options = argparse.Namespace(
    rng=slice(None),
    showattrs=0,
    verbose=0,
    dump=0,
    colinfo=0,
    idxinfo=0,
)


def dump_leaf(leaf):
    """Dump an HDF5 leaf node."""
    if options.verbose:
        print(repr(leaf))
    else:
        print(str(leaf))
    if options.showattrs:
        print(f"  {leaf.attrs!r}")
    if options.dump and not isinstance(leaf, tb.unimplemented.UnImplemented):
        print("  Data dump:")
        # print((leaf.read(options.rng.start, options.rng.stop,
        #        options.rng.step))
        # This is better for large objects
        if options.rng.start is None:
            start = 0
        else:
            start = options.rng.start
        if options.rng.stop is None:
            if leaf.shape != ():
                stop = leaf.shape[0]
        else:
            stop = options.rng.stop
        if options.rng.step is None:
            step = 1
        else:
            step = options.rng.step
        if leaf.shape == ():
            print("[SCALAR] %s" % (leaf[()]))
        else:
            for i in range(start, stop, step):
                print(f"[{i}] {leaf[i]}")

    if isinstance(leaf, tb.table.Table) and options.colinfo:
        # Show info of columns
        for colname in leaf.colnames:
            print(repr(leaf.cols._f_col(colname)))

    if isinstance(leaf, tb.table.Table) and options.idxinfo:
        # Show info of indexes
        for colname in leaf.colnames:
            col = leaf.cols._f_col(colname)
            if isinstance(col, tb.table.Column) and col.index is not None:
                idx = col.index
                print(repr(idx))


def dump_group(pgroup, sort=False):
    """Dump an HDF5 group."""
    node_kinds = pgroup._v_file._node_kinds[1:]
    what = pgroup._f_walk_groups()
    if sort:
        what = sorted(what, key=operator.attrgetter("_v_pathname"))
    for group in what:
        print(str(group))
        if options.showattrs:
            print(f"  {group._v_attrs!r}")
        for kind in node_kinds:
            for node in group._f_list_nodes(kind):
                if options.verbose or options.dump:
                    dump_leaf(node)
                else:
                    print(str(node))


def _get_parser():
    parser = argparse.ArgumentParser(
        description=(
            "The ptdump utility allows you look into the contents of your "
            "PyTables files. It lets you see not only the data but also "
            "the metadata (that is, the *structure* and additional "
            "information in the form of *attributes*)."
        )
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="dump more metainformation on nodes",
    )
    parser.add_argument(
        "-d",
        "--dump",
        action="store_true",
        help="dump data information on leaves",
    )
    parser.add_argument(
        "-a",
        "--showattrs",
        action="store_true",
        help="show attributes in nodes (only useful when -v or -d are active)",
    )
    parser.add_argument(
        "-s",
        "--sort",
        action="store_true",
        help="sort output by node name",
    )
    parser.add_argument(
        "-c",
        "--colinfo",
        action="store_true",
        help="""show info of columns in tables (only useful when -v or -d
        are active)""",
    )
    parser.add_argument(
        "-i",
        "--idxinfo",
        action="store_true",
        help="""show info of indexed columns (only useful when -v or -d are
        active)""",
    )
    parser.add_argument(
        "-R",
        "--range",
        dest="rng",
        metavar="RANGE",
        help=(
            "select a RANGE of rows (in the form 'start,stop,step') "
            "during the copy of *all* the leaves. "
            "Default values are 'None, None, 1', which means a copy of all "
            "the rows."
        ),
    )
    parser.add_argument(
        "src",
        metavar="filename[:nodepath]",
        help="name of the HDF5 file to dump",
    )

    return parser


def main():
    """Implement the main CLI interface."""
    parser = _get_parser()

    args = parser.parse_args(namespace=options)

    # Get the options
    if isinstance(args.rng, str):
        try:
            options.rng = eval("slice(" + args.rng + ")")
        except Exception:
            parser.error("Error when getting the range parameter.")
        else:
            args.dump = 1

    # Catch the files passed as the last arguments
    src = args.src.rsplit(":", 1)
    if len(src) == 1:
        filename, nodename = src[0], "/"
    else:
        filename, nodename = src
        if nodename == "":
            # case where filename == "filename:" instead of "filename:/"
            nodename = "/"

    try:
        h5file = tb.open_file(filename, "r")
    except Exception as e:
        return f"Cannot open input file: {e}"

    with h5file:
        # Check whether the specified node is a group or a leaf
        nodeobject = h5file.get_node(nodename)
        if isinstance(nodeobject, tb.group.Group):
            # Close the file again and reopen using the root_uep
            dump_group(nodeobject, args.sort)
        elif isinstance(nodeobject, tb.leaf.Leaf):
            # If it is not a Group, it must be a Leaf
            dump_leaf(nodeobject)
        else:
            # This should never happen
            print("Unrecognized object:", nodeobject)
