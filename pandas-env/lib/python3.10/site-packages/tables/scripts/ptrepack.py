"""This utility lets you repack your data files in a flexible way.

Pass the flag -h to this for help on usage.

"""

import argparse
import sys
import warnings
from pathlib import Path
from time import perf_counter as clock
from time import process_time as cpuclock

import tables as tb


# Global variables
verbose = False
regoldindexes = True
createsysattrs = True

numpy_aliases = [
    'numeric',
    'Numeric',
    'numarray',
    'NumArray',
    'CharArray',
]


def newdst_group(dstfileh, dstgroup, title, filters):
    group = dstfileh.root
    # Now, create the new group. This works even if dstgroup == '/'
    for nodename in dstgroup.split('/'):
        if nodename == '':
            continue
        # First try if possible intermediate groups does already exist.
        try:
            group2 = dstfileh.get_node(group, nodename)
        except tb.exceptions.NoSuchNodeError:
            # The group does not exist. Create it.
            group2 = dstfileh.create_group(group, nodename,
                                           title=title,
                                           filters=filters)
        group = group2
    return group


def recreate_indexes(table, dstfileh, dsttable):
    listoldindexes = table._listoldindexes
    if listoldindexes != []:
        if not regoldindexes:
            if verbose:
                print("[I]Not regenerating indexes for table: '%s:%s'" %
                      (dstfileh.filename, dsttable._v_pathname))
            return
        # Now, recreate the indexed columns
        if verbose:
            print("[I]Regenerating indexes for table: '%s:%s'" %
                  (dstfileh.filename, dsttable._v_pathname))
        for colname in listoldindexes:
            if verbose:
                print("[I]Indexing column: '%s'. Please wait..." % colname)
            colobj = dsttable.cols._f_col(colname)
            # We don't specify the filters for the indexes
            colobj.create_index(filters=None)


def copy_leaf(srcfile, dstfile, srcnode, dstnode, title,
              filters, copyuserattrs, overwritefile, overwrtnodes, stats,
              start, stop, step, chunkshape, sortby, check_CSI,
              propindexes, upgradeflavors, allow_padding):
    # Open the source file
    srcfileh = tb.open_file(srcfile, 'r', allow_padding=allow_padding)
    # Get the source node (that should exist)
    srcnode = srcfileh.get_node(srcnode)

    # Get the destination node and its parent
    last_slash = dstnode.rindex('/')
    if last_slash == len(dstnode)-1:
        # print("Detected a trailing slash in destination node. "
        #       "Interpreting it as a destination group.")
        dstgroup = dstnode[:-1]
    elif last_slash > 0:
        dstgroup = dstnode[:last_slash]
    else:
        dstgroup = "/"
    dstleaf = dstnode[last_slash + 1:]
    if dstleaf == "":
        dstleaf = srcnode.name
    # Check whether the destination group exists or not
    if Path(dstfile).is_file() and not overwritefile:
        dstfileh = tb.open_file(dstfile, 'a',
                                pytables_sys_attrs=createsysattrs,
                                allow_padding=allow_padding)
        try:
            dstgroup = dstfileh.get_node(dstgroup)
        except Exception:
            # The dstgroup does not seem to exist. Try creating it.
            dstgroup = newdst_group(dstfileh, dstgroup, title, filters)
        else:
            # The node exists, but it is really a group?
            if not isinstance(dstgroup, tb.group.Group):
                # No. Should we overwrite it?
                if overwrtnodes:
                    parent = dstgroup._v_parent
                    last_slash = dstgroup._v_pathname.rindex('/')
                    dstgroupname = dstgroup._v_pathname[last_slash + 1:]
                    dstgroup.remove()
                    dstgroup = dstfileh.create_group(parent, dstgroupname,
                                                     title=title,
                                                     filters=filters)
                else:
                    raise RuntimeError("Please check that the node names are "
                                       "not duplicated in destination, and "
                                       "if so, add the --overwrite-nodes "
                                       "flag if desired.")
    else:
        # The destination file does not exist or will be overwritten.
        dstfileh = tb.open_file(dstfile, 'w', title=title, filters=filters,
                                pytables_sys_attrs=createsysattrs,
                                allow_padding=allow_padding)
        dstgroup = newdst_group(dstfileh, dstgroup, title="", filters=filters)

    # Finally, copy srcnode to dstnode
    try:
        dstnode = srcnode.copy(
            dstgroup, dstleaf, filters=filters,
            copyuserattrs=copyuserattrs, overwrite=overwrtnodes,
            stats=stats, start=start, stop=stop, step=step,
            chunkshape=chunkshape,
            sortby=sortby, check_CSI=check_CSI, propindexes=propindexes)
    except Exception:
        (type_, value, traceback) = sys.exc_info()
        print("Problems doing the copy from '%s:%s' to '%s:%s'" %
              (srcfile, srcnode, dstfile, dstnode))
        print(f"The error was --> {type_}: {value}")
        print("The destination file looks like:\n", dstfileh)
        # Close all the open files:
        srcfileh.close()
        dstfileh.close()
        raise RuntimeError("Please check that the node names are not "
                           "duplicated in destination, and if so, add "
                           "the --overwrite-nodes flag if desired.")

    # Upgrade flavors in dstnode, if required
    if upgradeflavors:
        if srcfileh.format_version.startswith("1"):
            # Remove original flavor in case the source file has 1.x format
            dstnode.del_attr('FLAVOR')
        elif srcfileh.format_version < "2.1":
            if dstnode.get_attr('FLAVOR') in numpy_aliases:
                dstnode.set_attr('FLAVOR', tb.flavor.internal_flavor)

    # Recreate possible old indexes in destination node
    if srcnode._c_classid == "TABLE":
        recreate_indexes(srcnode, dstfileh, dstnode)

    # Close all the open files:
    srcfileh.close()
    dstfileh.close()


def copy_children(srcfile, dstfile, srcgroup, dstgroup, title,
                  recursive, filters, copyuserattrs, overwritefile,
                  overwrtnodes, stats, start, stop, step,
                  chunkshape, sortby, check_CSI, propindexes,
                  upgradeflavors, allow_padding, use_hardlinks=True):
    """Copy the children from source group to destination group"""
    # Open the source file with srcgroup as root_uep
    srcfileh = tb.open_file(srcfile, 'r', root_uep=srcgroup,
                            allow_padding=allow_padding)
    #  Assign the root to srcgroup
    srcgroup = srcfileh.root

    created_dstgroup = False
    # Check whether the destination group exists or not
    if Path(dstfile).is_file() and not overwritefile:
        dstfileh = tb.open_file(dstfile, 'a',
                                pytables_sys_attrs=createsysattrs,
                                allow_padding=allow_padding)
        try:
            dstgroup = dstfileh.get_node(dstgroup)
        except tb.exceptions.NoSuchNodeError:
            # The dstgroup does not seem to exist. Try creating it.
            dstgroup = newdst_group(dstfileh, dstgroup, title, filters)
            created_dstgroup = True
        else:
            # The node exists, but it is really a group?
            if not isinstance(dstgroup, tb.group.Group):
                # No. Should we overwrite it?
                if overwrtnodes:
                    parent = dstgroup._v_parent
                    last_slash = dstgroup._v_pathname.rindex('/')
                    dstgroupname = dstgroup._v_pathname[last_slash + 1:]
                    dstgroup.remove()
                    dstgroup = dstfileh.create_group(parent, dstgroupname,
                                                     title=title,
                                                     filters=filters)
                else:
                    raise RuntimeError("Please check that the node names are "
                                       "not duplicated in destination, and "
                                       "if so, add the --overwrite-nodes "
                                       "flag if desired.")
    else:
        # The destination file does not exist or will be overwritten.
        dstfileh = tb.open_file(dstfile, 'w', title=title, filters=filters,
                                pytables_sys_attrs=createsysattrs,
                                allow_padding=allow_padding)
        dstgroup = newdst_group(dstfileh, dstgroup, title="", filters=filters)
        created_dstgroup = True

    # Copy the attributes to dstgroup, if needed
    if created_dstgroup and copyuserattrs:
        srcgroup._v_attrs._f_copy(dstgroup)

    # Finally, copy srcgroup children to dstgroup
    try:
        srcgroup._f_copy_children(
            dstgroup, recursive=recursive, filters=filters,
            copyuserattrs=copyuserattrs, overwrite=overwrtnodes,
            stats=stats, start=start, stop=stop, step=step,
            chunkshape=chunkshape,
            sortby=sortby, check_CSI=check_CSI, propindexes=propindexes,
            use_hardlinks=use_hardlinks)
    except Exception:
        (type_, value, traceback) = sys.exc_info()
        print("Problems doing the copy from '%s:%s' to '%s:%s'" %
              (srcfile, srcgroup, dstfile, dstgroup))
        print(f"The error was --> {type_}: {value}")
        print("The destination file looks like:\n", dstfileh)
        # Close all the open files:
        srcfileh.close()
        dstfileh.close()
        raise RuntimeError("Please check that the node names are not "
                           "duplicated in destination, and if so, add the "
                           "--overwrite-nodes flag if desired. In "
                           "particular, pay attention that root_uep is not "
                           "fooling you.")

    # Upgrade flavors in dstnode, if required
    if upgradeflavors:
        for dstnode in dstgroup._f_walknodes("Leaf"):
            if srcfileh.format_version.startswith("1"):
                # Remove original flavor in case the source file has 1.x format
                dstnode.del_attr('FLAVOR')
            elif srcfileh.format_version < "2.1":
                if dstnode.get_attr('FLAVOR') in numpy_aliases:
                    dstnode.set_attr('FLAVOR', tb.flavor.internal_flavor)

    # Convert the remaining tables with old indexes (if any)
    for table in srcgroup._f_walknodes("Table"):
        dsttable = dstfileh.get_node(dstgroup, table._v_pathname)
        recreate_indexes(table, dstfileh, dsttable)

    # Close all the open files:
    srcfileh.close()
    dstfileh.close()


def _get_parser():
    parser = argparse.ArgumentParser(
        description='''This utility is very powerful and lets you copy any
        leaf, group or complete subtree into another file.
        During the copy process you are allowed to change the filter
        properties if you want so. Also, in the case of duplicated pathnames,
        you can decide if you want to overwrite already existing nodes on the
        destination file. Generally speaking, ptrepack can be useful in may
        situations, like replicating a subtree in another file, change the
        filters in objects and see how affect this to the compression degree
        or I/O performance, consolidating specific data in repositories or
        even *importing* generic HDF5 files and create true PyTables
        counterparts.''')

    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='show verbose information',
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', dest='overwritefile',
        help='overwrite destination file',
    )
    parser.add_argument(
        '-R', '--range', dest='rng', metavar='RANGE',
        help='''select a RANGE of rows (in the form "start,stop,step")
        during the copy of *all* the leaves.
        Default values are "None,None,1", which means a copy of all the
        rows.''',
    )
    parser.add_argument(
        '--non-recursive', action='store_false', default=True,
        dest='recursive',
        help='do not do a recursive copy. Default is to do it',
    )
    parser.add_argument(
        '--dest-title', dest='title', default='',
        help='title for the new file (if not specified, the source is copied)',
    )
    parser.add_argument(
        '--dont-create-sysattrs', action='store_false', default=True,
        dest='createsysattrs',
        help='do not create sys attrs (default is to do it)',
    )
    parser.add_argument(
        '--dont-copy-userattrs', action='store_false', default=True,
        dest='copyuserattrs',
        help='do not copy the user attrs (default is to do it)',
    )
    parser.add_argument(
        '--overwrite-nodes', action='store_true', dest='overwrtnodes',
        help='''overwrite destination nodes if they exist.
        Default is to not overwrite them''',
    )
    parser.add_argument(
        '--complevel', type=int, default=0,
        help='''set a compression level (0 for no compression, which is the
        default)''',
    )
    parser.add_argument(
        '--complib', choices=(
            "zlib", "lzo", "bzip2", "blosc", "blosc:blosclz",
            "blosc:lz4", "blosc:lz4hc", "blosc:zlib", "blosc:zstd",
            "blosc2", "blosc2:blosclz",
            "blosc2:lz4", "blosc2:lz4hc", "blosc2:zlib", "blosc2:zstd"
        ),
        default='zlib',
        help='''set the compression library to be used during the copy.
        Defaults to %(default)s''',
    )
    parser.add_argument(
        '--shuffle', type=int, choices=(0, 1),
        help='''activate or not the shuffle filter (default is active if
        complevel > 0)''',
    )
    parser.add_argument(
        '--bitshuffle', type=int, choices=(0, 1),
        help='activate or not the bitshuffle filter (not active by default)',
    )
    parser.add_argument(
        '--fletcher32', type=int, choices=(0, 1),
        help='''whether to activate or not the fletcher32 filter (not active
        by default)''',
    )
    parser.add_argument(
        '--keep-source-filters', action='store_true', dest='keepfilters',
        help='''use the original filters in source files.
        The default is not doing that if any of --complevel, --complib,
        --shuffle --bitshuffle or --fletcher32 option is specified''',
    )
    parser.add_argument(
        '--chunkshape', default='keep',
        help='''set a chunkshape.
        Possible options are: "keep" | "auto" | int | tuple.
        A value of "auto" computes a sensible value for the chunkshape of the
        leaves copied.  The default is to "keep" the original value''',
    )
    parser.add_argument(
        '--upgrade-flavors', action='store_true', dest='upgradeflavors',
        help='''when repacking PyTables 1.x or PyTables 2.x files, the flavor
        of leaves will be unset. With this, such a leaves will be serialized
        as objects with the internal flavor ('numpy' for 3.x series)''',
    )
    parser.add_argument(
        '--dont-regenerate-old-indexes', action='store_false', default=True,
        dest='regoldindexes',
        help='''disable regenerating old indexes.
        The default is to regenerate old indexes as they are found''',
    )
    parser.add_argument(
        '--sortby', metavar='COLUMN',
        help='''do a table copy sorted by the index in "column".
        For reversing the order, use a negative value in the "step" part of
        "RANGE" (see "-r" flag).  Only applies to table objects''',
    )
    parser.add_argument(
        '--checkCSI', action='store_true',
        help='force the check for a CSI index for the --sortby column',
    )
    parser.add_argument(
        '--propindexes', action='store_true',
        help='''propagate the indexes existing in original tables. The default
        is to not propagate them.  Only applies to table objects''',
    )
    parser.add_argument(
        '--dont-allow-padding', action='store_true',
        dest="dont_allow_padding",
        help='''remove the possible padding in compound types in source files.
        The default is to propagate it.  Only applies to table objects''',
    )
    parser.add_argument(
        'src', metavar='sourcefile:sourcegroup', help='source file/group',
    )
    parser.add_argument(
        'dst', metavar='destfile:destgroup', help='destination file/group',
    )

    return parser


def main():
    global verbose
    global regoldindexes
    global createsysattrs

    parser = _get_parser()
    args = parser.parse_args()

    # check arguments
    if args.rng:
        try:
            args.rng = eval("slice(" + args.rng + ")")
        except Exception:
            parser.error("Error when getting the range parameter.")

    if args.chunkshape.isdigit() or args.chunkshape.startswith('('):
        args.chunkshape = eval(args.chunkshape)

    if args.complevel < 0 or args.complevel > 9:
        parser.error(
            'invalid "complevel" value, it sould be in te range [0, 9]'
        )

    # Catch the files passed as the last arguments
    src = args.src.rsplit(':', 1)
    dst = args.dst.rsplit(':', 1)
    if len(src) == 1:
        srcfile, srcnode = src[0], "/"
    else:
        srcfile, srcnode = src
    if len(dst) == 1:
        dstfile, dstnode = dst[0], "/"
    else:
        dstfile, dstnode = dst

    if srcnode == "":
        # case where filename == "filename:" instead of "filename:/"
        srcnode = "/"

    if dstnode == "":
        # case where filename == "filename:" instead of "filename:/"
        dstnode = "/"

    # Ignore the warnings for tables that contains oldindexes
    # (these will be handled by the copying routines)
    warnings.filterwarnings("ignore", category=tb.exceptions.OldIndexWarning)

    # Ignore the flavors warnings during upgrading flavor operations
    if args.upgradeflavors:
        warnings.filterwarnings("ignore", category=tb.exceptions.FlavorWarning)

    # Build the Filters instance
    filter_params = (
        args.complevel,
        args.complib,
        args.shuffle,
        args.bitshuffle,
        args.fletcher32,
    )
    if (filter_params == (None,) * 4 or args.keepfilters):
        filters = None
    else:
        if args.complevel is None:
            args.complevel = 0
        if args.shuffle is None:
            if args.complevel > 0:
                args.shuffle = True
            else:
                args.shuffle = False
        if args.bitshuffle is None:
            args.bitshuffle = False
        if args.bitshuffle:
            # Shuffle and bitshuffle are mutually exclusive
            args.shuffle = False
        if args.complib is None:
            args.complib = "zlib"
        if args.fletcher32 is None:
            args.fletcher32 = False
        filters = tb.leaf.Filters(complevel=args.complevel,
                                  complib=args.complib, shuffle=args.shuffle,
                                  bitshuffle=args.bitshuffle,
                                  fletcher32=args.fletcher32)

    # The start, stop and step params:
    start, stop, step = None, None, 1  # Defaults
    if args.rng:
        start, stop, step = args.rng.start, args.rng.stop, args.rng.step

    # Set globals
    verbose = args.verbose
    regoldindexes = args.regoldindexes
    createsysattrs = args.createsysattrs

    # Some timing
    t1 = clock()
    cpu1 = cpuclock()
    # Copy the file
    if verbose:
        print("+=+" * 20)
        print("Recursive copy:", args.recursive)
        print("Applying filters:", filters)
        if args.sortby is not None:
            print("Sorting table(s) by column:", args.sortby)
            print("Forcing a CSI creation:", args.checkCSI)
        if args.propindexes:
            print("Recreating indexes in copied table(s)")
        print(f"Start copying {srcfile}:{srcnode} to {dstfile}:{dstnode}")
        print("+=+" * 20)

    allow_padding = not args.dont_allow_padding
    # Check whether the specified source node is a group or a leaf
    h5srcfile = tb.open_file(srcfile, 'r', allow_padding=allow_padding)
    srcnodeobject = h5srcfile.get_node(srcnode)

    # Close the file again
    h5srcfile.close()

    stats = {'groups': 0, 'leaves': 0, 'links': 0, 'bytes': 0, 'hardlinks': 0}
    if isinstance(srcnodeobject, tb.group.Group):
        copy_children(
            srcfile, dstfile, srcnode, dstnode,
            title=args.title, recursive=args.recursive, filters=filters,
            copyuserattrs=args.copyuserattrs, overwritefile=args.overwritefile,
            overwrtnodes=args.overwrtnodes, stats=stats,
            start=start, stop=stop, step=step, chunkshape=args.chunkshape,
            sortby=args.sortby, check_CSI=args.checkCSI,
            propindexes=args.propindexes,
            upgradeflavors=args.upgradeflavors,
            allow_padding=allow_padding,
            use_hardlinks=True)
    else:
        # If not a Group, it should be a Leaf
        copy_leaf(
            srcfile, dstfile, srcnode, dstnode,
            title=args.title, filters=filters,
            copyuserattrs=args.copyuserattrs,
            overwritefile=args.overwritefile, overwrtnodes=args.overwrtnodes,
            stats=stats, start=start, stop=stop, step=step,
            chunkshape=args.chunkshape,
            sortby=args.sortby, check_CSI=args.checkCSI,
            propindexes=args.propindexes,
            upgradeflavors=args.upgradeflavors,
            allow_padding=allow_padding,
        )

    # Gather some statistics
    t2 = clock()
    cpu2 = cpuclock()
    tcopy = t2 - t1
    cpucopy = cpu2 - cpu1
    if verbose:
        ngroups = stats['groups']
        nleaves = stats['leaves']
        nlinks = stats['links']
        nhardlinks = stats['hardlinks']
        nbytescopied = stats['bytes']
        nnodes = ngroups + nleaves + nlinks + nhardlinks

        print(
            "Groups copied:", ngroups,
            ", Leaves copied:", nleaves,
            ", Links copied:", nlinks,
            ", Hard links copied:", nhardlinks,
        )
        if args.copyuserattrs:
            print("User attrs copied")
        else:
            print("User attrs not copied")
        print(f"KBytes copied: {nbytescopied / 1024:.3f}")
        print(
            f"Time copying: {tcopy:.3f} s (real) {cpucopy:.3f} s "
            f"(cpu)  {cpucopy / tcopy:.0%}")
        print(f"Copied nodes/sec: {nnodes / tcopy:.1f}")
        print(f"Copied KB/s : {nbytescopied / tcopy / 1024:.0f}")
