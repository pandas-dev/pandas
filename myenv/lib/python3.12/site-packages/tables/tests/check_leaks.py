from pathlib import Path
from time import perf_counter as clock

import tables as tb

tref = clock()
trel = tref


def show_mem(explain):
    global tref, trel

    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmSize:"):
            vmsize = int(line.split()[1])
        elif line.startswith("VmRSS:"):
            vmrss = int(line.split()[1])
        elif line.startswith("VmData:"):
            vmdata = int(line.split()[1])
        elif line.startswith("VmStk:"):
            vmstk = int(line.split()[1])
        elif line.startswith("VmExe:"):
            vmexe = int(line.split()[1])
        elif line.startswith("VmLib:"):
            vmlib = int(line.split()[1])

    print("\nMemory usage: ******* %s *******" % explain)
    print(f"VmSize: {vmsize:>7} kB\tVmRSS: {vmrss:>7} kB")
    print(f"VmData: {vmdata:>7} kB\tVmStk: {vmstk:>7} kB")
    print(f"VmExe:  {vmexe:>7} kB\tVmLib: {vmlib:>7} kB")
    print("WallClock time:", clock() - tref, end=' ')
    print("  Delta time:", clock() - trel)
    trel = clock()


def write_group(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            fileh.create_group(fileh.root, 'group' + str(child),
                               "child: %d" % child)
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_group(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'group' + str(child))
            assert node is not None
            # flavor = node._v_attrs.CLASS
#         for child in fileh.walk_nodes():
#             pass
        show_mem("After reading metadata. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def write_array(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            fileh.create_array(fileh.root, 'array' + str(child),
                               [1, 1], "child: %d" % child)
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_array(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'array' + str(child))
            # flavor = node._v_attrs.FLAVOR
            data = node[:]  # Read data
            assert data is not None
        show_mem("After reading data. Iter %s" % i)
        # for child in range(nchildren):
        #     node = fileh.get_node(fileh.root, 'array' + str(child))
        #     flavor = node._v_attrs.FLAVOR
        #     # flavor = node._v_attrs
        # for child in fileh.walk_nodes():
        #     pass
        # show_mem("After reading metadata. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def write_carray(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            fileh.create_carray(fileh.root, 'array' + str(child),
                                tb.IntAtom(), (2,), "child: %d" % child)
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_carray(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'array' + str(child))
            # flavor = node._v_attrs.FLAVOR
            data = node[:]  # Read data
            assert data is not None
            # print("data-->", data)
        show_mem("After reading data. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def write_earray(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            ea = fileh.create_earray(fileh.root, 'array' + str(child),
                                     tb.IntAtom(), shape=(0,),
                                     title="child: %d" % child)
            ea.append([1, 2, 3])
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_earray(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'array' + str(child))
            # flavor = node._v_attrs.FLAVOR
            data = node[:]  # Read data
            assert data is not None
            # print("data-->", data)
        show_mem("After reading data. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def write_vlarray(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            vl = fileh.create_vlarray(fileh.root, 'array' + str(child),
                                      tb.IntAtom(), "child: %d" % child)
            vl.append([1, 2, 3])
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_vlarray(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'array' + str(child))
            # flavor = node._v_attrs.FLAVOR
            data = node[:]  # Read data
            assert data is not None
            # print("data-->", data)
        show_mem("After reading data. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def write_table(filename, nchildren, niter):

    class Record(tb.IsDescription):
        var1 = tb.IntCol(pos=1)
        var2 = tb.StringCol(length=1, pos=2)
        var3 = tb.FloatCol(pos=3)

    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            t = fileh.create_table(fileh.root, 'table' + str(child),
                                   Record, "child: %d" % child)
            t.append([[1, "2", 3.]])
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_table(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'table' + str(child))
            # klass = node._v_attrs.CLASS
            data = node[:]  # Read data
            assert data is not None
            # print("data-->", data)
        show_mem("After reading data. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def write_xtable(filename, nchildren, niter):

    class Record(tb.IsDescription):
        var1 = tb.IntCol(pos=1)
        var2 = tb.StringCol(length=1, pos=2)
        var3 = tb.FloatCol(pos=3)

    for i in range(niter):
        fileh = tb.open_file(filename, mode="w")
        for child in range(nchildren):
            t = fileh.create_table(fileh.root, 'table' + str(child),
                                   Record, "child: %d" % child)
            t.append([[1, "2", 3.]])
            t.cols.var1.create_index()
        show_mem("After creating. Iter %s" % i)
        fileh.close()
        show_mem("After close")


def read_xtable(filename, nchildren, niter):
    for i in range(niter):
        fileh = tb.open_file(filename, mode="r")
        for child in range(nchildren):
            node = fileh.get_node(fileh.root, 'table' + str(child))
            # klass = node._v_attrs.CLASS
            # data = node[:]  # Read data
            # print("data-->", data)
        show_mem("After reading data. Iter %s" % i)
        fileh.close()
        show_mem("After close")
        del node


if __name__ == '__main__':
    import pstats
    import argparse
    import profile as prof

    def _get_parser():
        parser = argparse.ArgumentParser(
            description='Check for PyTables memory leaks.')
        parser.add_argument('-v', '--verbose', action='store_true',
                            help='enable verbose mode')
        parser.add_argument('-p', '--profile', action='store_true',
                            help='profile')
        parser.add_argument('-a', '--array', action='store_true',
                            help='create/read arrays (default)')
        parser.add_argument('-c', '--carray', action='store_true',
                            help='create/read carrays')
        parser.add_argument('-e', '--earray', action='store_true',
                            help='create/read earrays')
        parser.add_argument('-l', '--vlarray', action='store_true',
                            help='create/read vlarrays')
        parser.add_argument('-t', '--table', action='store_true',
                            help='create/read tables')
        parser.add_argument('-x', '--indexed-table', action='store_true',
                            dest='xtable', help='create/read indexed-tables')
        parser.add_argument('-g', '--group', action='store_true',
                            help='create/read groups')
        parser.add_argument('-r', '--read', action='store_true',
                            help='only read test')
        parser.add_argument('-w', '--write', action='store_true',
                            help='only write test')
        parser.add_argument('-n', '--nchildren', type=int, default=1000,
                            help='number of children (%(default)d is the '
                                 'default)')
        parser.add_argument('-i', '--niter', type=int, default=3,
                            help='number of iterations (default: %(default)d)')

        parser.add_argument('filename', help='HDF5 file name')

        return parser

    parser = _get_parser()
    args = parser.parse_args()

    # set 'array' as default value if no ather option has been specified
    for name in ('carray', 'earray', 'vlarray', 'table', 'xtable', 'group'):
        if getattr(args, name):
            break
    else:
        args.array = True

    filename = args.filename
    nchildren = args.nchildren
    niter = args.niter

    if args.array:
        fwrite = 'write_array'
        fread = 'read_array'
    elif args.carray:
        fwrite = 'write_carray'
        fread = 'read_carray'
    elif args.earray:
        fwrite = 'write_earray'
        fread = 'read_earray'
    elif args.vlarray:
        fwrite = 'write_vlarray'
        fread = 'read_vlarray'
    elif args.table:
        fwrite = 'write_table'
        fread = 'read_table'
    elif args.xtable:
        fwrite = 'write_xtable'
        fread = 'read_xtable'
    elif args.group:
        fwrite = 'write_group'
        fread = 'read_group'

    show_mem("Before open")
    if args.write:
        if args.profile:
            prof.run(str(fwrite)+'(filename, nchildren, niter)',
                     'write_file.prof')
            stats = pstats.Stats('write_file.prof')
            stats.strip_dirs()
            stats.sort_stats('time', 'calls')
            if args.verbose:
                stats.print_stats()
            else:
                stats.print_stats(20)
        else:
            eval(fwrite+'(filename, nchildren, niter)')
    if args.read:
        if args.profile:
            prof.run(fread+'(filename, nchildren, niter)', 'read_file.prof')
            stats = pstats.Stats('read_file.prof')
            stats.strip_dirs()
            stats.sort_stats('time', 'calls')
            if args.verbose:
                print('profile -verbose')
                stats.print_stats()
            else:
                stats.print_stats(20)
        else:
            eval(fread+'(filename, nchildren, niter)')
