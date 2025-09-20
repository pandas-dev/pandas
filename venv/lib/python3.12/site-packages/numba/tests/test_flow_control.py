import itertools

import unittest
from numba import jit
from numba.core.controlflow import CFGraph, ControlFlowAnalysis
from numba.core import types
from numba.core.bytecode import FunctionIdentity, ByteCode, _fix_LOAD_GLOBAL_arg
from numba.tests.support import TestCase

enable_pyobj_flags = {}

no_pyobj_flags = {'_nrt': True, 'nopython': True}


def for_loop_usecase1(x, y):
    result = 0
    for i in range(x):
        result += i
    return result


def for_loop_usecase2(x, y):
    result = 0
    for i, j in enumerate(range(x, y, -1)):
        result += i * j
    return result


def for_loop_usecase3(x, y):
    result = 0
    lst = [x, y]
    for i in lst:
        result += i
    return result


def for_loop_usecase4(x, y):
    result = 0
    for i in range(10):
        for j in range(10):
            result += 1
    return result


def for_loop_usecase5(x, y):
    result = 0
    for i in range(x):
        result += 1
        if result > y:
            break
    return result


def for_loop_usecase6(x, y):
    result = 0
    for i in range(x):
        if i > y:
            continue
        result += 1
    return result


def for_loop_usecase7(x, y):
    for i in range(x):
        x = 0
        for j in range(x):
            return 1
        else:
            pass
    return 0


def for_loop_usecase8(x, y):
    result = 0
    for i in range(x, y, y - x + 1):
        result += 1
    return result


def for_loop_usecase9(x, y):
    z = 0
    for i in range(x):
        x = 0
        for j in range(x):
            if j == x / 2:
                z += j
                break
        else:
            z += y

    return z


def for_loop_usecase10(x, y):
    for i in range(x):
        if i == y:
            z = y
            break
    else:
        z = i * 2
    return z


def while_loop_usecase1(x, y):
    result = 0
    i = 0
    while i < x:
        result += i
        i += 1
    return result


def while_loop_usecase2(x, y):
    result = 0
    while result != x:
        result += 1
    return result


def while_loop_usecase3(x, y):
    result = 0
    i = 0
    j = 0
    while i < x:
        while j < y:
            result += i + j
            i += 1
            j += 1
    return result


def while_loop_usecase4(x, y):
    result = 0
    while True:
        result += 1
        if result > x:
            break
    return result


def while_loop_usecase5(x, y):
    result = 0
    while result < x:
        if result > y:
            result += 2
            continue
        result += 1
    return result


def ifelse_usecase1(x, y):
    if x > 0:
        pass
    elif y > 0:
        pass
    else:
        pass
    return True


def ifelse_usecase2(x, y):
    if x > y:
        return 1
    elif x == 0 or y == 0:
        return 2
    else:
        return 3


def ifelse_usecase3(x, y):
    if x > 0:
        if y > 0:
            return 1
        elif y < 0:
            return 1
        else:
            return 0
    elif x < 0:
        return 1
    else:
        return 0


def ifelse_usecase4(x, y):
    if x == y:
        return 1


def ternary_ifelse_usecase1(x, y):
    return True if x > y else False


def double_infinite_loop(x, y):
    L = x
    i = y

    while True:
        while True:
            if i == L - 1:
                break
            i += 1
        i += 1
        if i >= L:
            break

    return i, L


def try_except_usecase():
    try:
        pass
    except Exception:
        pass


class TestFlowControl(TestCase):

    def run_test(self, pyfunc, x_operands, y_operands,
                 flags=enable_pyobj_flags):
        cfunc = jit((types.intp, types.intp), **flags)(pyfunc)
        for x, y in itertools.product(x_operands, y_operands):
            pyerr = None
            cerr = None
            try:
                pyres = pyfunc(x, y)
            except Exception as e:
                pyerr = e

            try:
                cres = cfunc(x, y)
            except Exception as e:
                if pyerr is None:
                    raise
                cerr = e
                self.assertEqual(type(pyerr), type(cerr))
            else:
                if pyerr is not None:
                    self.fail("Invalid for pure-python but numba works\n" +
                              pyerr)
                self.assertEqual(pyres, cres)

    def test_for_loop1(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase1, [-10, 0, 10], [0], flags=flags)

    def test_for_loop1_npm(self):
        self.test_for_loop1(flags=no_pyobj_flags)

    def test_for_loop2(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase2, [-10, 0, 10], [-10, 0, 10],
                      flags=flags)

    def test_for_loop2_npm(self):
        self.test_for_loop2(flags=no_pyobj_flags)

    def test_for_loop3(self, flags=enable_pyobj_flags):
        """
        List requires pyobject
        """
        self.run_test(for_loop_usecase3, [1], [2],
                      flags=flags)

    def test_for_loop3_npm(self):
        self.test_for_loop3(flags=no_pyobj_flags)

    def test_for_loop4(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase4, [10], [10], flags=flags)

    def test_for_loop4_npm(self):
        self.test_for_loop4(flags=no_pyobj_flags)

    def test_for_loop5(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase5, [100], [50], flags=flags)

    def test_for_loop5_npm(self):
        self.test_for_loop5(flags=no_pyobj_flags)

    def test_for_loop6(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase6, [100], [50], flags=flags)

    def test_for_loop6_npm(self):
        self.test_for_loop6(flags=no_pyobj_flags)

    def test_for_loop7(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase7, [5], [0], flags=flags)

    def test_for_loop7_npm(self):
        self.test_for_loop7(flags=no_pyobj_flags)

    def test_for_loop8(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase8, [0, 1], [0, 2, 10], flags=flags)

    def test_for_loop8_npm(self):
        self.test_for_loop8(flags=no_pyobj_flags)

    def test_for_loop9(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase9, [0, 1], [0, 2, 10], flags=flags)

    def test_for_loop9_npm(self):
        self.test_for_loop9(flags=no_pyobj_flags)

    def test_for_loop10(self, flags=enable_pyobj_flags):
        self.run_test(for_loop_usecase10, [5], [2, 7], flags=flags)

    def test_for_loop10_npm(self):
        self.test_for_loop10(flags=no_pyobj_flags)

    def test_while_loop1(self, flags=enable_pyobj_flags):
        self.run_test(while_loop_usecase1, [10], [0], flags=flags)

    def test_while_loop1_npm(self):
        self.test_while_loop1(flags=no_pyobj_flags)

    def test_while_loop2(self, flags=enable_pyobj_flags):
        self.run_test(while_loop_usecase2, [10], [0], flags=flags)

    def test_while_loop2_npm(self):
        self.test_while_loop2(flags=no_pyobj_flags)

    def test_while_loop3(self, flags=enable_pyobj_flags):
        self.run_test(while_loop_usecase3, [10], [10], flags=flags)

    def test_while_loop3_npm(self):
        self.test_while_loop3(flags=no_pyobj_flags)

    def test_while_loop4(self, flags=enable_pyobj_flags):
        self.run_test(while_loop_usecase4, [10], [0], flags=flags)

    def test_while_loop4_npm(self):
        self.test_while_loop4(flags=no_pyobj_flags)

    def test_while_loop5(self, flags=enable_pyobj_flags):
        self.run_test(while_loop_usecase5, [0, 5, 10], [0, 5, 10], flags=flags)

    def test_while_loop5_npm(self):
        self.test_while_loop5(flags=no_pyobj_flags)

    def test_ifelse1(self, flags=enable_pyobj_flags):
        self.run_test(ifelse_usecase1, [-1, 0, 1], [-1, 0, 1], flags=flags)

    def test_ifelse1_npm(self):
        self.test_ifelse1(flags=no_pyobj_flags)

    def test_ifelse2(self, flags=enable_pyobj_flags):
        self.run_test(ifelse_usecase2, [-1, 0, 1], [-1, 0, 1], flags=flags)

    def test_ifelse2_npm(self):
        self.test_ifelse2(flags=no_pyobj_flags)

    def test_ifelse3(self, flags=enable_pyobj_flags):
        self.run_test(ifelse_usecase3, [-1, 0, 1], [-1, 0, 1], flags=flags)

    def test_ifelse3_npm(self):
        self.test_ifelse3(flags=no_pyobj_flags)

    def test_ifelse4(self, flags=enable_pyobj_flags):
        self.run_test(ifelse_usecase4, [-1, 0, 1], [-1, 0, 1], flags=flags)

    def test_ifelse4_npm(self):
        self.test_ifelse4(flags=no_pyobj_flags)

    def test_ternary_ifelse1(self, flags=enable_pyobj_flags):
        self.run_test(ternary_ifelse_usecase1, [-1, 0, 1], [-1, 0, 1],
                      flags=flags)

    def test_ternary_ifelse1_npm(self):
        self.test_ternary_ifelse1(flags=no_pyobj_flags)

    def test_double_infinite_loop(self, flags=enable_pyobj_flags):
        self.run_test(double_infinite_loop, [10], [0],
                      flags=flags)

    def test_double_infinite_loop_npm(self):
        self.test_double_infinite_loop(flags=no_pyobj_flags)


class TestCFGraph(TestCase):
    """
    Test the numba.controlflow.CFGraph class.
    """

    def from_adj_list(self, d, entry_point=0):
        """
        Build a CFGraph class from a dict of adjacency lists.
        """
        g = CFGraph()
        # Need to add all nodes before adding edges
        for node in d:
            g.add_node(node)
        for node, dests in d.items():
            for dest in dests:
                g.add_edge(node, dest)
        return g

    def loopless1(self):
        """
        A simple CFG corresponding to the following code structure:

            c = (... if ... else ...) + ...
            return b + c
        """
        g = self.from_adj_list({0: [18, 12], 12: [21], 18: [21], 21: []})
        g.set_entry_point(0)
        g.process()
        return g

    def loopless1_dead_nodes(self):
        """
        Same as loopless1(), but with added dead blocks (some of them
        in a loop).
        """
        g = self.from_adj_list(
            {0: [18, 12],
             12: [21],
             18: [21],
             21: [],
             91: [12, 0],
             92: [91, 93],
             93: [92],
             94: [],
             })
        g.set_entry_point(0)
        g.process()
        return g

    def loopless2(self):
        """
        A loopless CFG corresponding to the following code structure:

            c = (... if ... else ...) + ...
            if c:
                return ...
            else:
                return ...

        Note there are two exit points, and the entry point has been
        changed to a non-zero value.
        """
        g = self.from_adj_list(
            {99: [18, 12], 12: [21], 18: [21], 21: [42, 34], 34: [], 42: []})
        g.set_entry_point(99)
        g.process()
        return g

    def multiple_loops(self):
        """
        A CFG with multiple nested loops:

            for y in b:
                for x in a:
                    # This loop has two back edges
                    if b:
                        continue
                    else:
                        continue
            for z in c:
                if z:
                    return ...
        """
        g = self.from_adj_list({0: [7],
                                7: [10, 60],
                                10: [13],
                                13: [20],
                                20: [56, 23],
                                23: [32, 44],
                                32: [20],
                                44: [20],
                                56: [57],
                                57: [7],
                                60: [61],
                                61: [68],
                                68: [87, 71],
                                71: [80, 68],
                                80: [],
                                87: [88],
                                88: []}
                               )
        g.set_entry_point(0)
        g.process()
        return g

    def multiple_exits(self):
        """
        A CFG with three loop exits, one of which is also a function
        exit point, and another function exit point:

            for x in a:
                if a:
                    return b
                elif b:
                    break
            return c
        """
        g = self.from_adj_list(
            {0: [7],
             7: [10, 36],
             10: [19, 23],
             19: [],
             23: [29, 7],
             29: [37],
             36: [37],
             37: []
             })
        g.set_entry_point(0)
        g.process()
        return g

    def infinite_loop1(self):
        """
        A CFG with a infinite loop and an alternate exit point:

            if c:
                return
            while True:
                if a:
                    ...
                else:
                    ...
        """
        g = self.from_adj_list(
            {0: [10, 6], 6: [], 10: [13], 13: [26, 19], 19: [13], 26: [13]})
        g.set_entry_point(0)
        g.process()
        return g

    def infinite_loop2(self):
        """
        A CFG with no exit point at all:

            while True:
                if a:
                    ...
                else:
                    ...
        """
        g = self.from_adj_list({0: [3], 3: [16, 9], 9: [3], 16: [3]})
        g.set_entry_point(0)
        g.process()
        return g

    def test_simple_properties(self):
        g = self.loopless1()
        self.assertEqual(sorted(g.successors(0)), [(12, None), (18, None)])
        self.assertEqual(sorted(g.successors(21)), [])
        self.assertEqual(sorted(g.predecessors(0)), [])
        self.assertEqual(sorted(g.predecessors(21)), [(12, None), (18, None)])

    def test_exit_points(self):
        g = self.loopless1()
        self.assertEqual(sorted(g.exit_points()), [21])
        g = self.loopless1_dead_nodes()
        self.assertEqual(sorted(g.exit_points()), [21])
        g = self.loopless2()
        self.assertEqual(sorted(g.exit_points()), [34, 42])
        g = self.multiple_loops()
        self.assertEqual(sorted(g.exit_points()), [80, 88])
        g = self.infinite_loop1()
        self.assertEqual(sorted(g.exit_points()), [6])
        g = self.infinite_loop2()
        self.assertEqual(sorted(g.exit_points()), [])
        g = self.multiple_exits()
        self.assertEqual(sorted(g.exit_points()), [19, 37])

    def test_dead_nodes(self):
        g = self.loopless1()
        self.assertEqual(len(g.dead_nodes()), 0)
        self.assertEqual(sorted(g.nodes()),
                         [0, 12, 18, 21])
        g = self.loopless2()
        self.assertEqual(len(g.dead_nodes()), 0)
        self.assertEqual(sorted(g.nodes()),
                         [12, 18, 21, 34, 42, 99])
        g = self.multiple_loops()
        self.assertEqual(len(g.dead_nodes()), 0)
        g = self.infinite_loop1()
        self.assertEqual(len(g.dead_nodes()), 0)
        g = self.multiple_exits()
        self.assertEqual(len(g.dead_nodes()), 0)
        # Only this example has dead nodes
        g = self.loopless1_dead_nodes()
        self.assertEqual(sorted(g.dead_nodes()),
                         [91, 92, 93, 94])
        self.assertEqual(sorted(g.nodes()),
                         [0, 12, 18, 21])

    def test_descendents(self):
        g = self.loopless2()
        d = g.descendents(34)
        self.assertEqual(sorted(d), [])
        d = g.descendents(42)
        self.assertEqual(sorted(d), [])
        d = g.descendents(21)
        self.assertEqual(sorted(d), [34, 42])
        d = g.descendents(99)
        self.assertEqual(sorted(d), [12, 18, 21, 34, 42])
        g = self.infinite_loop1()
        d = g.descendents(26)
        self.assertEqual(sorted(d), [])
        d = g.descendents(19)
        self.assertEqual(sorted(d), [])
        d = g.descendents(13)
        self.assertEqual(sorted(d), [19, 26])
        d = g.descendents(10)
        self.assertEqual(sorted(d), [13, 19, 26])
        d = g.descendents(6)
        self.assertEqual(sorted(d), [])
        d = g.descendents(0)
        self.assertEqual(sorted(d), [6, 10, 13, 19, 26])

    def test_topo_order(self):
        g = self.loopless1()
        self.assertIn(g.topo_order(),
                      ([0, 12, 18, 21], [0, 18, 12, 21]))
        g = self.loopless2()
        self.assertIn(g.topo_order(),
                      ([99, 18, 12, 21, 34, 42], [99, 12, 18, 21, 34, 42]))
        g = self.infinite_loop2()
        self.assertIn(g.topo_order(),
                      ([0, 3, 9, 16], [0, 3, 16, 9]))
        g = self.infinite_loop1()
        self.assertIn(g.topo_order(),
                      ([0, 6, 10, 13, 19, 26], [0, 6, 10, 13, 26, 19],
                       [0, 10, 13, 19, 26, 6], [0, 10, 13, 26, 19, 6]))

    def test_topo_sort(self):
        def check_topo_sort(nodes, expected):
            self.assertIn(list(g.topo_sort(nodes)), expected)
            self.assertIn(list(g.topo_sort(nodes[::-1])), expected)
            self.assertIn(list(g.topo_sort(nodes, reverse=True))[::-1],
                          expected)
            self.assertIn(list(g.topo_sort(nodes[::-1], reverse=True))[::-1],
                          expected)
            self.random.shuffle(nodes)
            self.assertIn(list(g.topo_sort(nodes)), expected)
            self.assertIn(list(g.topo_sort(nodes, reverse=True))[::-1],
                          expected)

        g = self.loopless2()
        check_topo_sort([21, 99, 12, 34], ([99, 12, 21, 34],))
        # NOTE: topo_sort() is not stable
        check_topo_sort([18, 12, 42, 99],
                        ([99, 12, 18, 42], [99, 18, 12, 42]))
        g = self.multiple_exits()
        check_topo_sort([19, 10, 7, 36],
                        ([7, 10, 19, 36], [7, 10, 36, 19], [7, 36, 10, 19]))

    def check_dominators(self, got, expected):
        self.assertEqual(sorted(got), sorted(expected))
        for node in sorted(got):
            self.assertEqual(sorted(got[node]), sorted(expected[node]),
                             "mismatch for %r" % (node,))

    def test_dominators_loopless(self):
        def eq_(d, l):
            self.assertEqual(sorted(doms[d]), l)
        for g in [self.loopless1(), self.loopless1_dead_nodes()]:
            doms = g.dominators()
            eq_(0, [0])
            eq_(12, [0, 12])
            eq_(18, [0, 18])
            eq_(21, [0, 21])
        g = self.loopless2()
        doms = g.dominators()
        eq_(99, [99])
        eq_(12, [12, 99])
        eq_(18, [18, 99])
        eq_(21, [21, 99])
        eq_(34, [21, 34, 99])
        eq_(42, [21, 42, 99])

    def test_dominators_loops(self):
        g = self.multiple_exits()
        doms = g.dominators()
        self.check_dominators(doms,
                              {0: [0],
                               7: [0, 7],
                               10: [0, 7, 10],
                               19: [0, 7, 10, 19],
                               23: [0, 7, 10, 23],
                               29: [0, 7, 10, 23, 29],
                               36: [0, 7, 36],
                               37: [0, 7, 37],
                               })
        g = self.multiple_loops()
        doms = g.dominators()
        self.check_dominators(doms,
                              {0: [0],
                               7: [0, 7],
                               10: [0, 10, 7],
                               13: [0, 10, 13, 7],
                               20: [0, 10, 20, 13, 7],
                               23: [0, 20, 23, 7, 10, 13],
                               32: [32, 0, 20, 23, 7, 10, 13],
                               44: [0, 20, 23, 7, 10, 44, 13],
                               56: [0, 20, 7, 56, 10, 13],
                               57: [0, 20, 7, 56, 57, 10, 13],
                               60: [0, 60, 7],
                               61: [0, 60, 61, 7],
                               68: [0, 68, 60, 61, 7],
                               71: [0, 68, 71, 7, 60, 61],
                               80: [80, 0, 68, 71, 7, 60, 61],
                               87: [0, 68, 87, 7, 60, 61],
                               88: [0, 68, 87, 88, 7, 60, 61]
                               })
        g = self.infinite_loop1()
        doms = g.dominators()
        self.check_dominators(doms,
                              {0: [0],
                               6: [0, 6],
                               10: [0, 10],
                               13: [0, 10, 13],
                               19: [0, 10, 19, 13],
                               26: [0, 10, 13, 26],
                               })

    def test_post_dominators_loopless(self):
        def eq_(d, l):
            self.assertEqual(sorted(doms[d]), l)
        for g in [self.loopless1(), self.loopless1_dead_nodes()]:
            doms = g.post_dominators()
            eq_(0, [0, 21])
            eq_(12, [12, 21])
            eq_(18, [18, 21])
            eq_(21, [21])
        g = self.loopless2()
        doms = g.post_dominators()
        eq_(34, [34])
        eq_(42, [42])
        eq_(21, [21])
        eq_(18, [18, 21])
        eq_(12, [12, 21])
        eq_(99, [21, 99])

    def test_post_dominators_loops(self):
        g = self.multiple_exits()
        doms = g.post_dominators()
        self.check_dominators(doms,
                              {0: [0, 7],
                               7: [7],
                               10: [10],
                               19: [19],
                               23: [23],
                               29: [29, 37],
                               36: [36, 37],
                               37: [37],
                               })
        g = self.multiple_loops()
        doms = g.post_dominators()
        self.check_dominators(doms,
                              {0: [0, 60, 68, 61, 7],
                               7: [60, 68, 61, 7],
                               10: [68, 7, 10, 13, 20, 56, 57, 60, 61],
                               13: [68, 7, 13, 20, 56, 57, 60, 61],
                               20: [20, 68, 7, 56, 57, 60, 61],
                               23: [68, 7, 20, 23, 56, 57, 60, 61],
                               32: [32, 68, 7, 20, 56, 57, 60, 61],
                               44: [68, 7, 44, 20, 56, 57, 60, 61],
                               56: [68, 7, 56, 57, 60, 61],
                               57: [57, 60, 68, 61, 7],
                               60: [60, 68, 61],
                               61: [68, 61],
                               68: [68],
                               71: [71],
                               80: [80],
                               87: [88, 87],
                               88: [88]
                               })

    def test_post_dominators_infinite_loops(self):
        # Post-dominators with infinite loops need special care
        # (the ordinary algorithm won't work).
        g = self.infinite_loop1()
        doms = g.post_dominators()
        self.check_dominators(doms,
                              {0: [0],
                               6: [6],
                               10: [10, 13],
                               13: [13],
                               19: [19],
                               26: [26],
                               })
        g = self.infinite_loop2()
        doms = g.post_dominators()
        self.check_dominators(doms,
                              {0: [0, 3],
                               3: [3],
                               9: [9],
                               16: [16],
                               })

    def test_dominator_tree(self):
        def check(graph, expected):
            domtree = graph.dominator_tree()
            self.assertEqual(domtree, expected)

        check(self.loopless1(),
              {0: {12, 18, 21}, 12: set(), 18: set(), 21: set()})
        check(self.loopless2(),
              {12: set(), 18: set(), 21: {34, 42}, 34: set(), 42: set(),
               99: {18, 12, 21}})
        check(self.loopless1_dead_nodes(),
              {0: {12, 18, 21}, 12: set(), 18: set(), 21: set()})
        check(self.multiple_loops(),
              {0: {7}, 7: {10, 60}, 60: {61}, 61: {68}, 68: {71, 87},
               87: {88}, 88: set(), 71: {80}, 80: set(), 10: {13},
               13: {20}, 20: {56, 23}, 23: {32, 44}, 44: set(),
               32: set(), 56: {57}, 57: set()})
        check(self.multiple_exits(),
              {0: {7}, 7: {10, 36, 37}, 36: set(), 10: {19, 23},
               23: {29}, 29: set(), 37: set(), 19: set()})
        check(self.infinite_loop1(),
              {0: {10, 6}, 6: set(), 10: {13}, 13: {26, 19}, 19: set(),
               26: set()})
        check(self.infinite_loop2(),
              {0: {3}, 3: {16, 9}, 9: set(), 16: set()})

    def test_immediate_dominators(self):
        def check(graph, expected):
            idoms = graph.immediate_dominators()
            self.assertEqual(idoms, expected)

        check(self.loopless1(),
              {0: 0, 12: 0, 18: 0, 21: 0})
        check(self.loopless2(),
              {18: 99, 12: 99, 21: 99, 42: 21, 34: 21, 99: 99})
        check(self.loopless1_dead_nodes(),
              {0: 0, 12: 0, 18: 0, 21: 0})
        check(self.multiple_loops(),
              {0: 0, 7: 0, 10: 7, 13: 10, 20: 13, 23: 20,
               32: 23, 44: 23, 56: 20, 57: 56, 60: 7, 61: 60,
               68: 61, 71: 68, 80: 71, 87: 68, 88: 87})
        check(self.multiple_exits(),
              {0:0, 7: 0, 10: 7, 19: 10, 23: 10, 29: 23, 36: 7, 37: 7})
        check(self.infinite_loop1(),
              {0: 0, 6: 0, 10: 0, 13: 10, 19: 13, 26: 13})
        check(self.infinite_loop2(),
              {0: 0, 3: 0, 9: 3, 16: 3})

    def test_dominance_frontier(self):
        def check(graph, expected):
            df = graph.dominance_frontier()
            self.assertEqual(df, expected)

        check(self.loopless1(),
              {0: set(), 12: {21}, 18: {21}, 21: set()})
        check(self.loopless2(),
              {18: {21}, 12: {21}, 21: set(), 42: set(), 34: set(), 99: set()})
        check(self.loopless1_dead_nodes(),
              {0: set(), 12: {21}, 18: {21}, 21: set()})
        check(self.multiple_loops(),
              {0: set(), 7: {7}, 10: {7}, 13: {7}, 20: {20, 7}, 23: {20},
               32: {20}, 44: {20}, 56: {7}, 57: {7}, 60: set(), 61: set(),
               68: {68}, 71: {68}, 80: set(), 87: set(), 88: set()})
        check(self.multiple_exits(),
              {0: set(), 7: {7}, 10: {37, 7}, 19: set(),
               23: {37, 7}, 29: {37}, 36: {37}, 37: set()})
        check(self.infinite_loop1(),
              {0: set(), 6: set(), 10: set(), 13: {13}, 19: {13}, 26: {13}})
        check(self.infinite_loop2(),
              {0: set(), 3: {3}, 9: {3}, 16: {3}})

    def test_backbone_loopless(self):
        for g in [self.loopless1(), self.loopless1_dead_nodes()]:
            self.assertEqual(sorted(g.backbone()), [0, 21])
        g = self.loopless2()
        self.assertEqual(sorted(g.backbone()), [21, 99])

    def test_backbone_loops(self):
        g = self.multiple_loops()
        self.assertEqual(sorted(g.backbone()), [0, 7, 60, 61, 68])
        g = self.infinite_loop1()
        self.assertEqual(sorted(g.backbone()), [0])
        g = self.infinite_loop2()
        self.assertEqual(sorted(g.backbone()), [0, 3])

    def test_loops(self):
        for g in [self.loopless1(), self.loopless1_dead_nodes(),
                  self.loopless2()]:
            self.assertEqual(len(g.loops()), 0)

        g = self.multiple_loops()
        # Loop headers
        self.assertEqual(sorted(g.loops()), [7, 20, 68])
        outer1 = g.loops()[7]
        inner1 = g.loops()[20]
        outer2 = g.loops()[68]
        self.assertEqual(outer1.header, 7)
        self.assertEqual(sorted(outer1.entries), [0])
        self.assertEqual(sorted(outer1.exits), [60])
        self.assertEqual(sorted(outer1.body),
                         [7, 10, 13, 20, 23, 32, 44, 56, 57])
        self.assertEqual(inner1.header, 20)
        self.assertEqual(sorted(inner1.entries), [13])
        self.assertEqual(sorted(inner1.exits), [56])
        self.assertEqual(sorted(inner1.body), [20, 23, 32, 44])
        self.assertEqual(outer2.header, 68)
        self.assertEqual(sorted(outer2.entries), [61])
        self.assertEqual(sorted(outer2.exits), [80, 87])
        self.assertEqual(sorted(outer2.body), [68, 71])
        for node in [0, 60, 61, 80, 87, 88]:
            self.assertEqual(g.in_loops(node), [])
        for node in [7, 10, 13, 56, 57]:
            self.assertEqual(g.in_loops(node), [outer1])
        for node in [20, 23, 32, 44]:
            self.assertEqual(g.in_loops(node), [inner1, outer1])
        for node in [68, 71]:
            self.assertEqual(g.in_loops(node), [outer2])

        g = self.infinite_loop1()
        # Loop headers
        self.assertEqual(sorted(g.loops()), [13])
        loop = g.loops()[13]
        self.assertEqual(loop.header, 13)
        self.assertEqual(sorted(loop.entries), [10])
        self.assertEqual(sorted(loop.exits), [])
        self.assertEqual(sorted(loop.body), [13, 19, 26])
        for node in [0, 6, 10]:
            self.assertEqual(g.in_loops(node), [])
        for node in [13, 19, 26]:
            self.assertEqual(g.in_loops(node), [loop])

        g = self.infinite_loop2()
        # Loop headers
        self.assertEqual(sorted(g.loops()), [3])
        loop = g.loops()[3]
        self.assertEqual(loop.header, 3)
        self.assertEqual(sorted(loop.entries), [0])
        self.assertEqual(sorted(loop.exits), [])
        self.assertEqual(sorted(loop.body), [3, 9, 16])
        for node in [0]:
            self.assertEqual(g.in_loops(node), [])
        for node in [3, 9, 16]:
            self.assertEqual(g.in_loops(node), [loop])

        g = self.multiple_exits()
        # Loop headers
        self.assertEqual(sorted(g.loops()), [7])
        loop = g.loops()[7]
        self.assertEqual(loop.header, 7)
        self.assertEqual(sorted(loop.entries), [0])
        self.assertEqual(sorted(loop.exits), [19, 29, 36])
        self.assertEqual(sorted(loop.body), [7, 10, 23])
        for node in [0, 19, 29, 36]:
            self.assertEqual(g.in_loops(node), [])
        for node in [7, 10, 23]:
            self.assertEqual(g.in_loops(node), [loop])

    def test_loop_dfs_pathological(self):
        # The follow adjlist is an export from the reproducer in #6186
        g = self.from_adj_list({
            0: {38, 14},
            14: {38, 22},
            22: {38, 30},
            30: {42, 38},
            38: {42},
            42: {64, 50},
            50: {64, 58},
            58: {128},
            64: {72, 86},
            72: {80, 86},
            80: {128},
            86: {108, 94},
            94: {108, 102},
            102: {128},
            108: {128, 116},
            116: {128, 124},
            124: {128},
            128: {178, 174},
            174: {178},
            178: {210, 206},
            206: {210},
            210: {248, 252},
            248: {252},
            252: {282, 286},
            282: {286},
            286: {296, 326},
            296: {330},
            326: {330},
            330: {370, 340},
            340: {374},
            370: {374},
            374: {380, 382},
            380: {382},
            382: {818, 390},
            390: {456, 458},
            456: {458},
            458: {538, 566},
            538: {548, 566},
            548: set(),
            566: {586, 572},
            572: {586},
            586: {708, 596},
            596: {608},
            608: {610},
            610: {704, 620},
            620: {666, 630},
            630: {636, 646},
            636: {666, 646},
            646: {666},
            666: {610},
            704: {706},
            706: {818},
            708: {720},
            720: {722},
            722: {816, 732},
            732: {778, 742},
            742: {748, 758},
            748: {778, 758},
            758: {778},
            778: {722},
            816: {818},
            818: set(),
        })
        g.set_entry_point(0)
        g.process()
        stats = {}
        # Compute backedges and store the iteration count for testing
        back_edges = g._find_back_edges(stats=stats)
        self.assertEqual(back_edges, {(666, 610), (778, 722)})
        self.assertEqual(stats['iteration_count'], 155)

    def test_equals(self):

        def get_new():
            g = self.from_adj_list({0: [18, 12], 12: [21], 18: [21], 21: []})
            g.set_entry_point(0)
            g.process()
            return g

        x = get_new()
        y = get_new()

        # identical
        self.assertEqual(x, y)

        # identical but defined in a different order
        g = self.from_adj_list({0: [12, 18], 18: [21], 21: [], 12: [21]})
        g.set_entry_point(0)
        g.process()
        self.assertEqual(x, g)

        # different entry point
        z = get_new()
        z.set_entry_point(18)
        z.process()
        self.assertNotEqual(x, z)

        # extra node/edge, same entry point
        z = self.from_adj_list({0: [18, 12], 12: [21], 18: [21], 21: [22],
                                22: []})
        z.set_entry_point(0)
        z.process()
        self.assertNotEqual(x, z)

        # same nodes, different edges
        a = self.from_adj_list({0: [18, 12], 12: [0], 18: []})
        a.set_entry_point(0)
        a.process()
        z = self.from_adj_list({0: [18, 12], 12: [18], 18: []})
        z.set_entry_point(0)
        z.process()
        self.assertNotEqual(a, z)


class TestRealCodeDomFront(TestCase):
    """Test IDOM and DOMFRONT computation on real python bytecode.
    Note: there will be less testing on IDOM (esp in loop) because of
    the extra blocks inserted by the interpreter.  But, testing on DOMFRONT
    (which depends on IDOM) is easier.

    Testing is done by associating names to basicblock by using globals of
    the pattern "SET_BLOCK_<name>", which are scanned by
    `.get_cfa_and_namedblocks` into *namedblocks* dictionary.  That way, we
    can check that a block of a certain name is a IDOM or DOMFRONT of another
    named block.
    """
    def cfa(self, bc):
        cfa = ControlFlowAnalysis(bc)
        cfa.run()
        return cfa

    def get_cfa_and_namedblocks(self, fn):
        fid = FunctionIdentity.from_function(fn)
        bc = ByteCode(func_id=fid)
        cfa = self.cfa(bc)
        namedblocks = self._scan_namedblocks(bc, cfa)

        #### To debug, uncomment below
        # print(bc.dump())
        # print("IDOMS")
        # for k, v in sorted(cfa.graph.immediate_dominators().items()):
        #     print('{} -> {}'.format(k, v))
        # print("DOMFRONT")
        # for k, vs in sorted(cfa.graph.dominance_frontier().items()):
        #     print('{} -> {}'.format(k, vs))
        # print(namedblocks)
        # cfa.graph.render_dot().view()

        return cfa, namedblocks

    def _scan_namedblocks(self, bc, cfa):
        """Scan namedblocks as denoted by a LOAD_GLOBAL bytecode referring
        to global variables with the pattern "SET_BLOCK_<name>", where "<name>"
        would be the name for the current block.
        """
        namedblocks = {}
        blocks = sorted([x.offset for x in cfa.iterblocks()])
        prefix = 'SET_BLOCK_'

        for inst in bc:
            # Find LOAD_GLOBAL that refers to "SET_BLOCK_<name>"
            if inst.opname == 'LOAD_GLOBAL':
                gv = bc.co_names[_fix_LOAD_GLOBAL_arg(inst.arg)]
                if gv.startswith(prefix):
                    name = gv[len(prefix):]
                    # Find the block where this instruction resides
                    for s, e in zip(blocks, blocks[1:] + [blocks[-1] + 1]):
                        if s <= inst.offset < e:
                            break
                    else:
                        raise AssertionError('unreachable loop')
                    blkno = s
                    namedblocks[name] = blkno
        return namedblocks

    def test_loop(self):
        def foo(n):
            c = 0
            SET_BLOCK_A                     # noqa: F821
            i = 0
            while SET_BLOCK_B0:             # noqa: F821
                SET_BLOCK_B1                # noqa: F821
                c += 1
                i += 1
            SET_BLOCK_C                     # noqa: F821
            return c

        cfa, blkpts = self.get_cfa_and_namedblocks(foo)

        # Py3.10 turns while loop into if(...) { do {...} while(...) }.
        # Also, `SET_BLOCK_B0` is duplicated. As a result, the second B0
        # is picked up by `blkpts`.
        domfront = cfa.graph.dominance_frontier()
        self.assertFalse(domfront[blkpts['A']])
        self.assertFalse(domfront[blkpts['C']])

    def test_loop_nested_and_break(self):
        def foo(n):
            SET_BLOCK_A                     # noqa: F821
            while SET_BLOCK_B0:             # noqa: F821
                SET_BLOCK_B1                # noqa: F821
                while SET_BLOCK_C0:         # noqa: F821
                    SET_BLOCK_C1            # noqa: F821
                    if SET_BLOCK_D0:        # noqa: F821
                        SET_BLOCK_D1        # noqa: F821
                        break
                    elif n:
                        SET_BLOCK_D2        # noqa: F821
                    SET_BLOCK_E             # noqa: F821
                SET_BLOCK_F                 # noqa: F821
            SET_BLOCK_G                     # noqa: F821

        cfa, blkpts = self.get_cfa_and_namedblocks(foo)
        self.assertEqual(blkpts['D0'], blkpts['C1'])

        # Py3.10 changes while loop into if-do-while
        domfront = cfa.graph.dominance_frontier()
        self.assertFalse(domfront[blkpts['A']])
        self.assertFalse(domfront[blkpts['G']])
        # 2 domfront members for C1
        # C0 because of the loop; F because of the break.
        self.assertEqual({blkpts['F']}, domfront[blkpts['D1']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['D2']])

    def test_if_else(self):
        def foo(a, b):
            c = 0
            SET_BLOCK_A           # noqa: F821
            if a < b:
                SET_BLOCK_B       # noqa: F821
                c = 1
            elif SET_BLOCK_C0:    # noqa: F821
                SET_BLOCK_C1      # noqa: F821
                c = 2
            else:
                SET_BLOCK_D       # noqa: F821
                c = 3

            SET_BLOCK_E           # noqa: F821
            if a % b == 0:
                SET_BLOCK_F       # noqa: F821
                c += 1
            SET_BLOCK_G           # noqa: F821
            return c

        cfa, blkpts = self.get_cfa_and_namedblocks(foo)

        idoms = cfa.graph.immediate_dominators()
        self.assertEqual(blkpts['A'], idoms[blkpts['B']])
        self.assertEqual(blkpts['A'], idoms[blkpts['C0']])
        self.assertEqual(blkpts['C0'], idoms[blkpts['C1']])
        self.assertEqual(blkpts['C0'], idoms[blkpts['D']])
        self.assertEqual(blkpts['A'], idoms[blkpts['E']])
        self.assertEqual(blkpts['E'], idoms[blkpts['F']])
        self.assertEqual(blkpts['E'], idoms[blkpts['G']])

        domfront = cfa.graph.dominance_frontier()
        self.assertFalse(domfront[blkpts['A']])
        self.assertFalse(domfront[blkpts['E']])
        self.assertFalse(domfront[blkpts['G']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['B']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['C0']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['C1']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['D']])
        self.assertEqual({blkpts['G']}, domfront[blkpts['F']])

    def test_if_else_nested(self):
        def foo():
            if SET_BLOCK_A0:                # noqa: F821
                SET_BLOCK_A1                # noqa: F821
                if SET_BLOCK_B0:            # noqa: F821
                    SET_BLOCK_B1            # noqa: F821
                    a = 0
                else:
                    if SET_BLOCK_C0:        # noqa: F821
                        SET_BLOCK_C1        # noqa: F821
                        a = 1
                    else:
                        SET_BLOCK_C2        # noqa: F821
                        a = 2
                    SET_BLOCK_D             # noqa: F821
                SET_BLOCK_E                 # noqa: F821
            SET_BLOCK_F                     # noqa: F821
            return a

        cfa, blkpts = self.get_cfa_and_namedblocks(foo)

        idoms = cfa.graph.immediate_dominators()
        self.assertEqual(blkpts['A0'], idoms[blkpts['A1']])
        self.assertEqual(blkpts['A1'], idoms[blkpts['B1']])
        self.assertEqual(blkpts['A1'], idoms[blkpts['C0']])
        self.assertEqual(blkpts['C0'], idoms[blkpts['D']])
        self.assertEqual(blkpts['A1'], idoms[blkpts['E']])
        self.assertEqual(blkpts['A0'], idoms[blkpts['F']])

        domfront = cfa.graph.dominance_frontier()
        self.assertFalse(domfront[blkpts['A0']])
        self.assertFalse(domfront[blkpts['F']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['B1']])
        self.assertEqual({blkpts['D']}, domfront[blkpts['C1']])
        self.assertEqual({blkpts['E']}, domfront[blkpts['D']])
        self.assertEqual({blkpts['F']}, domfront[blkpts['E']])

    def test_infinite_loop(self):
        def foo():
            SET_BLOCK_A                     # noqa: F821
            while True:  # infinite loop
                if SET_BLOCK_B:             # noqa: F821
                    SET_BLOCK_C             # noqa: F821
                    return
                SET_BLOCK_D                 # noqa: F821
            SET_BLOCK_E                     # noqa: F821

        cfa, blkpts = self.get_cfa_and_namedblocks(foo)

        idoms = cfa.graph.immediate_dominators()
        # Py3.10 optimizes away the infinite loop and removes SET_BLOCK_E from
        # the bytecode.
        self.assertNotIn('E', blkpts)
        self.assertEqual(blkpts['B'], idoms[blkpts['C']])
        self.assertEqual(blkpts['B'], idoms[blkpts['D']])

        domfront = cfa.graph.dominance_frontier()
        self.assertFalse(domfront[blkpts['A']])
        self.assertFalse(domfront[blkpts['C']])
        self.assertEqual({blkpts['B']}, domfront[blkpts['B']])
        self.assertEqual({blkpts['B']}, domfront[blkpts['D']])


if __name__ == '__main__':
    unittest.main()
