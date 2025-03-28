import unittest
import warnings
from contextlib import contextmanager

import numpy as np

import llvmlite.binding as llvm

from numba import njit, types
from numba.core.errors import NumbaInvalidConfigWarning
from numba.core.codegen import _parse_refprune_flags
from numba.tests.support import override_config, TestCase


@contextmanager
def set_refprune_flags(flags):
    with override_config('LLVM_REFPRUNE_FLAGS', flags):
        yield


class TestRefOpPruning(TestCase):

    _numba_parallel_test_ = False

    def check(self, func, *argtys, **prune_types):
        """
        Asserts the the func compiled with argument types "argtys" reports
        refop pruning statistics. The **prune_types** kwargs list each kind
        of pruning and whether the stat should be zero (False) or >0 (True).

        Note: The exact statistic varies across platform.

        NOTE: Tests using this `check` method need to run in subprocesses as
        `njit` sets up the module pass manager etc once and the overrides have
        no effect else.
        """

        with override_config('LLVM_REFPRUNE_PASS', '1'):
            cres = njit((*argtys,))(func).overloads[(*argtys,)]

        pstats = cres.metadata.get('prune_stats', None)
        self.assertIsNotNone(pstats)

        for k, v in prune_types.items():
            stat = getattr(pstats, k, None)
            self.assertIsNotNone(stat)
            msg = f'failed checking {k}'
            if v:
                self.assertGreater(stat, 0, msg=msg)
            else:
                self.assertEqual(stat, 0, msg=msg)

    @TestCase.run_test_in_subprocess
    def test_basic_block_1(self):
        # some nominally involved control flow and ops, there's only basic_block
        # opportunities present here.
        def func(n):
            a = np.zeros(n)
            acc = 0
            if n > 4:
                b = a[1:]
                acc += b[1]
            else:
                c = a[:-1]
                acc += c[0]
            return acc

        self.check(func, (types.intp), basicblock=True)

    @TestCase.run_test_in_subprocess
    def test_diamond_1(self):
        # most basic?! diamond
        def func(n):
            a = np.ones(n)
            x = 0
            if n > 2:
                x = a.sum()
            return x + 1

        # disable fanout pruning
        with set_refprune_flags('per_bb,diamond'):
            self.check(func, (types.intp), basicblock=True, diamond=True,
                       fanout=False, fanout_raise=False)

    @TestCase.run_test_in_subprocess
    def test_diamond_2(self):
        # more complex diamonds
        def func(n):
            con = []
            for i in range(n):
                con.append(np.arange(i))
            c = 0.0
            for arr in con:
                c += arr.sum() / (1 + arr.size)
            return c

        # disable fanout pruning
        with set_refprune_flags('per_bb,diamond'):
            self.check(func, (types.intp), basicblock=True, diamond=True,
                       fanout=False, fanout_raise=False)

    @TestCase.run_test_in_subprocess
    def test_fanout_1(self):
        # most basic?! fan-out
        def func(n):
            a = np.zeros(n)
            b = np.zeros(n)
            x = (a, b)
            acc = 0.
            for i in x:
                acc += i[0]
            return acc

        self.check(func, (types.intp), basicblock=True, fanout=True)

    @TestCase.run_test_in_subprocess
    def test_fanout_2(self):
        # fanout with raise
        def func(n):
            a = np.zeros(n)
            b = np.zeros(n)
            x = (a, b)
            for i in x:
                if n:
                    raise ValueError
            return x

        with set_refprune_flags('per_bb,fanout'):
            self.check(func, (types.intp), basicblock=True, diamond=False,
                       fanout=True, fanout_raise=False)

    @TestCase.run_test_in_subprocess
    def test_fanout_3(self):
        # fanout with raise
        def func(n):
            ary = np.arange(n)
            # basically an impl of array.sum
            c = 0
            # The raise is from StopIteration of next(iterator) implicit in
            # the for loop
            for v in np.nditer(ary):
                c += v.item()
            return 1

        with set_refprune_flags('per_bb,fanout_raise'):
            self.check(func, (types.intp), basicblock=True, diamond=False,
                       fanout=False, fanout_raise=True)


class TestRefPruneFlags(TestCase):
    def setUp(self):
        warnings.simplefilter('error', NumbaInvalidConfigWarning)

    def tearDown(self):
        warnings.resetwarnings()

    def test_warn_invalid_flags(self):
        with set_refprune_flags('abc,per_bb,cde'):
            with self.assertWarns(NumbaInvalidConfigWarning) as cm:
                optval = _parse_refprune_flags()
            self.assertEqual(len(cm.warnings), 2)
            self.assertIn('abc', str(cm.warnings[0].message))
            self.assertIn('cde', str(cm.warnings[1].message))
            self.assertEqual(optval, llvm.RefPruneSubpasses.PER_BB)

    def test_valid_flag(self):
        with set_refprune_flags('per_bb, diamond, fanout,fanout_raise'):
            optval = _parse_refprune_flags()
            self.assertEqual(optval, llvm.RefPruneSubpasses.ALL)

    def test_the_all_flag(self):
        with set_refprune_flags('all'):
            optval = _parse_refprune_flags()
            self.assertEqual(optval, llvm.RefPruneSubpasses.ALL)

    def test_some_flags(self):
        with set_refprune_flags('per_bb, fanout'):
            optval = _parse_refprune_flags()
            enumcls = llvm.RefPruneSubpasses
            self.assertEqual(optval, enumcls.PER_BB | enumcls.FANOUT)


if __name__ == "__main__":
    unittest.main()
