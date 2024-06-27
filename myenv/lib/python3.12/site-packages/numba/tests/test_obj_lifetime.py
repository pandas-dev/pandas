import collections
import weakref
import gc
import operator
from itertools import takewhile

import unittest
from numba import njit, jit
from numba.core.compiler import CompilerBase, DefaultPassBuilder
from numba.core.untyped_passes import PreserveIR
from numba.core.typed_passes import IRLegalization
from numba.core import types, ir
from numba.tests.support import TestCase, override_config, SerialMixin


class _Dummy(object):

    def __init__(self, recorder, name):
        self.recorder = recorder
        self.name = name
        recorder._add_dummy(self)

    def __add__(self, other):
        assert isinstance(other, _Dummy)
        return _Dummy(self.recorder, "%s + %s" % (self.name, other.name))

    def __iter__(self):
        return _DummyIterator(self.recorder, "iter(%s)" % self.name)


class _DummyIterator(_Dummy):

    count = 0

    def __next__(self):
        if self.count >= 3:
            raise StopIteration
        self.count += 1
        return _Dummy(self.recorder, "%s#%s" % (self.name, self.count))

    next = __next__


class RefRecorder(object):
    """
    An object which records events when instances created through it
    are deleted.  Custom events can also be recorded to aid in
    diagnosis.
    """

    def __init__(self):
        self._counts = collections.defaultdict(int)
        self._events = []
        self._wrs = {}

    def make_dummy(self, name):
        """
        Make an object whose deletion will be recorded as *name*.
        """
        return _Dummy(self, name)

    def _add_dummy(self, dummy):
        wr = weakref.ref(dummy, self._on_disposal)
        self._wrs[wr] = dummy.name

    __call__ = make_dummy

    def mark(self, event):
        """
        Manually append *event* to the recorded events.
        *event* can be formatted using format().
        """
        count = self._counts[event] + 1
        self._counts[event] = count
        self._events.append(event.format(count=count))

    def _on_disposal(self, wr):
        name = self._wrs.pop(wr)
        self._events.append(name)

    @property
    def alive(self):
        """
        A list of objects which haven't been deleted yet.
        """
        return [wr() for wr in self._wrs]

    @property
    def recorded(self):
        """
        A list of recorded events.
        """
        return self._events


def simple_usecase1(rec):
    a = rec('a')
    b = rec('b')
    c = rec('c')
    a = b + c
    rec.mark('--1--')
    d = a + a   # b + c + b + c
    rec.mark('--2--')
    return d

def simple_usecase2(rec):
    a = rec('a')
    b = rec('b')
    rec.mark('--1--')
    x = a
    y = x
    a = None
    return y

def looping_usecase1(rec):
    a = rec('a')
    b = rec('b')
    c = rec('c')
    x = b
    for y in a:
        x = x + y
        rec.mark('--loop bottom--')
    rec.mark('--loop exit--')
    x = x + c
    return x

def looping_usecase2(rec):
    a = rec('a')
    b = rec('b')
    cum = rec('cum')
    for x in a:
        rec.mark('--outer loop top--')
        cum = cum + x
        z = x + x
        rec.mark('--inner loop entry #{count}--')
        for y in b:
            rec.mark('--inner loop top #{count}--')
            cum = cum + y
            rec.mark('--inner loop bottom #{count}--')
        rec.mark('--inner loop exit #{count}--')
        if cum:
            cum = y + z
        else:
            # Never gets here, but let the Numba compiler see a `break` opcode
            break
        rec.mark('--outer loop bottom #{count}--')
    else:
        rec.mark('--outer loop else--')
    rec.mark('--outer loop exit--')
    return cum

def generator_usecase1(rec):
    a = rec('a')
    b = rec('b')
    yield a
    yield b

def generator_usecase2(rec):
    a = rec('a')
    b = rec('b')
    for x in a:
        yield x
    yield b


class MyError(RuntimeError):
    pass

def do_raise(x):
    raise MyError(x)

def raising_usecase1(rec):
    a = rec('a')
    b = rec('b')
    d = rec('d')
    if a:
        do_raise("foo")
        c = rec('c')
        c + a
    c + b

def raising_usecase2(rec):
    a = rec('a')
    b = rec('b')
    if a:
        c = rec('c')
        do_raise(b)
    a + c

def raising_usecase3(rec):
    a = rec('a')
    b = rec('b')
    if a:
        raise MyError(b)


def del_before_definition(rec):
    """
    This test reveal a bug that there is a del on uninitialized variable
    """
    n = 5
    for i in range(n):
        rec.mark(str(i))
        n = 0
        for j in range(n):
            return 0
        else:
            if i < 2:
                continue
            elif i == 2:
                for j in range(i):
                    return i
                rec.mark('FAILED')
            rec.mark('FAILED')
        rec.mark('FAILED')
    rec.mark('OK')
    return -1


def inf_loop_multiple_back_edge(rec):
    """
    test to reveal bug of invalid liveness when infinite loop has multiple
    backedge.
    """
    while True:
        rec.mark("yield")
        yield
        p = rec('p')
        if p:
            rec.mark('bra')
            pass


class TestObjLifetime(TestCase):
    """
    Test lifetime of Python objects inside jit-compiled functions.
    """

    def compile(self, pyfunc):
        # Note: looplift must be disabled. The test require the function
        #       control-flow to be unchanged.
        cfunc = jit((types.pyobject,), forceobj=True, looplift=False)(pyfunc)
        return cfunc

    def compile_and_record(self, pyfunc, raises=None):
        rec = RefRecorder()
        cfunc = self.compile(pyfunc)
        if raises is not None:
            with self.assertRaises(raises):
                cfunc(rec)
        else:
            cfunc(rec)
        return rec

    def assertRecordOrder(self, rec, expected):
        """
        Check that the *expected* markers occur in that order in *rec*'s
        recorded events.
        """
        actual = []
        recorded = rec.recorded
        remaining = list(expected)
        # Find out in which order, if any, the expected events were recorded
        for d in recorded:
            if d in remaining:
                actual.append(d)
                # User may or may not expect duplicates, handle them properly
                remaining.remove(d)
        self.assertEqual(actual, expected,
                         "the full list of recorded events is: %r" % (recorded,))

    def test_simple1(self):
        rec = self.compile_and_record(simple_usecase1)
        self.assertFalse(rec.alive)
        self.assertRecordOrder(rec, ['a', 'b', '--1--'])
        self.assertRecordOrder(rec, ['a', 'c', '--1--'])
        self.assertRecordOrder(rec, ['--1--', 'b + c', '--2--'])

    def test_simple2(self):
        rec = self.compile_and_record(simple_usecase2)
        self.assertFalse(rec.alive)
        self.assertRecordOrder(rec, ['b', '--1--', 'a'])

    def test_looping1(self):
        rec = self.compile_and_record(looping_usecase1)
        self.assertFalse(rec.alive)
        # a and b are unneeded after the loop, check they were disposed of
        self.assertRecordOrder(rec, ['a', 'b', '--loop exit--', 'c'])
        # check disposal order of iterator items and iterator
        self.assertRecordOrder(rec, ['iter(a)#1', '--loop bottom--',
                                     'iter(a)#2', '--loop bottom--',
                                     'iter(a)#3', '--loop bottom--',
                                     'iter(a)', '--loop exit--',
                                     ])

    def test_looping2(self):
        rec = self.compile_and_record(looping_usecase2)
        self.assertFalse(rec.alive)
        # `a` is disposed of after its iterator is taken
        self.assertRecordOrder(rec, ['a', '--outer loop top--'])
        # Check disposal of iterators
        self.assertRecordOrder(rec, ['iter(a)', '--outer loop else--',
                                     '--outer loop exit--'])
        self.assertRecordOrder(rec, ['iter(b)', '--inner loop exit #1--',
                                     'iter(b)', '--inner loop exit #2--',
                                     'iter(b)', '--inner loop exit #3--',
                                     ])
        # Disposal of in-loop variable `x`
        self.assertRecordOrder(rec, ['iter(a)#1', '--inner loop entry #1--',
                                     'iter(a)#2', '--inner loop entry #2--',
                                     'iter(a)#3', '--inner loop entry #3--',
                                     ])
        # Disposal of in-loop variable `z`
        self.assertRecordOrder(rec, ['iter(a)#1 + iter(a)#1',
                                     '--outer loop bottom #1--',
                                     ])

    def exercise_generator(self, genfunc):
        cfunc = self.compile(genfunc)
        # Exhaust the generator
        rec = RefRecorder()
        with self.assertRefCount(rec):
            gen = cfunc(rec)
            next(gen)
            self.assertTrue(rec.alive)
            list(gen)
            self.assertFalse(rec.alive)
        # Instantiate the generator but never iterate
        rec = RefRecorder()
        with self.assertRefCount(rec):
            gen = cfunc(rec)
            del gen
            gc.collect()
            self.assertFalse(rec.alive)
        # Stop iterating before exhaustion
        rec = RefRecorder()
        with self.assertRefCount(rec):
            gen = cfunc(rec)
            next(gen)
            self.assertTrue(rec.alive)
            del gen
            gc.collect()
            self.assertFalse(rec.alive)

    def test_generator1(self):
        self.exercise_generator(generator_usecase1)

    def test_generator2(self):
        self.exercise_generator(generator_usecase2)

    def test_del_before_definition(self):
        rec = self.compile_and_record(del_before_definition)
        self.assertEqual(rec.recorded, ['0', '1', '2'])

    def test_raising1(self):
        with self.assertRefCount(do_raise):
            rec = self.compile_and_record(raising_usecase1, raises=MyError)
            self.assertFalse(rec.alive)

    def test_raising2(self):
        with self.assertRefCount(do_raise):
            rec = self.compile_and_record(raising_usecase2, raises=MyError)
            self.assertFalse(rec.alive)

    def test_raising3(self):
        with self.assertRefCount(MyError):
            rec = self.compile_and_record(raising_usecase3, raises=MyError)
            self.assertFalse(rec.alive)

    def test_inf_loop_multiple_back_edge(self):
        cfunc = self.compile(inf_loop_multiple_back_edge)
        rec = RefRecorder()
        iterator = iter(cfunc(rec))
        next(iterator)
        self.assertEqual(rec.alive, [])
        next(iterator)
        self.assertEqual(rec.alive, [])
        next(iterator)
        self.assertEqual(rec.alive, [])
        self.assertEqual(rec.recorded,
                         ['yield', 'p', 'bra', 'yield', 'p', 'bra', 'yield'])


class TestExtendingVariableLifetimes(SerialMixin, TestCase):
    # Test for `numba.config.EXTEND_VARIABLE_LIFETIMES` which moves the ir.Del
    # nodes to just before a block's terminator, i.e. their lifetime is extended
    # beyond the point of last use.

    def test_lifetime_basic(self):

        def get_ir(extend_lifetimes):
            class IRPreservingCompiler(CompilerBase):

                def define_pipelines(self):
                    pm = DefaultPassBuilder.define_nopython_pipeline(self.state)
                    pm.add_pass_after(PreserveIR, IRLegalization)
                    pm.finalize()
                    return [pm]

            @njit(pipeline_class=IRPreservingCompiler)
            def foo():
                a = 10
                b = 20
                c = a + b
                # a and b are now unused, standard behaviour is ir.Del for them here
                d = c / c
                return d

            with override_config('EXTEND_VARIABLE_LIFETIMES', extend_lifetimes):
                foo()
                cres = foo.overloads[foo.signatures[0]]
                func_ir = cres.metadata['preserved_ir']

            return func_ir


        def check(func_ir, expect):
            # assert single block
            self.assertEqual(len(func_ir.blocks), 1)
            blk = next(iter(func_ir.blocks.values()))

            # check sequencing
            for expect_class, got_stmt in zip(expect, blk.body):
                self.assertIsInstance(got_stmt, expect_class)

        del_after_use_ir = get_ir(False)
        # should be 3 assigns (a, b, c), 2 del (a, b), assign (d), del (c)
        # assign for cast d to return, del (d), return
        expect = [*((ir.Assign,) * 3), ir.Del, ir.Del, ir.Assign, ir.Del,
                  ir.Assign, ir.Del, ir.Return]
        check(del_after_use_ir, expect)

        del_at_block_end_ir = get_ir(True)
        # should be 4 assigns (a, b, c, d), assign for cast d to return,
        # 4 dels (a, b, c, d) then the return.
        expect = [*((ir.Assign,) * 4), ir.Assign, *((ir.Del,) * 4), ir.Return]
        check(del_at_block_end_ir, expect)

    def test_dbg_extend_lifetimes(self):

        def get_ir(**options):
            class IRPreservingCompiler(CompilerBase):

                def define_pipelines(self):
                    pm = DefaultPassBuilder.define_nopython_pipeline(self.state)
                    pm.add_pass_after(PreserveIR, IRLegalization)
                    pm.finalize()
                    return [pm]

            @njit(pipeline_class=IRPreservingCompiler, **options)
            def foo():
                a = 10
                b = 20
                c = a + b
                # a and b are now unused, standard behaviour is ir.Del for them here
                d = c / c
                return d

            foo()
            cres = foo.overloads[foo.signatures[0]]
            func_ir = cres.metadata['preserved_ir']

            return func_ir

        # _dbg_extend_lifetimes is on when debug=True
        ir_debug = get_ir(debug=True)
        # explicitly turn on _dbg_extend_lifetimes
        ir_debug_ext = get_ir(debug=True, _dbg_extend_lifetimes=True)
        # explicitly turn off _dbg_extend_lifetimes
        ir_debug_no_ext = get_ir(debug=True, _dbg_extend_lifetimes=False)

        def is_del_grouped_at_the_end(fir):
            [blk] = fir.blocks.values()
            # Mark all statements that are ir.Del
            inst_is_del = [isinstance(stmt, ir.Del) for stmt in blk.body]
            # Get the leading segment that are not dels
            not_dels = list(takewhile(operator.not_, inst_is_del))
            # Compute the starting position of the dels
            begin = len(not_dels)
            # Get the remaining segment that are all dels
            all_dels = list(takewhile(operator.truth, inst_is_del[begin:]))
            # Compute the ending position of the dels
            end = begin + len(all_dels)
            # If the dels are all grouped at the end (before the terminator),
            # the end position will be the last position of the list
            return end == len(inst_is_del) - 1

        self.assertTrue(is_del_grouped_at_the_end(ir_debug))
        self.assertTrue(is_del_grouped_at_the_end(ir_debug_ext))
        self.assertFalse(is_del_grouped_at_the_end(ir_debug_no_ext))


if __name__ == "__main__":
    unittest.main()
