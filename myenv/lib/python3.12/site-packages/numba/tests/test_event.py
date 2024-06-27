import unittest
import string

import numpy as np

from numba import njit, jit, literal_unroll
from numba.core import event as ev
from numba.tests.support import TestCase, override_config


class TestEvent(TestCase):

    def setUp(self):
        # Trigger compilation to ensure all listeners are initialized
        njit(lambda: None)()
        self.__registered_listeners = len(ev._registered)

    def tearDown(self):
        # Check there is no lingering listeners
        self.assertEqual(len(ev._registered), self.__registered_listeners)

    def test_recording_listener(self):
        @njit
        def foo(x):
            return x + x

        with ev.install_recorder("numba:compile") as rec:
            foo(1)

        self.assertIsInstance(rec, ev.RecordingListener)
        # Check there must be at least two events.
        # Because there must be a START and END for the compilation of foo()
        self.assertGreaterEqual(len(rec.buffer), 2)

    def test_compiler_lock_event(self):
        @njit
        def foo(x):
            return x + x

        foo(1)
        md = foo.get_metadata(foo.signatures[0])
        lock_duration = md['timers']['compiler_lock']
        self.assertIsInstance(lock_duration, float)
        self.assertGreater(lock_duration, 0)

    def test_llvm_lock_event(self):
        @njit
        def foo(x):
            return x + x

        foo(1)
        md = foo.get_metadata(foo.signatures[0])
        lock_duration = md['timers']['llvm_lock']
        self.assertIsInstance(lock_duration, float)
        self.assertGreater(lock_duration, 0)

    def test_run_pass_event(self):
        @njit
        def foo(x):
            return x + x

        with ev.install_recorder("numba:run_pass") as recorder:
            foo(2)

        self.assertGreater(len(recorder.buffer), 0)
        for _, event in recorder.buffer:
            # Check that all fields are there
            data = event.data
            self.assertIsInstance(data['name'], str)
            self.assertIsInstance(data['qualname'], str)
            self.assertIsInstance(data['module'], str)
            self.assertIsInstance(data['flags'], str)
            self.assertIsInstance(data['args'], str)
            self.assertIsInstance(data['return_type'], str)

    def test_install_listener(self):
        ut = self

        class MyListener(ev.Listener):
            def on_start(self, event):
                ut.assertEqual(event.status, ev.EventStatus.START)
                ut.assertEqual(event.kind, "numba:compile")
                ut.assertIs(event.data["dispatcher"], foo)
                dispatcher = event.data["dispatcher"]
                ut.assertIs(dispatcher, foo)
                # Check that the compiling signature is NOT in the overloads
                ut.assertNotIn(event.data["args"], dispatcher.overloads)

            def on_end(self, event):
                ut.assertEqual(event.status, ev.EventStatus.END)
                ut.assertEqual(event.kind, "numba:compile")
                dispatcher = event.data["dispatcher"]
                ut.assertIs(dispatcher, foo)
                # Check that the compiling signature is in the overloads
                ut.assertIn(event.data["args"], dispatcher.overloads)

        @njit
        def foo(x):
            return x

        listener = MyListener()
        with ev.install_listener("numba:compile", listener) as yielded:
            foo(1)

        # Check that the yielded value is the same listener
        self.assertIs(listener, yielded)

    def test_global_register(self):
        ut = self

        class MyListener(ev.Listener):
            def on_start(self, event):
                ut.assertEqual(event.status, ev.EventStatus.START)
                ut.assertEqual(event.kind, "numba:compile")
                # Check it is the same dispatcher
                dispatcher = event.data["dispatcher"]
                ut.assertIs(dispatcher, foo)
                # Check that the compiling signature is NOT in the overloads
                ut.assertNotIn(event.data["args"], dispatcher.overloads)

            def on_end(self, event):
                ut.assertEqual(event.status, ev.EventStatus.END)
                ut.assertEqual(event.kind, "numba:compile")
                # Check it is the same dispatcher
                dispatcher = event.data["dispatcher"]
                ut.assertIs(dispatcher, foo)
                # Check that the compiling signature is in the overloads
                ut.assertIn(event.data["args"], dispatcher.overloads)

        @njit
        def foo(x):
            return x

        listener = MyListener()
        ev.register("numba:compile", listener)
        foo(1)
        ev.unregister("numba:compile", listener)

    def test_lifted_dispatcher(self):
        @jit(forceobj=True)
        def foo():
            object()   # to trigger loop-lifting
            c = 0
            for i in range(10):
                c += i
            return c

        with ev.install_recorder("numba:compile") as rec:
            foo()

        # Check that there are 4 events.
        # Two for `foo()` and two for the lifted loop.
        self.assertGreaterEqual(len(rec.buffer), 4)

        cres = foo.overloads[foo.signatures[0]]
        [ldisp] = cres.lifted

        lifted_cres = ldisp.overloads[ldisp.signatures[0]]
        self.assertIsInstance(
            lifted_cres.metadata["timers"]["compiler_lock"],
            float,
        )
        self.assertIsInstance(
            lifted_cres.metadata["timers"]["llvm_lock"],
            float,
        )

    def test_timing_properties(self):
        a = tuple(string.ascii_lowercase)

        @njit
        def bar(x):
            acc = 0
            for i in literal_unroll(a):
                if i in {'1': x}:
                    acc += 1
                else:
                    acc += np.sqrt(x[0, 0])
            return np.sin(x), acc

        @njit
        def foo(x):
            return bar(np.zeros((x, x)))

        with override_config('LLVM_PASS_TIMINGS', True):
            foo(1)

        def get_timers(fn, prop):
            md = fn.get_metadata(fn.signatures[0])
            return md[prop]

        foo_timers = get_timers(foo, 'timers')
        bar_timers = get_timers(bar, 'timers')
        foo_llvm_timer = get_timers(foo, 'llvm_pass_timings')
        bar_llvm_timer = get_timers(bar, 'llvm_pass_timings')

        # Check: time spent in bar() must be longer than in foo()
        self.assertLess(bar_timers['llvm_lock'],
                        foo_timers['llvm_lock'])
        self.assertLess(bar_timers['compiler_lock'],
                        foo_timers['compiler_lock'])

        # Check: time spent in LLVM itself must be less than in the LLVM lock
        self.assertLess(foo_llvm_timer.get_total_time(),
                        foo_timers['llvm_lock'])
        self.assertLess(bar_llvm_timer.get_total_time(),
                        bar_timers['llvm_lock'])

        # Check: time spent in LLVM lock must be less than in compiler
        self.assertLess(foo_timers['llvm_lock'],
                        foo_timers['compiler_lock'])
        self.assertLess(bar_timers['llvm_lock'],
                        bar_timers['compiler_lock'])


if __name__ == "__main__":
    unittest.main()
