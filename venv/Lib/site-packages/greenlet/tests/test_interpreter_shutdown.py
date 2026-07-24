# -*- coding: utf-8 -*-
"""
Tests for greenlet behavior during interpreter shutdown (Py_FinalizeEx).

These tests are organized into four groups:

  A. Core safety (smoke): no crashes with active greenlets at shutdown.
  B. Cleanup semantics: GreenletExit / finally still works during
     normal thread exit (the standard production path).
  C. Atexit "still works" tests: getcurrent() / greenlet construction
     during atexit handlers registered AFTER greenlet import (i.e.
     BEFORE greenlet's cleanup handler in LIFO order) must return
     valid objects — verifies the guards don't over-block.
  D. TDD-certified regression tests: getcurrent() must return None
     when called AFTER greenlet's cleanup (GC finalization phase
     or late atexit phase).  These tests fail on greenlet 3.3.2
     and pass with the fix across Python 3.10-3.14.
"""
import sys
import subprocess
import unittest
import textwrap

from greenlet.tests import TestCase


class TestInterpreterShutdown(TestCase): # pylint:disable=too-many-public-methods

    def _run_shutdown_script(self, script_body):
        """
        Run a Python script in a subprocess that exercises greenlet
        during interpreter shutdown. Returns (returncode, stdout, stderr).
        """
        full_script = textwrap.dedent(script_body)
        result = subprocess.run(
            [sys.executable, '-c', full_script],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        return result.returncode, result.stdout, result.stderr

    # -----------------------------------------------------------------
    # Group A: Core safety — no crashes with active greenlets at exit
    # -----------------------------------------------------------------

    def test_active_greenlet_at_shutdown_no_crash(self):

        # An active (suspended) greenlet that is deallocated during
        # interpreter shutdown should not crash the process.

        # Before the fix, this would SIGSEGV on Python < 3.11 because
        # _green_dealloc_kill_started_non_main_greenlet tried to call
        # g_switch() during Py_FinalizeEx.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet

            def worker():
                greenlet.getcurrent().parent.switch("from worker")
                return "done"

            g = greenlet.greenlet(worker)
            result = g.switch()
            assert result == "from worker", result
            print("OK: exiting with active greenlet")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: exiting with active greenlet", stdout)

    def test_multiple_active_greenlets_at_shutdown(self):
        # Multiple suspended greenlets at shutdown should all be cleaned
        # up without crashing.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet

            def worker(name):
                greenlet.getcurrent().parent.switch(f"hello from {name}")
                return "done"

            greenlets = []
            for i in range(10):
                g = greenlet.greenlet(worker)
                result = g.switch(f"g{i}")
                greenlets.append(g)

            print(f"OK: {len(greenlets)} active greenlets at shutdown")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 10 active greenlets at shutdown", stdout)

    def test_nested_greenlets_at_shutdown(self):
        # Nested (chained parent) greenlets at shutdown should not crash.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet

            def inner():
                greenlet.getcurrent().parent.switch("inner done")

            def outer():
                g_inner = greenlet.greenlet(inner)
                g_inner.switch()
                greenlet.getcurrent().parent.switch("outer done")

            g = greenlet.greenlet(outer)
            result = g.switch()
            assert result == "outer done", result
            print("OK: nested greenlets at shutdown")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: nested greenlets at shutdown", stdout)

    def test_threaded_greenlets_at_shutdown(self):
        # Greenlets in worker threads that are still referenced at
        # shutdown should not crash.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet
            import threading

            results = []

            def thread_worker():
                def greenlet_func():
                    greenlet.getcurrent().parent.switch("from thread greenlet")
                    return "done"

                g = greenlet.greenlet(greenlet_func)
                val = g.switch()
                results.append((g, val))

            threads = []
            for _ in range(3):
                t = threading.Thread(target=thread_worker)
                t.start()
                threads.append(t)

            for t in threads:
                t.join()

            print(f"OK: {len(results)} threaded greenlets at shutdown")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 3 threaded greenlets at shutdown", stdout)

    # -----------------------------------------------------------------
    # Group B: Cleanup semantics — thread exit
    # -----------------------------------------------------------------
    #
    # These tests verify that GreenletExit / try-finally still work
    # correctly during normal thread exit (the standard production
    # path, e.g. uWSGI worker threads finishing a request).  This is
    # NOT interpreter shutdown; the guards do not fire here.

    def test_greenlet_cleanup_during_thread_exit(self):
        # When a thread exits normally while holding active greenlets,
        # GreenletExit IS thrown and cleanup code runs.  This is the
        # standard cleanup path used in production (e.g. uWSGI worker
        # threads finishing a request).

        rc, stdout, stderr = self._run_shutdown_script("""\
            import os
            import threading
            import greenlet

            _write = os.write

            def thread_func():
                def worker(_w=_write,
                           _GreenletExit=greenlet.GreenletExit):
                    try:
                        greenlet.getcurrent().parent.switch("suspended")
                    except _GreenletExit:
                        _w(1, b"CLEANUP: GreenletExit caught\\n")
                        raise

                g = greenlet.greenlet(worker)
                g.switch()
                # Thread exits with active greenlet -> thread-state
                # cleanup triggers GreenletExit

            t = threading.Thread(target=thread_func)
            t.start()
            t.join()
            print("OK: thread cleanup done")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: thread cleanup done", stdout)
        self.assertIn("CLEANUP: GreenletExit caught", stdout)

    def test_finally_block_during_thread_exit(self):
        # try/finally blocks in active greenlets run correctly when the
        # owning thread exits.
        rc, stdout, stderr = self._run_shutdown_script("""\
            import os
            import threading
            import greenlet

            _write = os.write

            def thread_func():
                def worker(_w=_write):
                    try:
                        greenlet.getcurrent().parent.switch("suspended")
                    finally:
                        _w(1, b"FINALLY: cleanup executed\\n")

                g = greenlet.greenlet(worker)
                g.switch()

            t = threading.Thread(target=thread_func)
            t.start()
            t.join()
            print("OK: thread cleanup done")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: thread cleanup done", stdout)
        self.assertIn("FINALLY: cleanup executed", stdout)

    def test_many_greenlets_with_cleanup_at_shutdown(self):
        # Stress test: many active greenlets with cleanup code at shutdown.
        # Ensures no crashes regardless of deallocation order.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import sys
            import greenlet

            cleanup_count = 0

            def worker(idx):
                global cleanup_count
                try:
                    greenlet.getcurrent().parent.switch(f"ready-{idx}")
                except greenlet.GreenletExit:
                    cleanup_count += 1
                    raise

            greenlets = []
            for i in range(50):
                g = greenlet.greenlet(worker)
                result = g.switch(i)
                greenlets.append(g)

            print(f"OK: {len(greenlets)} greenlets about to shut down")
            # Note: we can't easily print cleanup_count during shutdown
            # since it happens after the main module's code runs.
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 50 greenlets about to shut down", stdout)

    def test_deeply_nested_greenlets_at_shutdown(self):
        # Deeply nested greenlet parent chains at shutdown.
        # Tests that the deallocation order doesn't cause issues.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet

            def level(depth, max_depth):
                if depth < max_depth:
                    g = greenlet.greenlet(level)
                    g.switch(depth + 1, max_depth)
                greenlet.getcurrent().parent.switch(f"depth-{depth}")

            g = greenlet.greenlet(level)
            result = g.switch(0, 10)
            print(f"OK: nested to depth 10, got {result}")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: nested to depth 10", stdout)

    def test_greenlet_with_traceback_at_shutdown(self):
        # A greenlet that has an active exception context when it's
        # suspended should not crash during shutdown cleanup.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet

            def worker():
                try:
                    raise ValueError("test error")
                except ValueError:
                    # Suspend while an exception is active on the stack
                    greenlet.getcurrent().parent.switch("suspended with exc")
                return "done"

            g = greenlet.greenlet(worker)
            result = g.switch()
            assert result == "suspended with exc"
            print("OK: greenlet with active exception at shutdown")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: greenlet with active exception at shutdown", stdout)

    # -----------------------------------------------------------------
    # Group C: getcurrent() / construction / gettrace() / settrace()
    # during atexit — registered AFTER greenlet import
    # -----------------------------------------------------------------

    def test_getcurrent_during_atexit_no_crash(self):
        # getcurrent() in an atexit handler registered AFTER greenlet
        # import must return a valid greenlet (not None), because LIFO
        # ordering means this handler runs BEFORE greenlet's cleanup.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet

            def call_getcurrent_at_exit():
                try:
                    g = greenlet.getcurrent()
                    if g is not None and type(g).__name__ == 'greenlet':
                        print(f"OK: getcurrent returned valid greenlet")
                    elif g is None:
                        print("FAIL: getcurrent returned None (over-blocked)")
                    else:
                        print(f"FAIL: unexpected {g!r}")
                except Exception as e:
                    print(f"OK: getcurrent raised {type(e).__name__}: {e}")

            atexit.register(call_getcurrent_at_exit)
            print("OK: atexit registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: atexit registered", stdout)
        self.assertIn("OK: getcurrent returned valid greenlet", stdout,
                      "getcurrent() should return a valid greenlet when called "
                      "before greenlet's cleanup handler (LIFO ordering)")

    def test_gettrace_during_atexit_no_crash(self):
        # Calling greenlet.gettrace() during atexit must not crash.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet

            def check_at_exit():
                try:
                    result = greenlet.gettrace()
                    print(f"OK: gettrace returned {result!r}")
                except Exception as e:
                    print(f"OK: gettrace raised {type(e).__name__}: {e}")

            atexit.register(check_at_exit)
            print("OK: registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: registered", stdout)

    def test_settrace_during_atexit_no_crash(self):
        # Calling greenlet.settrace() during atexit must not crash.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet

            def check_at_exit():
                try:
                    greenlet.settrace(lambda *args: None)
                    print("OK: settrace succeeded")
                except Exception as e:
                    print(f"OK: settrace raised {type(e).__name__}: {e}")

            atexit.register(check_at_exit)
            print("OK: registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: registered", stdout)

    def test_getcurrent_with_active_greenlets_during_atexit(self):
        # getcurrent() during atexit (registered after import) with active
        # greenlets must still return a valid greenlet, since LIFO means
        # this runs before greenlet's cleanup.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet

            def worker():
                greenlet.getcurrent().parent.switch("ready")

            greenlets = []
            for i in range(5):
                g = greenlet.greenlet(worker)
                result = g.switch()
                greenlets.append(g)

            def check_at_exit():
                try:
                    g = greenlet.getcurrent()
                    if g is not None and type(g).__name__ == 'greenlet':
                        print(f"OK: getcurrent returned valid greenlet")
                    elif g is None:
                        print("FAIL: getcurrent returned None (over-blocked)")
                    else:
                        print(f"FAIL: unexpected {g!r}")
                except Exception as e:
                    print(f"OK: getcurrent raised {type(e).__name__}: {e}")

            atexit.register(check_at_exit)
            print(f"OK: {len(greenlets)} active greenlets, atexit registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 5 active greenlets, atexit registered", stdout)
        self.assertIn("OK: getcurrent returned valid greenlet", stdout,
                      "getcurrent() should return a valid greenlet when called "
                      "before greenlet's cleanup handler (LIFO ordering)")

    def test_greenlet_construction_during_atexit_no_crash(self):
        # Constructing a new greenlet during atexit (registered after
        # import) must succeed, since this runs before greenlet's cleanup.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet

            def create_greenlets_at_exit():
                try:
                    def noop():
                        pass
                    g = greenlet.greenlet(noop)
                    if g is not None:
                        print(f"OK: created greenlet successfully")
                    else:
                        print("FAIL: greenlet() returned None")
                except Exception as e:
                    print(f"OK: construction raised {type(e).__name__}: {e}")

            atexit.register(create_greenlets_at_exit)
            print("OK: atexit registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: atexit registered", stdout)
        self.assertIn("OK: created greenlet successfully", stdout)

    def test_greenlet_construction_with_active_greenlets_during_atexit(self):
        # Constructing new greenlets during atexit when other active
        # greenlets already exist (maximizes the chance of a non-empty
        # deleteme list).

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet

            def worker():
                greenlet.getcurrent().parent.switch("ready")

            greenlets = []
            for i in range(10):
                g = greenlet.greenlet(worker)
                g.switch()
                greenlets.append(g)

            def create_at_exit():
                try:
                    new_greenlets = []
                    for i in range(5):
                        g = greenlet.greenlet(lambda: None)
                        new_greenlets.append(g)
                    print(f"OK: created {len(new_greenlets)} greenlets at exit")
                except Exception as e:
                    print(f"OK: raised {type(e).__name__}: {e}")

            atexit.register(create_at_exit)
            print(f"OK: {len(greenlets)} active greenlets, atexit registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 10 active greenlets, atexit registered", stdout)

    def test_greenlet_construction_with_cross_thread_deleteme_during_atexit(self):
        # Create greenlets in a worker thread, transfer them to the main
        # thread, then drop them — populating the deleteme list. Then
        # construct a new greenlet during atexit. On Python < 3.11
        # clear_deleteme_list() could previously crash if the
        # PythonAllocator vector copy failed during early Py_FinalizeEx;
        # using std::swap eliminates that allocation.
        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import greenlet
            import threading

            cross_thread_refs = []

            def thread_worker():
                # Create greenlets in this thread
                def gl_body():
                    greenlet.getcurrent().parent.switch("ready")
                for _ in range(20):
                    g = greenlet.greenlet(gl_body)
                    g.switch()
                    cross_thread_refs.append(g)

            t = threading.Thread(target=thread_worker)
            t.start()
            t.join()

            # Dropping these references in the main thread
            # causes them to be added to the main thread's
            # deleteme list (deferred cross-thread dealloc).
            cross_thread_refs.clear()

            def create_at_exit():
                try:
                    g = greenlet.greenlet(lambda: None)
                    print(f"OK: created greenlet at exit {g!r}")
                except Exception as e:
                    print(f"OK: raised {type(e).__name__}: {e}")

            atexit.register(create_at_exit)
            print("OK: cross-thread setup done, atexit registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: cross-thread setup done, atexit registered", stdout)


    # -----------------------------------------------------------------
    # Group D.1: TDD-certified — getcurrent() during GC finalization
    #
    # These tests use gc.disable() + reference cycles to force __del__
    # to run during Py_FinalizeEx's GC pass, where Py_IsFinalizing()
    # is True.  Without the fix, getcurrent() returns a live greenlet
    # (unguarded); with the fix, it returns None.
    #
    # TDD verification (greenlet 3.3.2 = RED, patched = GREEN):
    #   Python 3.10: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.11: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.12: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.13: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.14: RED (UNGUARDED) → GREEN (GUARDED)
    # -----------------------------------------------------------------

    def test_getcurrent_returns_none_during_gc_finalization(self):
        # greenlet.getcurrent() must raise an exception when called from a
        # __del__ method during Py_FinalizeEx's GC collection pass.

        # On Python >= 3.11, _Py_IsFinalizing() is True during this
        # phase.  Without the Py_IsFinalizing() guard in mod_getcurrent,
        # this would return a greenlet — the same unguarded code path
        # that leads to SIGSEGV in production (uWSGI worker recycling).
        rc, stdout, stderr = self._run_shutdown_script("""\
            import gc
            import os
            import greenlet

            gc.disable()

            class CleanupChecker:
                def __del__(self):
                    try:
                        try:
                            greenlet.getcurrent()
                        except RuntimeError:
                            os.write(1, b"GUARDED: getcurrent=None\\n")
                        else:
                            os.write(1, b"UNGUARDED: getcurrent="
                                     + type(cur).__name__.encode() + b"\\n")
                    except Exception as e:
                        os.write(1, b"EXCEPTION: " + str(e).encode() + b"\\n")

            # Reference cycle: only collected during Py_FinalizeEx GC pass
            a = CleanupChecker()
            b = {"ref": a}
            a._cycle = b
            del a, b

            print("OK: deferred cycle created")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: deferred cycle created", stdout)
        self.assertIn("GUARDED: getcurrent=None", stdout)

    def test_getcurrent_returns_none_during_gc_finalization_with_active_greenlets(self):
        # Same as above but with active greenlets at shutdown, which
        # increases the amount of C++ destructor work during finalization.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import gc
            import os
            import greenlet

            gc.disable()

            class CleanupChecker:
                def __del__(self):
                    try:
                        try:
                            greenlet.getcurrent()
                        except RuntimeError:
                            os.write(1, b"GUARDED: getcurrent=None\\n")
                        else:
                            os.write(1, b"UNGUARDED: getcurrent="
                                     + type(cur).__name__.encode() + b"\\n")
                    except Exception as e:
                        os.write(1, b"EXCEPTION: " + str(e).encode() + b"\\n")

            # Create active greenlets
            greenlets = []
            for _ in range(10):
                def worker():
                    greenlet.getcurrent().parent.switch("suspended")
                g = greenlet.greenlet(worker)
                g.switch()
                greenlets.append(g)

            # Reference cycle deferred to Py_FinalizeEx
            a = CleanupChecker()
            b = {"ref": a}
            a._cycle = b
            del a, b

            print(f"OK: {len(greenlets)} active greenlets, cycle deferred")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 10 active greenlets, cycle deferred", stdout)
        self.assertIn("GUARDED: getcurrent=None", stdout)

    def test_getcurrent_returns_none_during_gc_finalization_cross_thread(self):
        # Combines cross-thread greenlet deallocation (deleteme list)
        # with the GC finalization check.  This simulates the production
        # scenario where uWSGI worker threads create greenlets that are
        # transferred to the main thread, then cleaned up during
        # Py_FinalizeEx.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import gc
            import os
            import threading
            import greenlet

            gc.disable()

            class CleanupChecker:
                def __del__(self):
                    try:
                        try:
                            greenlet.getcurrent()
                        except RuntimeError:
                            os.write(1, b"GUARDED: getcurrent=None\\n")
                        else:
                            os.write(1, b"UNGUARDED: getcurrent="
                                     + type(cur).__name__.encode() + b"\\n")
                    except Exception as e:
                        os.write(1, b"EXCEPTION: " + str(e).encode() + b"\\n")

            # Create cross-thread greenlet references
            cross_refs = []
            def thread_fn():
                for _ in range(20):
                    def body():
                        greenlet.getcurrent().parent.switch("x")
                    g = greenlet.greenlet(body)
                    g.switch()
                    cross_refs.append(g)
            t = threading.Thread(target=thread_fn)
            t.start()
            t.join()
            cross_refs.clear()

            # Reference cycle deferred to Py_FinalizeEx
            a = CleanupChecker()
            b = {"ref": a}
            a._cycle = b
            del a, b

            print("OK: cross-thread cleanup + cycle deferred")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: cross-thread cleanup + cycle deferred", stdout)
        self.assertIn("GUARDED: getcurrent=None", stdout)


    # -----------------------------------------------------------------
    # Group D.2: TDD-certified — getcurrent() during atexit phase
    #
    # These tests register the checker BEFORE importing greenlet.
    # Python's atexit is LIFO, so greenlet's handler (registered at
    # import) runs FIRST and sets g_greenlet_shutting_down=1; then
    # the checker runs SECOND and observes getcurrent() → None.
    #
    # This covers the atexit phase where _Py_IsFinalizing() is still
    # False on ALL Python versions — the exact window that causes
    # SIGSEGV in production (uWSGI worker recycling → Py_FinalizeEx).
    #
    # TDD verification (greenlet 3.3.2 = RED, patched = GREEN):
    #   Python 3.10: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.11: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.12: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.13: RED (UNGUARDED) → GREEN (GUARDED)
    #   Python 3.14: RED (UNGUARDED) → GREEN (GUARDED)
    # -----------------------------------------------------------------

    def test_getcurrent_returns_none_during_atexit_phase(self):
        # greenlet.getcurrent() must NOT return None when called from an
        # atexit handler that runs AFTER greenlet's own atexit handler.

        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import os

            def late_checker():
                try:
                    import greenlet
                    cur = greenlet.getcurrent()
                    if cur is None:
                        os.write(1, b"GUARDED: getcurrent=None\\n")
                    else:
                        os.write(1, b"UNGUARDED: getcurrent="
                                 + type(cur).__name__.encode() + b"\\n")
                except Exception as e:
                    os.write(1, b"EXCEPTION: " + str(e).encode() + b"\\n")

            # Register BEFORE importing greenlet.  LIFO order:
            # greenlet's handler (registered at import) runs FIRST,
            # late_checker runs SECOND — seeing the flag already set.
            atexit.register(late_checker)

            import greenlet
            print("OK: atexit registered before greenlet import")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: atexit registered before greenlet import", stdout)
        self.assertIn("UNGUARDED", stdout)


    def test_getcurrent_returns_none_during_atexit_phase_with_active_greenlets(self):
        # Same as above but with active greenlets
        rc, stdout, stderr = self._run_shutdown_script("""\
            import atexit
            import os

            def late_checker():
                try:
                    import greenlet
                    cur = greenlet.getcurrent()
                    if cur is None:
                        os.write(1, b"GUARDED: getcurrent=None\\n")
                    else:
                        os.write(1, b"UNGUARDED: getcurrent="
                                 + type(cur).__name__.encode() + b"\\n")
                except Exception as e:
                    os.write(1, b"EXCEPTION: " + str(e).encode() + b"\\n")

            atexit.register(late_checker)

            import greenlet

            greenlets = []
            for _ in range(10):
                def worker():
                    greenlet.getcurrent().parent.switch("parked")
                g = greenlet.greenlet(worker)
                g.switch()
                greenlets.append(g)

            print(f"OK: {len(greenlets)} active greenlets, atexit registered")
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn("OK: 10 active greenlets, atexit registered", stdout)
        self.assertIn("UNGUARDED", stdout)

    def test_api_getcurrent_no_system_error_at_module_gc_time(self):
        # If we use the C API directly to return a greenlet AFTER
        # atexit threads have been run, we don't crash, we get a
        # specific error. We arrange for this by putting a __del__ on
        # an object that lives in greenlet's own (extension module)
        # dict; this is cleaned out sometime during the module cleanup
        # steps.
        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet
            from greenlet.tests import _test_extension

            class WithDel:
                # must cache the method we want, because by the time we
                # run, module globals may have been cleaned up.
                def __del__(self, gc=_test_extension.getcurrent_api):
                    print('Destructor running')
                    gc() # Should print an unraisable RuntimeException

            greenlet._greenlet.with_del = WithDel()
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn('Destructor running', stdout)
        self.assertIn('RuntimeError: greenlet is being finalized', stderr)


    def test_switch_no_error_at_module_gc_time(self):
        # Switching to a greenlet we've captured during
        # module tear down doesn't cause a crash
        rc, stdout, stderr = self._run_shutdown_script("""\
            import greenlet
            from greenlet.tests import _test_extension

            gs = []
            # must cache the objects we want, because by the time we
            # run, module globals may have been cleaned up.
            def do_it(gs=gs):
                print('current', gs)
                gs[0].parent.switch(1)


            gs.append(greenlet.greenlet(do_it))
            gs.append(greenlet.greenlet(do_it))
            gs[1].switch()

            class WithDel:
                def __del__(self, gs=gs):
                    print('Destructor running')
                    r = gs[0].switch()
                    print('Result', r)

            greenlet._greenlet.with_del = WithDel()
        """)
        self.assertEqual(rc, 0, f"Process crashed (rc={rc}):\n{stdout}{stderr}")
        self.assertIn('Destructor running', stdout)
        self.assertIn('Result 1', stdout)


if __name__ == '__main__':
    unittest.main()
