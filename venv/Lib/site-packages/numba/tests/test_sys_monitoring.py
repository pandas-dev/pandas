import dis
import queue
import sys
import threading
import traceback
import unittest
from collections import Counter
from unittest.mock import Mock, call
from numba.tests.support import TestCase
from numba import jit, objmode
from numba.core.utils import PYVERSION
from numba.core.serialize import _numba_unpickle


def _enable_sysmon(disp):
    """Decorator to enable sys.monitoring on the dispatcher"""
    disp._enable_sysmon = True
    return disp


def generate_usecase():
    @_enable_sysmon
    @jit('int64(int64)',)
    def foo(x):
        return x + 1

    def call_foo(x):
        return 2 * foo(x + 5)

    return foo, call_foo


if PYVERSION in ((3, 12), (3, 13), (3, 14)):
    PY_START = sys.monitoring.events.PY_START
    PY_RETURN = sys.monitoring.events.PY_RETURN
    RAISE = sys.monitoring.events.RAISE
    PY_UNWIND = sys.monitoring.events.PY_UNWIND
    STOP_ITERATION = sys.monitoring.events.STOP_ITERATION
    NO_EVENTS = sys.monitoring.events.NO_EVENTS


TOOL2MONITORTYPE = {0 : "Debugger",
                    1 : "Coverage",
                    2 : "Profiler",
                    5 : "Optimizer"}


@unittest.skipUnless(PYVERSION >= (3, 12), "needs Python 3.12+")
class TestMonitoring(TestCase):
    # Tests the interaction of the Numba dispatcher with `sys.monitoring`.
    #
    # Note that it looks like a lot of these try..finally type patterns could
    # be written using a context manager, this is true, but it is not written
    # like that deliberately as a context manager adds implementation details
    # onto the stack which makes it harder to debug tests.

    def setUp(self):
        # First... check if there's other monitoring stuff registered (e.g. test
        # is running under cProfile or coverage), skip if so.
        monitor_kinds = []
        for i in range(6): # there are 5 tool IDs
            if sys.monitoring.get_tool(i) is not None:
                monitor_kinds.append(TOOL2MONITORTYPE[i])

        if monitor_kinds:
            msg = ("Cannot run monitoring tests when other monitors are "
                   "active, found monitor(s) of type: "
                   f"{', '.join(monitor_kinds)}")
            self.skipTest(msg)

        # set up some standard functions and answers for use throughout
        self.foo, self.call_foo = generate_usecase()
        self.arg = 10
        self.foo_result = self.arg + 5 + 1
        self.call_foo_result = 2 * self.foo_result
        # pretend to be a profiler in the majority of these unit tests
        self.tool_id = sys.monitoring.PROFILER_ID

    def gather_mock_calls_multithreads(self, mockcalls):
        # Gather mock-calls for the `self.foo` and `self.call_foo`
        matched = Counter()
        target_codeobjs = {self.call_foo.__code__, self.foo.__code__}
        for cb_args in mockcalls._mock_call_args_list:
            (codeobj, *args) = cb_args.args
            if codeobj in target_codeobjs:
                matched[(codeobj, *args)] += 1
        return matched

    def check_py_start_calls_multithreads(self, allcalls):
        # Checks that PY_START calls were correctly captured for a
        # `self.call_foo(self.arg)` call in multithreads
        matched = self.gather_mock_calls_multithreads(allcalls[PY_START])
        self.assertEqual(len(matched), 2) # two types of call

        # Find the resume op, this is where the code for `call_foo` "starts"
        inst = [x for x in dis.get_instructions(self.call_foo)
                if x.opname == "RESUME"]
        offset = inst[0].offset
        self.assertEqual(matched[self.call_foo.__code__, offset], 2)
        self.assertEqual(matched[self.foo.__code__, 0], 2)
        self.assertEqual(matched.total(), 4)

    def check_py_start_calls(self, allcalls):
        # Checks that PY_START calls were correctly captured for a
        # `self.call_foo(self.arg)` call.
        mockcalls = allcalls[PY_START]
        self.assertEqual(mockcalls.call_count, 2)
        # Find the resume op, this is where the code for `call_foo` "starts"
        inst = [x for x in dis.get_instructions(self.call_foo)
                if x.opname == "RESUME"]
        offset = inst[0].offset
        # Numba always reports the start location as offset 0.
        calls = (call(self.call_foo.__code__, offset),
                 call(self.foo.__code__, 0))
        mockcalls.assert_has_calls(calls)

    def check_py_return_calls_multithreads(self, allcalls):
        # Checks that PY_RETURN calls were correctly captured for a
        # `self.call_foo(self.arg)` call.
        matched = self.gather_mock_calls_multithreads(allcalls[PY_RETURN])
        offset = [x for x in dis.get_instructions(self.call_foo)][-1].offset
        self.assertEqual(matched[self.foo.__code__, 0, self.foo_result], 2)
        self.assertEqual(
            matched[self.call_foo.__code__, offset, self.call_foo_result], 2
        )
        self.assertEqual(matched.total(), 4)

    def check_py_return_calls(self, allcalls):
        # Checks that PY_RETURN calls were correctly captured for a
        # `self.call_foo(self.arg)` call.
        mockcalls = allcalls[PY_RETURN]
        self.assertEqual(mockcalls.call_count, 2)
        # These are in the order the returns were encountered. Return from `foo`
        # occurred first, followed by return from `call_foo`.
        # NOTE: it is a known issue that Numba reports the PY_RETURN event as
        # occurring at offset 0. At present there's no information about the
        # location that the return occurred propagating from the machine code
        # back to the dispatcher (where the monitoring events are handled).
        offset = [x for x in dis.get_instructions(self.call_foo)][-1].offset
        calls = [call(self.foo.__code__, 0, self.foo_result),
                 call(self.call_foo.__code__, offset, self.call_foo_result)]
        mockcalls.assert_has_calls(calls)

    def run_with_events(self, function, args, events, tool_id=None,
                        barrier=None):
        # Runs function with args with monitoring set for events on `tool_id`
        # (if present, else just uses the default of "PROFILER_ID") returns a
        # dictionary event->callback.
        try:
            if tool_id is None:
                _tool_id = self.tool_id
            else:
                _tool_id = tool_id
            sys.monitoring.use_tool_id(_tool_id, "custom_monitor")
            callbacks = {}
            event_bitmask = 0
            for event in events:
                callback = Mock()
                sys.monitoring.register_callback(_tool_id, event, callback)
                callbacks[event] = callback
                event_bitmask |= event
            # only start monitoring once callbacks are registered
            sys.monitoring.set_events(_tool_id, event_bitmask)
            if barrier is not None:
                # Wait for all threads to have enabled events.
                barrier()
            function(*args)
        finally:
            # clean up state
            if barrier is not None:
                # Wait for all threads to finish `function()`
                # This makes sure all threads have a chance to see the events
                # from other threads.
                barrier()
            sys.monitoring.set_events(_tool_id, NO_EVENTS)
            for event in events:
                sys.monitoring.register_callback(_tool_id, event, None)
            sys.monitoring.free_tool_id(_tool_id)
        return callbacks

    def test_start_event(self):
        # test event PY_START
        cb = self.run_with_events(self.call_foo, (self.arg,), (PY_START,))
        # Check...
        self.assertEqual(len(cb), 1)
        self.check_py_start_calls(cb)

    def test_return_event(self):
        # test event PY_RETURN
        cb = self.run_with_events(self.call_foo, (self.arg,), (PY_RETURN,))
        # Check...
        self.assertEqual(len(cb), 1)
        self.check_py_return_calls(cb)

    def test_call_event_chain(self):
        # test event PY_START and PY_RETURN monitored at the same time
        cb = self.run_with_events(self.call_foo, (self.arg,),
                                  (PY_START, PY_RETURN))
        # Check...
        self.assertEqual(len(cb), 2)
        self.check_py_return_calls(cb)
        self.check_py_start_calls(cb)

    # --------------------------------------------------------------------------
    # NOTE: About the next two tests...
    # Numba doesn't support "local event" level monitoring, it's implemented
    # in CPython via adjusting the code object bytecode to use
    # "instrumented" opcodes. When the interpreter encounters an
    # instrumented opcode it triggers the event handling pathways. As Numba
    # doesn't interpret the bytecode instruction-at-a-time there's not
    # really any way to support this. Two things to check...
    # 1. The an instrumented code object doesn't trigger events in
    #    the dispatcher.
    # 2. That Numba can compile instrumented functions (it should be able
    #    to without any problem as the instrumented bytecode should not
    #    leak into `.co_code`.).

    def test_instrumented_code_does_not_trigger_numba_events(self):
        # 1. from above.
        @jit('int64(int64)',)
        def foo(x):
            return x + 3

        try:
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            callbacks = {}
            event_bitmask = 0
            events = (PY_START, PY_RETURN)
            for event in events:
                callback = Mock()
                sys.monitoring.register_callback(tool_id, event, callback)
                callbacks[event] = callback
                event_bitmask |= event

            sys.monitoring.set_local_events(tool_id, foo.__code__,
                                            event_bitmask)
            result = foo(self.arg)
        finally:
            for event in events:
                sys.monitoring.register_callback(tool_id, event, None)
            sys.monitoring.set_local_events(tool_id, foo.__code__, 0)
            sys.monitoring.free_tool_id(tool_id)

        # check
        self.assertEqual(result, foo.py_func(self.arg))
        self.assertEqual(len(callbacks), 2)
        callbacks[PY_START].assert_not_called()
        callbacks[PY_RETURN].assert_not_called()

    def test_instrumented_code_can_be_compiled(self):
        # 2. from above.

        def foo(x):
            return x + 1

        try:
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.set_local_events(tool_id, foo.__code__, PY_START)
            sys.monitoring.register_callback(tool_id, PY_START, Mock())
            # test compile
            result = jit(foo)(self.arg)
            self.assertEqual(result, foo(self.arg))
        finally:
            sys.monitoring.register_callback(tool_id, PY_START, None)
            sys.monitoring.set_local_events(tool_id, foo.__code__, 0)
            sys.monitoring.free_tool_id(tool_id)

    def test_unhandled_events_are_ignored(self):
        # Check an unhandled event e.g. PY_YIELD isn't reported.
        def generate(dec):
            @dec('void()')
            def producer():
                yield 10

            @dec('int64()')
            def consumer():
                p = producer()
                return next(p)

            return consumer

        event = sys.monitoring.events.PY_YIELD
        # check that pure python reports
        wrapper = lambda sig: lambda fn: fn
        py_consumer = generate(wrapper)
        py_cb = self.run_with_events(py_consumer, (),  (event,))
        py_cb[event].assert_called_once()
        # check the numba does not report
        nb_consumer = generate(jit)
        nb_cb = self.run_with_events(nb_consumer, (),  (event,))
        nb_cb[event].assert_not_called()

    def test_event_with_no_callback_runs(self):
        # This checks the situation where an event is being monitored but
        # there's no callback associated with the event. In the dispatcher C
        # code the loop over tools will be entered, but nothing will get called
        # as the "instrument" is missing (NULL).
        try:
            event = PY_START
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.set_events(tool_id, event)
            # NO CALLBACK IS REGISTERED!
            active_events = sys.monitoring.get_events(tool_id)
            self.assertEqual(active_events, event)
            result = self.call_foo(self.arg)
            active_events = sys.monitoring.get_events(tool_id)
            self.assertEqual(active_events, event)
            self.assertEqual(result, self.call_foo_result)
        finally:
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.free_tool_id(tool_id)

    def test_disable_from_callback(self):
        # Event callbacks can disable a _local_ event at a specific location to
        # prevent it triggering in the future by returning
        # `sys.monitoring.DISABLE`. As this only applies to local events, doing
        # this should have absolutely no impact for the global events that Numba
        # supports.

        callback = Mock(return_value=sys.monitoring.DISABLE)

        try:
            event = PY_START
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.set_events(tool_id, event)
            sys.monitoring.register_callback(tool_id, event, callback)
            active_events = sys.monitoring.get_events(tool_id)
            self.assertEqual(active_events, event)
            result = self.call_foo(self.arg)
            active_events = sys.monitoring.get_events(tool_id)
            self.assertEqual(active_events, event)
            self.assertEqual(result, self.call_foo_result)
            callback.assert_called()
        finally:
            # It is necessary to restart events that have been disabled. The
            # "disabled" state of the `PY_START` event for the tool
            # `self.tool_id` "leaks" into subsequent tests. These subsequent
            # tests then end up failing as events that should have been
            # triggered are not triggered due to the state leak! It's not really
            # clear why this happens, if it is part of the design or a side
            # effect of the design, or if this behaviour is simply a bug in
            # CPython itself.
            sys.monitoring.restart_events()
            sys.monitoring.register_callback(tool_id, event, None)
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.free_tool_id(tool_id)

    def test_mutation_from_objmode(self):
        try:
            # Check that it's possible to enable an event (mutate the event
            # state)from an `objmode` block. Monitoring for PY_RETURN is set in
            # objmode once the function starts executing.
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            event = PY_RETURN
            # register the callback... note that the event isn't switched on yet
            callback = Mock()
            sys.monitoring.register_callback(tool_id, event, callback)

            def objmode_enable_event(switch_on_event):
                if switch_on_event:
                    sys.monitoring.set_events(tool_id, event)

            @_enable_sysmon
            @jit('int64(int64)')
            def foo(enable):
                with objmode:
                    objmode_enable_event(enable)
                return enable + 7

            # this should not trigger the return callback
            foo(0)
            callback.assert_not_called()

            # this should trigger the return callback
            foo(1)
            # switch off the event so the callback mock is protected from
            # mutation.
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            # check what happened
            callback.assert_called()
            # 2 calls, 1 is the return from the objmode_enable_event, the other
            # is the return from foo.
            self.assertEqual(callback.call_count, 2)
        finally:
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.register_callback(tool_id, event, None)
            sys.monitoring.free_tool_id(tool_id)

    def test_multiple_tool_id(self):
        # Check that multiple tools will work across different combinations of
        # events that Numba dispatcher supports, namely:
        # (NO_EVENTS, PY_START, PY_RETURN).

        # the use of NO_EVENTS is superfluous, it is to demonstrate usage.
        tool_ids_2_events = {sys.monitoring.DEBUGGER_ID: (NO_EVENTS,),
                             sys.monitoring.COVERAGE_ID: (PY_START,),
                             sys.monitoring.PROFILER_ID: (PY_RETURN,),
                             sys.monitoring.OPTIMIZER_ID:
                                 (PY_START, PY_RETURN,),}

        all_callbacks = {}
        try:
            for tool_id, events in tool_ids_2_events.items():
                sys.monitoring.use_tool_id(tool_id, f"custom_monitor_{tool_id}")
                event_bitmask = 0
                callbacks = {}
                all_callbacks[tool_id] = callbacks
                for event in events:
                    callback = Mock()
                    # Can't set an event for NO_EVENTS!
                    if event != NO_EVENTS:
                        sys.monitoring.register_callback(tool_id, event,
                                                         callback)
                    callbacks[event] = callback
                    event_bitmask |= event
                # only start monitoring once callbacks are registered
            for tool_id in tool_ids_2_events.keys():
                sys.monitoring.set_events(tool_id, event_bitmask)
            self.call_foo(self.arg)
        finally:
            # clean up state
            for tool_id, events in tool_ids_2_events.items():
                for event in events:
                    # Can't remove an event for NO_EVENTS!
                    if event != NO_EVENTS:
                        sys.monitoring.register_callback(tool_id, event, None)
                sys.monitoring.set_events(tool_id, NO_EVENTS)
                sys.monitoring.free_tool_id(tool_id)

        # Now check all_callbacks...

        # check debugger tool slot
        dbg_tool = all_callbacks[sys.monitoring.DEBUGGER_ID]
        self.assertEqual(len(dbg_tool), 1) # one event to capture
        callback = dbg_tool[NO_EVENTS]
        callback.assert_not_called()

        # check coverage tool slot
        cov_tool = all_callbacks[sys.monitoring.COVERAGE_ID]
        self.assertEqual(len(cov_tool), 1) # one event to capture
        self.check_py_start_calls(cov_tool)

        # check profiler tool slot
        prof_tool = all_callbacks[sys.monitoring.PROFILER_ID]
        self.assertEqual(len(prof_tool), 1) # one event to capture
        self.check_py_return_calls(prof_tool)

        # check optimiser tool slot
        opt_tool = all_callbacks[sys.monitoring.OPTIMIZER_ID]
        self.assertEqual(len(opt_tool), 2) # two events to capture
        self.check_py_start_calls(opt_tool)
        self.check_py_return_calls(opt_tool)

    def test_raising_under_monitoring(self):
        # Check that Numba can raise an exception whilst monitoring is running
        # and that 1. `RAISE` is issued 2. `PY_UNWIND` is issued, 3. that
        # `PY_RETURN` is not issued.

        ret_callback = Mock()
        raise_callback = Mock()
        unwind_callback = Mock()

        msg = 'exception raised'

        @_enable_sysmon
        @jit('()')
        def foo():
            raise ValueError(msg)

        store_raised = None
        try:
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.register_callback(tool_id, PY_RETURN, ret_callback)
            sys.monitoring.register_callback(tool_id, RAISE, raise_callback)
            sys.monitoring.register_callback(tool_id, PY_UNWIND,
                                             unwind_callback)
            sys.monitoring.set_events(tool_id, PY_RETURN | RAISE | PY_UNWIND)
            try:
                foo()
            except ValueError as raises:
                store_raised = raises
            # switch off monitoring
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            # check that the ret_callback was called once (by Numba unpickle to
            # fetch the exception info out of the stored bytes).
            ret_callback.assert_called_once()
            # and that elements that are feasible to check about the call are
            # as expected
            the_call = ret_callback.call_args_list[0]
            self.assertEqual(the_call.args[0], _numba_unpickle.__code__)
            self.assertEqual(the_call.args[2][0], ValueError)
            self.assertEqual(the_call.args[2][1][0], msg)

            # check that the RAISE event callback was triggered
            raise_callback.assert_called()
            numba_unpickle_call = raise_callback.call_args_list[0]
            self.assertEqual(numba_unpickle_call.args[0],
                             _numba_unpickle.__code__)
            self.assertIsInstance(numba_unpickle_call.args[2], KeyError)
            foo_call = raise_callback.call_args_list[1]
            self.assertEqual(foo_call.args[0], foo.py_func.__code__)
            self.assertIsInstance(foo_call.args[2], ValueError)
            self.assertIn(msg, str(foo_call.args[2]))

            # check that PY_UNWIND event callback was called
            unwind_callback.assert_called_once()
            unwind_call = unwind_callback.call_args_list[0]
            self.assertEqual(unwind_call.args[0], foo.py_func.__code__)
            self.assertIsInstance(unwind_call.args[2], ValueError)
            self.assertIn(msg, str(unwind_call.args[2]))
        finally:
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.register_callback(tool_id, PY_RETURN, None)
            sys.monitoring.register_callback(tool_id, RAISE, None)
            sys.monitoring.register_callback(tool_id, PY_UNWIND, None)
            sys.monitoring.free_tool_id(tool_id)

        self.assertIn(msg, str(store_raised))

    def test_stop_iteration_under_monitoring(self):
        # Check that Numba can raise an StopIteration exception whilst
        # monitoring is running and that:
        # 1. RAISE is issued for an explicitly raised StopIteration exception.
        # 2. PY_RETURN is issued appropriately for the unwinding stack
        # 3. STOP_ITERATION is not issued as there is no implicit StopIteration
        #    raised.

        return_callback = Mock()
        raise_callback = Mock()
        stopiter_callback = Mock()

        msg = 'exception raised'

        @_enable_sysmon
        @jit('()')
        def foo():
            raise StopIteration(msg)

        store_raised = None
        try:
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.register_callback(tool_id, PY_RETURN,
                                             return_callback)
            sys.monitoring.register_callback(tool_id, RAISE,
                                             raise_callback)
            sys.monitoring.register_callback(tool_id, STOP_ITERATION,
                                             stopiter_callback)
            sys.monitoring.set_events(tool_id,
                                      PY_RETURN | STOP_ITERATION | RAISE)
            try:
                foo()
            except StopIteration as raises:
                store_raised = raises
            # switch off monitoring
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            # check that the return_callback was called once (by Numba unpickle
            # to fetch the exception info out of the stored bytes).
            return_callback.assert_called_once()
            # and that elements that are feasible to check about the call are
            # as expected
            the_call = return_callback.call_args_list[0]
            self.assertEqual(the_call.args[0], _numba_unpickle.__code__)
            self.assertEqual(the_call.args[2][0], StopIteration)
            self.assertEqual(the_call.args[2][1][0], msg)

            # check that the RAISE event callback was triggered
            raise_callback.assert_called()
            # check that it's 3 long (numba unpickle, jit(foo), the test method)
            self.assertEqual(raise_callback.call_count, 3)

            # check the numba pickle call
            numba_unpickle_call = raise_callback.call_args_list[0]
            self.assertEqual(numba_unpickle_call.args[0],
                             _numba_unpickle.__code__)
            self.assertIsInstance(numba_unpickle_call.args[2], KeyError)

            # check the jit(foo) call
            foo_call = raise_callback.call_args_list[1]
            self.assertEqual(foo_call.args[0], foo.py_func.__code__)
            self.assertIsInstance(foo_call.args[2], StopIteration)
            self.assertIn(msg, str(foo_call.args[2]))

            # check the test method call
            meth_call = raise_callback.call_args_list[2]
            test_method_code = sys._getframe().f_code
            self.assertEqual(meth_call.args[0], test_method_code)
            self.assertIsInstance(meth_call.args[2], StopIteration)
            self.assertIn(msg, str(meth_call.args[2]))

            # check that the STOP_ITERATION event was not triggered
            stopiter_callback.assert_not_called()
        finally:
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.register_callback(tool_id, PY_RETURN, None)
            sys.monitoring.register_callback(tool_id, STOP_ITERATION, None)
            sys.monitoring.register_callback(tool_id, RAISE, None)
            sys.monitoring.free_tool_id(tool_id)

        self.assertIn(msg, str(store_raised))

    def test_raising_callback_unwinds_from_jit_on_success_path(self):
        # An event callback can legitimately raise an exception, this test
        # makes sure Numba's dispatcher handles it ok on the "successful path",
        # i.e. the JIT compiled function didn't raise an exception at runtime.

        msg = "deliberately broken callback"

        callback = Mock(side_effect=ValueError(msg))

        store_raised = None
        try:
            event = PY_START
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.set_events(tool_id, event)
            sys.monitoring.register_callback(tool_id, event, callback)
            self.foo(self.arg)
        except ValueError as raises:
            store_raised = raises
        finally:
            sys.monitoring.register_callback(tool_id, event, None)
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.free_tool_id(tool_id)

        callback.assert_called_once()
        self.assertIn(msg, str(store_raised))

    def test_raising_callback_unwinds_from_jit_on_raising_path(self):
        # An event callback can legitimately raise an exception, this test
        # makes sure Numba's dispatcher handles it ok on the
        # "unsuccessful path", i.e. the JIT compiled function raised an
        # exception at runtime. This test checks the RAISE event, as the
        # callback itself raises, it overrides the exception coming from the
        # JIT compiled function.

        msg_callback = "deliberately broken callback"
        msg_execution = "deliberately broken execution"

        callback = Mock(side_effect=ValueError(msg_callback))

        class LocalException(Exception):
            pass

        @_enable_sysmon
        @jit("()")
        def raising():
            raise LocalException(msg_execution)

        store_raised = None
        try:
            event = RAISE
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.set_events(tool_id, event)
            sys.monitoring.register_callback(tool_id, event, callback)
            raising()
        except ValueError as raises:
            store_raised = raises
        finally:
            sys.monitoring.register_callback(tool_id, event, None)
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.free_tool_id(tool_id)

        callback.assert_called()
        # Called 3x (numba unpickle, ValueError in callback, the test method)
        self.assertEqual(callback.call_count, 3)

        # check the numba unpickle call
        numba_unpickle_call = callback.call_args_list[0]
        self.assertEqual(numba_unpickle_call.args[0], _numba_unpickle.__code__)
        self.assertIsInstance(numba_unpickle_call.args[2], KeyError)

        # check the jit(raising) call
        raising_call = callback.call_args_list[1]
        self.assertEqual(raising_call.args[0], raising.py_func.__code__)
        self.assertIs(raising_call.args[2], callback.side_effect)

        # check the test method call
        meth_call = callback.call_args_list[2]
        test_method_code = sys._getframe().f_code
        self.assertEqual(meth_call.args[0], test_method_code)
        self.assertIs(meth_call.args[2], callback.side_effect)

        # check the stored exception is the expected exception
        self.assertIs(store_raised, callback.side_effect)

    def test_raising_callback_unwinds_from_jit_on_unwind_path(self):
        # An event callback can legitimately raise an exception, this test
        # makes sure Numba's dispatcher handles it ok on the
        # "unsuccessful path", i.e. the JIT compiled function raised an
        # exception at runtime. This test checks the PY_UNWIND event. CPython
        # seems to not notice the PY_UNWIND coming from the exception arising
        # from the raise in the event callback, it just has the PY_UNWIND from
        # the raise in the JIT compiled function.

        msg_callback = "deliberately broken callback"
        msg_execution = "deliberately broken execution"

        callback = Mock(side_effect=ValueError(msg_callback))

        class LocalException(Exception):
            pass

        @_enable_sysmon
        @jit("()")
        def raising():
            raise LocalException(msg_execution)

        store_raised = None
        try:
            event = PY_UNWIND
            tool_id = self.tool_id
            sys.monitoring.use_tool_id(tool_id, "custom_monitor")
            sys.monitoring.set_events(tool_id, event)
            sys.monitoring.register_callback(tool_id, event, callback)
            raising()
        except ValueError as raises:
            store_raised = raises
        finally:
            sys.monitoring.register_callback(tool_id, event, None)
            sys.monitoring.set_events(tool_id, NO_EVENTS)
            sys.monitoring.free_tool_id(tool_id)

        callback.assert_called_once()

        # check the jit(raising) call
        raising_call = callback.call_args_list[0]
        self.assertEqual(raising_call.args[0], raising.py_func.__code__)
        self.assertEqual(type(raising_call.args[2]), LocalException)
        self.assertEqual(str(raising_call.args[2]), msg_execution)

        # check the stored_raise
        self.assertIs(store_raised, callback.side_effect)

    def test_monitoring_multiple_threads(self):
        # Two threads, different tools and events registered on each thread.
        # Each test creates a global event capturing. The threads use barriers
        # to wait for each other to start and stop capturing. This way they
        # see the events from each other. One thread is capturing PY_START
        # and the other is capturing PY_RETURN.
        barrier = threading.Barrier(2)

        def barrier_cb():
            barrier.wait()

        def t1_work(self, q):
            try:
                # test event PY_START on a "debugger tool"
                cb = self.run_with_events(self.call_foo, (self.arg,),
                                          (PY_START,),
                                          tool_id=sys.monitoring.DEBUGGER_ID,
                                          barrier=barrier_cb)
                # Check...
                self.assertEqual(len(cb), 1)
                self.check_py_start_calls_multithreads(cb)
            except Exception as e:
                q.put(''.join(traceback.format_exception(e)))

        def t2_work(self, q):
            try:
                # test event PY_RETURN on a "coverage tool"
                cb = self.run_with_events(self.call_foo, (self.arg,),
                                          (PY_RETURN,),
                                          tool_id=sys.monitoring.COVERAGE_ID,
                                          barrier=barrier_cb)
                # Check...
                self.assertEqual(len(cb), 1)
                self.check_py_return_calls_multithreads(cb)
            except Exception as e:
                q.put(''.join(traceback.format_exception(e)))

        q1 = queue.Queue()
        t1 = threading.Thread(target=t1_work, args=(self, q1))
        q2 = queue.Queue()
        t2 = threading.Thread(target=t2_work, args=(self, q2))

        threads = (t1, t2)
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # make sure there were no exceptions
        def assert_empty_queue(q):
            if q.qsize() != 0:
                while not q.empty():
                    print(q.get())
                self.fail("queue supposed to be empty")

        assert_empty_queue(q1)
        assert_empty_queue(q2)


@unittest.skipUnless(PYVERSION >= (3, 12), "needs Python 3.12+")
class TestMonitoringSelfTest(TestCase):

    def test_skipping_of_tests_if_monitoring_in_use(self):
        # check that the unit tests in the TestMonitoring class above will skip
        # if there are other monitoring tools registered in the thread (in this
        # case cProfile is used to cause that effect).
        r = self.subprocess_test_runner(TestMonitoring.__module__,
                                        'TestMonitoring',
                                        'test_start_event',
                                        flags={'-m': 'cProfile'})
        self.assertIn("skipped=1", str(r))


@unittest.skipUnless(PYVERSION >= (3, 12), "needs Python 3.12+")
class TestMonitoringEnvVarControl(TestCase):
    @TestCase.run_test_in_subprocess(
        envvars={"NUMBA_ENABLE_SYS_MONITORING": ''})
    def test_default_off(self):
        @jit
        def foo(x):
            return x + 1

        self.assertFalse(foo._enable_sysmon)

    @TestCase.run_test_in_subprocess(
        envvars={"NUMBA_ENABLE_SYS_MONITORING": '0'})
    def test_override_off(self):
        @jit
        def foo(x):
            return x + 1

        self.assertFalse(foo._enable_sysmon)

    @TestCase.run_test_in_subprocess(
        envvars={"NUMBA_ENABLE_SYS_MONITORING": '1'})
    def test_override_on(self):
        @jit
        def foo(x):
            return x + 1

        self.assertTrue(foo._enable_sysmon)


if __name__ == '__main__':
    unittest.main()
