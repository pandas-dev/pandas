"""Patch asyncio to allow nested event loops."""

import asyncio
import asyncio.events as events
import os
import sys
import threading
from contextlib import contextmanager, suppress
from heapq import heappop

_run_close_loop = True

class _NestAsyncio2:
    '''Internal class of `nest_asyncio2`.
     
    Mainly for holding the original properties to support unapply() and nest_asyncio2.run().
    '''
    pass

def apply(
    loop=None,
    *,
    run_close_loop: bool = False,
    error_on_mispatched: bool = False
):
    '''Patch asyncio to make its event loop reentrant.
    
    - `run_close_loop`: Close the event loop created by `asyncio.run()`, if any.
      See README for details.
    - `error_on_mispatched`:
      - `False` (default): Warn if asyncio is already patched by `nest_asyncio` on Python 3.12+.
      - `True`: Raise `RuntimeError` if asyncio is already patched by `nest_asyncio`.
    '''
    global _run_close_loop
    
    _patch_asyncio(error_on_mispatched=error_on_mispatched)
    _patch_policy()
    _patch_tornado()

    loop = loop or _get_event_loop()
    if loop is not None:
        _patch_loop(loop)

    _run_close_loop &= run_close_loop

if sys.version_info < (3, 12, 0):
    def _get_event_loop():
        return asyncio.get_event_loop()
elif sys.version_info < (3, 14, 0):
    def _get_event_loop():
        # Python 3.12~3.13:
        # Calling get_event_loop() will result in ResourceWarning: unclosed event loop
        loop = events._get_running_loop()
        if loop is None:
            policy = events.get_event_loop_policy()
            loop = policy._local._loop
        return loop
else:
    def _get_event_loop():
        # Python 3.14: Raises a RuntimeError if there is no current event loop.
        try:
            return asyncio.get_event_loop()
        except RuntimeError:
            return None

if sys.version_info < (3, 12, 0):    
    def run(main, *, debug=False):
        loop = asyncio.get_event_loop()
        loop.set_debug(debug)
        task = asyncio.ensure_future(main)
        try:
            return loop.run_until_complete(task)
        finally:
            if not task.done():
                task.cancel()
                with suppress(asyncio.CancelledError):
                    loop.run_until_complete(task)
else:
    def run(main, *, debug=False, loop_factory=None):
        new_event_loop = False
        set_event_loop = None
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            # if sys.version_info < (3, 16, 0):
            #     policy = asyncio.events._get_event_loop_policy()
            #     try:
            #         loop = policy.get_event_loop()
            #     except RuntimeError:
            #         loop = loop_factory()
            # else:
            #     loop = loop_factory()
            if not _run_close_loop:
                # Not running
                loop = _get_event_loop()
                if loop is None:
                    if loop_factory is None:
                        loop_factory = asyncio.new_event_loop
                    loop = loop_factory()
                    asyncio.set_event_loop(loop)
            else:
                if loop_factory is None:
                    loop = asyncio.new_event_loop()
                    # Not running
                    set_event_loop = _get_event_loop()
                    asyncio.set_event_loop(loop)
                else:
                    loop = loop_factory()
                new_event_loop = True
        _patch_loop(loop)

        loop.set_debug(debug)
        task = asyncio.ensure_future(main, loop=loop)
        try:
            return loop.run_until_complete(task)
        finally:
            if not task.done():
                task.cancel()
                with suppress(asyncio.CancelledError):
                    loop.run_until_complete(task)
            if set_event_loop:
                # asyncio.Runner just set_event_loop(None) but we are nested
                asyncio.set_event_loop(set_event_loop)
            if new_event_loop:
                # Avoid ResourceWarning: unclosed event loop
                loop.close()

def _patch_asyncio(*, error_on_mispatched: bool = False):
    """Patch asyncio module to use pure Python tasks and futures."""

    def _get_event_loop(stacklevel=3):
        loop = events._get_running_loop()
        if loop is None:
            loop = events.get_event_loop_policy().get_event_loop()
        return loop

    # Use module level _current_tasks, all_tasks and patch run method.
    if hasattr(asyncio, '_nest_patched'):
        if not hasattr(asyncio, '_nest_asyncio2'):
            if error_on_mispatched:
                raise RuntimeError('asyncio is already patched by nest_asyncio')
            elif sys.version_info >= (3, 12, 0):
                import warnings
                warnings.warn('asyncio is already patched by nest_asyncio. You may encounter bugs related to asyncio')
        return
    
    # Using _PyTask on Python 3.14+ will break current_task() (and all_tasks(),
    # _swap_current_task())
    # Even we replace it with _py_current_task(), it only works with _PyTask, but
    # the external loop is probably using _CTask.
    # https://github.com/python/cpython/pull/129899
    if sys.version_info >= (3, 6, 0) and sys.version_info < (3, 14, 0):
        asyncio.Task = asyncio.tasks._CTask = asyncio.tasks.Task = \
            asyncio.tasks._PyTask
        asyncio.Future = asyncio.futures._CFuture = asyncio.futures.Future = \
            asyncio.futures._PyFuture
    if sys.version_info < (3, 7, 0):
        asyncio.tasks._current_tasks = asyncio.tasks.Task._current_tasks
        asyncio.all_tasks = asyncio.tasks.Task.all_tasks
    # The same as asyncio.get_event_loop() on at least Python 3.14
    if sys.version_info >= (3, 9, 0) and sys.version_info < (3, 14, 0):
        events._get_event_loop = events.get_event_loop = \
            asyncio.get_event_loop = _get_event_loop
    asyncio.run = run
    asyncio._nest_patched = True
    asyncio._nest_asyncio2 = _NestAsyncio2()


def _patch_policy():
    """Patch the policy to always return a patched loop."""

    # Python 3.14:
    # get_event_loop() raises a RuntimeError if there is no current event loop.
    # So there is no need to _patch_loop() in it.
    # Patching new_event_loop() may be better, but policy is going to be removed...
    # Removed in Python 3.16
    # https://github.com/python/cpython/issues/127949
    if sys.version_info >= (3, 14, 0):
        return

    def get_event_loop(self):
        if self._local._loop is None:
            loop = self.new_event_loop()
            _patch_loop(loop)
            self.set_event_loop(loop)
        return self._local._loop

    if sys.version_info < (3, 14, 0):
        policy = events.get_event_loop_policy()
    else:
        policy = events._get_event_loop_policy()
    policy.__class__.get_event_loop = get_event_loop


def _patch_loop(loop):
    """Patch loop to make it reentrant."""

    def run_forever(self):
        with manage_run(self), manage_asyncgens(self):
            while True:
                self._run_once()
                if self._stopping:
                    break
        self._stopping = False

    def run_until_complete(self, future):
        with manage_run(self):
            f = asyncio.ensure_future(future, loop=self)
            if f is not future:
                f._log_destroy_pending = False
            while not f.done():
                self._run_once()
                if self._stopping:
                    break
            if not f.done():
                raise RuntimeError(
                    'Event loop stopped before Future completed.')

            # When a task completes inside _run_once(), Task.__step calls
            # Future.set_result(), which schedules all done callbacks via call_soon().
            # These callbacks are added to loop._ready but the `while not f.done()` loop
            # exits immediately — before the next _run_once() that would process them.
            # https://github.com/Chaoses-Ib/nest-asyncio2/issues/3

            # Process any callbacks scheduled during task completion
            # (e.g. task done callbacks added via call_soon in set_result).
            if self._ready:
                self._run_once()

            return f.result()

    def _run_once(self):
        """
        Simplified re-implementation of asyncio's _run_once that
        runs handles as they become ready.
        """
        ready = self._ready
        scheduled = self._scheduled
        while scheduled and scheduled[0]._cancelled:
            heappop(scheduled)

        timeout = (
            0 if ready or self._stopping
            else min(max(
                scheduled[0]._when - self.time(), 0), 86400) if scheduled
            else None)
        event_list = self._selector.select(timeout)
        self._process_events(event_list)

        end_time = self.time() + self._clock_resolution
        while scheduled and scheduled[0]._when < end_time:
            handle = heappop(scheduled)
            ready.append(handle)

        for _ in range(len(ready)):
            if not ready:
                break
            handle = ready.popleft()
            if not handle._cancelled:
                # preempt the current task so that that checks in
                # Task.__step do not raise
                if sys.version_info < (3, 14, 0):
                    curr_task = curr_tasks.pop(self, None)
                else:
                    # Work with both C and Py
                    try:
                        curr_task = asyncio.tasks._swap_current_task(self, None)
                    except KeyError:
                        curr_task = None

                try:
                    handle._run()
                finally:
                    # restore the current task
                    if curr_task is not None:
                        if sys.version_info < (3, 14, 0):
                            curr_tasks[self] = curr_task
                        else:
                            # Work with both C and Py
                            asyncio.tasks._swap_current_task(self, curr_task)

        handle = None

    @contextmanager
    def manage_run(self):
        """Set up the loop for running."""
        self._check_closed()
        old_thread_id = self._thread_id
        old_running_loop = events._get_running_loop()
        try:
            self._thread_id = threading.get_ident()
            events._set_running_loop(self)
            self._num_runs_pending += 1
            if self._is_proactorloop:
                if self._self_reading_future is None:
                    self.call_soon(self._loop_self_reading)
            yield
        finally:
            self._thread_id = old_thread_id
            events._set_running_loop(old_running_loop)
            self._num_runs_pending -= 1
            if self._is_proactorloop:
                if (self._num_runs_pending == 0
                        and self._self_reading_future is not None):
                    ov = self._self_reading_future._ov
                    self._self_reading_future.cancel()
                    if ov is not None:
                        self._proactor._unregister(ov)
                    self._self_reading_future = None

    @contextmanager
    def manage_asyncgens(self):
        if not hasattr(sys, 'get_asyncgen_hooks'):
            # Python version is too old.
            return
        old_agen_hooks = sys.get_asyncgen_hooks()
        try:
            self._set_coroutine_origin_tracking(self._debug)
            if self._asyncgens is not None:
                sys.set_asyncgen_hooks(
                    firstiter=self._asyncgen_firstiter_hook,
                    finalizer=self._asyncgen_finalizer_hook)
            yield
        finally:
            self._set_coroutine_origin_tracking(False)
            if self._asyncgens is not None:
                sys.set_asyncgen_hooks(*old_agen_hooks)

    def _check_running(self):
        """Do not throw exception if loop is already running."""
        pass

    if hasattr(loop, '_nest_patched'):
        return
    if not isinstance(loop, asyncio.BaseEventLoop):
        raise ValueError('Can\'t patch loop of type %s' % type(loop))
    cls = loop.__class__
    cls.run_forever = run_forever
    cls.run_until_complete = run_until_complete
    cls._run_once = _run_once
    cls._check_running = _check_running
    cls._check_runnung = _check_running  # typo in Python 3.7 source
    cls._num_runs_pending = 1 if loop.is_running() else 0
    cls._is_proactorloop = (
        os.name == 'nt' and issubclass(cls, asyncio.ProactorEventLoop))
    if sys.version_info < (3, 7, 0):
        cls._set_coroutine_origin_tracking = cls._set_coroutine_wrapper
    curr_tasks = asyncio.tasks._current_tasks \
        if sys.version_info >= (3, 7, 0) else asyncio.Task._current_tasks
    cls._nest_patched = True
    cls._nest_asyncio2 = _NestAsyncio2()


def _patch_tornado():
    """
    If tornado is imported before nest_asyncio, make tornado aware of
    the pure-Python asyncio Future.
    """
    if 'tornado' in sys.modules:
        import tornado.concurrent as tc  # type: ignore
        tc.Future = asyncio.Future
        if asyncio.Future not in tc.FUTURES:
            tc.FUTURES += (asyncio.Future,)
