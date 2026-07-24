"""Out-of-process payload for the free-threaded switch/critical-section deadlock.

See ``freethread_switch_deadlock.py`` in the repo root for the full write-up.
Before the fix, a greenlet switch swapped C stacks but left any held
``PyCriticalSection`` locks (tracked in ``tstate->critical_section``) locked.
asyncio's ``Task.__step`` holds such a lock on the running task across the step,
so switching into a child greenlet and touching that task from there blocked
forever. A fixed build (and any GIL-enabled build) prints the sentinel; a
regressed free-threaded build deadlocks, and the watchdog turns that hang into a
non-zero exit so the test fails loudly instead of stalling the whole suite.
"""
import asyncio
import faulthandler
import sys

import greenlet

# The deadlock is immediate when present; the generous timeout is only so a
# genuinely regressed build still exits on the slowest CI, never on a good one.
faulthandler.dump_traceback_later(15, exit=True)


async def main():
    task = asyncio.current_task()  # its running __step holds a lock on `task`

    def in_child_fiber():
        # Fresh C stack; the task's lock is still held by the fiber we left.
        # A regressed build never returns from this first call.
        task.add_done_callback(lambda _: None)
        task.remove_done_callback(lambda _: None)

    greenlet.greenlet(in_child_fiber).switch()


asyncio.run(main())
faulthandler.cancel_dump_traceback_later()
print("SWITCH CS OK")
sys.exit(0)
