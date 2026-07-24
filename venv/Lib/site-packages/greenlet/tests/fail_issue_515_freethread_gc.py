"""Reproducer for the free-threaded GC crash in issue #515.

https://github.com/python-greenlet/greenlet/issues/515

Before the fix, ``set_initial_state`` copied the parent thread state's
``_PyCStackRef`` list head into each new greenlet. Those nodes sit on the
parent greenlet's C stack, so once the child starts running on its own stack
and overwrites that region their ``next`` pointers dangle. The free-threaded
collector walks ``c_stack_refs`` for every thread, so a GC on any thread ends
up chasing those dangling nodes (inside ``gc_collect_main`` ->
``gc_visit_thread_stacks``) while such a child is active, and segfaults.

The fiddly part is getting a *non-empty* list to inherit at all. The
interpreter only holds a ``_PyCStackRef`` transiently, while it is resolving an
attribute -- ``_Py_LoadAttr_StackRefSteal`` pins ``self`` across ``obj.attr()``
-- so the greenlet has to be primed from inside that window. Priming it from a
property getter does exactly that; it is also the same deep-attribute context
Playwright's sync bridge happens to start its fiber from. A churn thread then
keeps the collector busy while the dispatcher runs on its stale list.

A fixed build prints "ISSUE 515 OK"; a regressed one segfaults in about a
second. With-GIL builds are unaffected -- ordinary refcounting keeps the nodes
alive -- so there it simply prints the OK line.
"""
import gc
import os
import sys
import threading

import greenlet

# The crash lands within a second, so a few seconds leaves plenty of margin on
# a slow machine without dragging out a fixed build.
DURATION = float(os.environ.get("REPRO_SECONDS", "4"))
NWORKERS = int(os.environ.get("REPRO_THREADS", "6"))


def _busy():
    # Keep the dispatcher doing work (resolving str methods) so a collection on
    # the churn thread catches it while it is the active greenlet.
    return "abc".upper().lower().strip().title()


class Bridge:
    def __init__(self, stop):
        self.stop = stop
        self.main = greenlet.getcurrent()
        self.disp = greenlet.greenlet(self._loop)

    def _loop(self):
        while not self.stop[0]:
            for _ in range(50):
                _busy()
            self.main.switch()

    @property
    def _prime(self):
        # Reaching this getter means we are mid attribute-resolution, so a
        # _PyCStackRef for ``self`` is live right now. Priming the dispatcher
        # from here is what makes it inherit a non-empty c_stack_refs head.
        self.disp.switch()
        return None

    def start(self):
        # Reading the property (rather than calling a method) is deliberate: the
        # getter runs *during* attribute resolution, so it primes the dispatcher
        # while the _PyCStackRef for ``self`` is still held. We return the value
        # only so this isn't a bare, pointless-looking expression statement.
        return self._prime

    def resume(self):
        self.disp.switch()


def worker(stop):
    bridge = Bridge(stop)
    bridge.start()
    # start() has returned, so the stack the inherited nodes point at is being
    # reused -- the dispatcher's list is stale now. Keep it active to trip the GC.
    while not stop[0]:
        bridge.resume()


def churn(stop):
    while not stop[0]:
        # A cycle forces a real (stop-the-world) collection rather than just a
        # refcount drop.
        nodes = [{"i": i, "next": None} for i in range(500)]
        for a, b in zip(nodes, nodes[1:]):
            a["next"] = b
        nodes[-1]["next"] = nodes[0]
        del nodes
        gc.collect()


def main():
    stop = [False]
    threads = [threading.Thread(target=churn, args=(stop,), daemon=True)]
    threads += [threading.Thread(target=worker, args=(stop,), daemon=True)
                for _ in range(NWORKERS)]
    for t in threads:
        t.start()
    threading.Event().wait(DURATION)
    stop[0] = True
    for t in threads:
        t.join(timeout=5)
    print("ISSUE 515 OK")


if __name__ == "__main__":
    main()
    sys.exit(0)
