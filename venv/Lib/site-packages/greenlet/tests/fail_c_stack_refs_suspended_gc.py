"""Regression test for GC of a suspended greenlet's C-stack refs (issue #515).

When a greenlet suspends while the interpreter is holding a ``_PyCStackRef``
(for example, mid attribute resolution), that deferred reference sits on the
greenlet's C stack. The free-threaded collector only walks the running
thread's list in ``gc_visit_thread_stacks``, so unless greenlet visits the
suspended one from ``tp_traverse`` an object reachable only through it is freed
early and used after free once the greenlet resumes.

Reproducing that needs the pinned object to be deferred-refcounted (so its
refcount doesn't keep it alive) and reachable only through the C-stack ref. A
class is deferred-refcounted, and a metaclass ``__get__`` makes a class act as a
descriptor -- so ``obj.attr``, where ``attr`` is such a class, pins the class in
a ``_PyCStackRef`` across the (Python) ``__get__``. We switch away from inside
``__get__``, drop every other reference to the class, and collect: on a
regressed free-threaded build the class is gone by the time we resume; the fix
keeps it alive because ``tp_traverse`` visited the snapshot taken at suspend.

Prints "C STACK REFS GC OK" on a fixed build (and on with-GIL builds, where
ordinary refcounting keeps the class alive); a regressed free-threaded build
reports the class was collected and exits non-zero.
"""
import gc
import sys
import weakref

import greenlet

parent = greenlet.getcurrent()
observed = {}


class Meta(type):
    def __get__(cls, obj, objtype=None):
        # type(cls) is Meta, which defines __get__, so ``cls`` is a descriptor;
        # the getattr machinery pins it in a _PyCStackRef across this call. A
        # class is deferred-refcounted on a free-threaded build, so only that
        # (deferred) C-stack ref will be keeping it alive in a moment.
        #
        # Drop the locals the descriptor protocol gave us (the class is ``cls``):
        # greenlet visits a suspended greenlet's frames, so a leftover frame ref
        # would keep the class alive on its own and mask the bug.
        del cls, obj, objtype
        child.switch()
        return 42


pinned = Meta('pinned', (), {})                    # a class => deferred-refcounted
holder = type('holder', (), {'attr': pinned})      # holder.attr invokes Meta.__get__
box = [pinned]
del pinned


def child_work():
    # ``parent`` is suspended inside Meta.__get__, holding a C-stack ref to the
    # class. Drop every *other* reference to it, then collect: now only greenlet
    # visiting the suspended C-stack ref can keep the class alive.
    ref = weakref.ref(box[0])
    box[0] = None
    del holder.attr
    for _ in range(5):
        gc.collect()
    observed['alive'] = ref() is not None
    parent.switch()


child = greenlet.greenlet(child_work)
result = holder().attr                             # Meta.__get__ -> switch -> gc -> back
assert result == 42, result

alive = observed.get('alive')
print(f"py={sys.version.split()[0]} gil={getattr(sys, '_is_gil_enabled', lambda: True)()} "
      f"greenlet={greenlet.__version__} alive={alive}", flush=True)
if not alive:
    raise SystemExit("REGRESSED: a class reachable only through a suspended "
                     "greenlet's C-stack ref was collected early")
print("C STACK REFS GC OK", flush=True)
