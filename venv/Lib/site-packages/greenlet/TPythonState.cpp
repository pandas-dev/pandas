#ifndef GREENLET_PYTHON_STATE_CPP
#define GREENLET_PYTHON_STATE_CPP

#include <Python.h>
#include "TGreenlet.hpp"

namespace greenlet {

PythonState::PythonState()
    : _top_frame()
#if GREENLET_USE_CFRAME
    ,cframe(nullptr)
    ,use_tracing(0)
#endif
#if GREENLET_PY314
    ,py_recursion_depth(0)
    ,current_executor(nullptr)
    ,stackpointer(nullptr)
    #ifdef Py_GIL_DISABLED
    ,c_stack_refs(nullptr)
    #endif
#elif GREENLET_PY312
    ,py_recursion_depth(0)
    ,c_recursion_depth(0)
#else
    ,recursion_depth(0)
#endif
#if GREENLET_PY313
    ,delete_later(nullptr)
    ,critical_section(0)
#else
    ,trash_delete_nesting(0)
#endif
#if GREENLET_PY311
    ,current_frame(nullptr)
    ,datastack_chunk(nullptr)
    ,datastack_top(nullptr)
    ,datastack_limit(nullptr)
#endif
{
#if GREENLET_USE_CFRAME
    /*
      The PyThreadState->cframe pointer usually points to memory on
      the stack, alloceted in a call into PyEval_EvalFrameDefault.

      Initially, before any evaluation begins, it points to the
      initial PyThreadState object's ``root_cframe`` object, which is
      statically allocated for the lifetime of the thread.

      A greenlet can last for longer than a call to
      PyEval_EvalFrameDefault, so we can't set its ``cframe`` pointer
      to be the current ``PyThreadState->cframe``; nor could we use
      one from the greenlet parent for the same reason. Yet a further
      no: we can't allocate one scoped to the greenlet and then
      destroy it when the greenlet is deallocated, because inside the
      interpreter the _PyCFrame objects form a linked list, and that too
      can result in accessing memory beyond its dynamic lifetime (if
      the greenlet doesn't actually finish before it dies, its entry
      could still be in the list).

      Using the ``root_cframe`` is problematic, though, because its
      members are never modified by the interpreter and are set to 0,
      meaning that its ``use_tracing`` flag is never updated. We don't
      want to modify that value in the ``root_cframe`` ourself: it
      *shouldn't* matter much because we should probably never get
      back to the point where that's the only cframe on the stack;
      even if it did matter, the major consequence of an incorrect
      value for ``use_tracing`` is that if its true the interpreter
      does some extra work --- however, it's just good code hygiene.

      Our solution: before a greenlet runs, after its initial
      creation, it uses the ``root_cframe`` just to have something to
      put there. However, once the greenlet is actually switched to
      for the first time, ``g_initialstub`` (which doesn't actually
      "return" while the greenlet is running) stores a new _PyCFrame on
      its local stack, and copies the appropriate values from the
      currently running _PyCFrame; this is then made the _PyCFrame for the
      newly-minted greenlet. ``g_initialstub`` then proceeds to call
      ``glet.run()``, which results in ``PyEval_...`` adding the
      _PyCFrame to the list. Switches continue as normal. Finally, when
      the greenlet finishes, the call to ``glet.run()`` returns and
      the _PyCFrame is taken out of the linked list and the stack value
      is now unused and free to expire.

      XXX: I think we can do better. If we're deallocing in the same
      thread, can't we traverse the list and unlink our frame?
      Can we just keep a reference to the thread state in case we
      dealloc in another thread? (Is that even possible if we're still
      running and haven't returned from g_initialstub?)
    */
    this->cframe = &PyThreadState_GET()->root_cframe;
#endif
}

#if GREENLET_PY314 && defined(Py_GIL_DISABLED)
void PythonState::capture_c_stack_refs(const PyThreadState* tstate) noexcept
{
    // Runs from operator<< while our C stack is still live and coherent, so we
    // can walk tstate's _PyCStackRef list and take a strong reference to every
    // object it holds. tp_traverse visits these once we're suspended, because
    // by then the nodes themselves have been relocated into the heap stack copy
    // and the saved list head no longer points at them. Strong references (not
    // _Py_VISIT_STACKREF, whose _PyGC_VisitStackRef isn't exported before 3.15);
    // a std::vector rather than a Python list/tuple because operator<< must not
    // allocate a GC-tracked object mid-switch. Rebuilt from scratch each time;
    // the list is empty at a typical switch, so this is usually just an empty
    // loop.
    this->c_stack_ref_snapshot.clear();
    for (const _PyCStackRef* node = ((_PyThreadStateImpl*)tstate)->c_stack_refs;
         node != nullptr; node = node->next) {
        if (!PyStackRef_IsNullOrInt(node->ref)) {
            this->c_stack_ref_snapshot.push_back(
                OwnedObject::owning(PyStackRef_AsPyObjectBorrow(node->ref)));
        }
    }
}
#endif


inline void PythonState::may_switch_away() noexcept
{
#if GREENLET_PY311
    // PyThreadState_GetFrame is probably going to have to allocate a
    // new frame object. That may trigger garbage collection. Because
    // we call this during the early phases of a switch (it doesn't
    // matter to which greenlet, as this has a global effect), if a GC
    // triggers a switch away, two things can happen, both bad:
    // - We might not get switched back to, halting forward progress.
    //   this is pathological, but possible.
    // - We might get switched back to with a different set of
    //   arguments or a throw instead of a switch. That would corrupt
    //   our state (specifically, PyErr_Occurred() and this->args()
    //   would no longer agree).
    //
    // Thus, when we call this API, we need to have GC disabled.
    // This method serves as a bottleneck we call when maybe beginning
    // a switch. In this way, it is always safe -- no risk of GC -- to
    // use ``_GetFrame()`` whenever we need to, just as it was in
    // <=3.10 (because subsequent calls will be cached and not
    // allocate memory).

    GCDisabledGuard no_gc;
    Py_XDECREF(PyThreadState_GetFrame(PyThreadState_GET()));
#endif
}

void PythonState::operator<<(const PyThreadState *const tstate) noexcept
{
    this->_context.steal(tstate->context);
#if GREENLET_USE_CFRAME
    /*
      IMPORTANT: ``cframe`` is a pointer into the STACK. Thus, because
      the call to ``slp_switch()`` changes the contents of the stack,
      you cannot read from ``ts_current->cframe`` after that call and
      necessarily get the same values you get from reading it here.
      Anything you need to restore from now to then must be saved in a
      global/threadlocal variable (because we can't use stack
      variables here either). For things that need to persist across
      the switch, use `will_switch_from`.
    */
    this->cframe = tstate->cframe;
  #if !GREENLET_PY312
    this->use_tracing = tstate->cframe->use_tracing;
  #endif
#endif // GREENLET_USE_CFRAME
#if GREENLET_PY311
  #if GREENLET_PY314
    this->py_recursion_depth = tstate->py_recursion_limit - tstate->py_recursion_remaining;
    this->current_executor = tstate->current_executor;
    #ifdef Py_GIL_DISABLED
    this->c_stack_refs = ((_PyThreadStateImpl*)tstate)->c_stack_refs;
    // Capture the deferred references now, while our C stack is still live, so
    // tp_traverse can keep them from being collected while we're suspended.
    this->capture_c_stack_refs(tstate);
    #endif
  #elif GREENLET_PY312
    this->py_recursion_depth = tstate->py_recursion_limit - tstate->py_recursion_remaining;
    this->c_recursion_depth = Py_C_RECURSION_LIMIT - tstate->c_recursion_remaining;
  #else // not 312
    this->recursion_depth = tstate->recursion_limit - tstate->recursion_remaining;
  #endif // GREENLET_PY312
  #if GREENLET_PY313
    this->current_frame = tstate->current_frame;
  #elif GREENLET_USE_CFRAME
    this->current_frame = tstate->cframe->current_frame;
  #endif
    this->datastack_chunk = tstate->datastack_chunk;
    this->datastack_top = tstate->datastack_top;
    this->datastack_limit = tstate->datastack_limit;

    PyFrameObject *frame = PyThreadState_GetFrame((PyThreadState *)tstate);
    Py_XDECREF(frame);  // PyThreadState_GetFrame gives us a new
                        // reference.
    this->_top_frame.steal(frame);
  #if GREENLET_PY314
    if (this->top_frame()) {
        this->stackpointer = this->_top_frame->f_frame->stackpointer;
    }
    else {
        this->stackpointer = nullptr;
    }
  #endif
  #if GREENLET_PY313
    // By contract of _PyTrash_thread_deposit_object,
    // the ``delete_later`` object has a refcount of 0.
    // We take a strong reference to it.
    //
    // Now, ``delete_later`` is managed as a
    // linked list whose objects are unconditionally deallocated
    // WITHOUT calling DECREF on them, so it's not clear what that is
    // actually accomplishing. That is, if another object is pushed on
    // the list and then the list is deallocated, this object will
    // still be deallocated. This strong reference serves as a form of
    // resurrection, meaning that when operator>> DECREFs it, we might
    // enter its ``tp_dealloc`` function again.
    //
    // In practice, it's quite difficult to arrange for this to be
    // a non-null value during a greenlet switch.
    // ``greenlet.tests.test_greenlet_trash`` tries, but under 3.14,
    // at least, fails to do so.
    this->delete_later = Py_XNewRef(tstate->delete_later);
#ifdef Py_GIL_DISABLED
    // Switching greenlets swaps C stacks, which to the free-threaded runtime is
    // the same predicament as detaching the thread: the PyCriticalSection nodes
    // chained off tstate->critical_section live on the stack we're leaving, and
    // their PyMutexes would stay locked behind our back. The greenlet we switch
    // to could then block forever taking one of those same locks -- e.g. an
    // asyncio event dispatched onto another fiber re-enters a Task/Future that
    // the suspended fiber is mid-step on. So drop the locks here the way
    // _PyThreadState_Detach() does and let operator>> re-take them on resume.
    if (tstate->critical_section != 0) {
        _PyCriticalSection_SuspendAll(const_cast<PyThreadState*>(tstate));
    }
#endif
    this->critical_section = tstate->critical_section;
  #elif GREENLET_PY312
    this->trash_delete_nesting = tstate->trash.delete_nesting;
  #else // not 312 or 3.13+
    this->trash_delete_nesting = tstate->trash_delete_nesting;
  #endif // GREENLET_PY312
#else // Not 311
    this->recursion_depth = tstate->recursion_depth;
    this->_top_frame.steal(tstate->frame);
    this->trash_delete_nesting = tstate->trash_delete_nesting;
#endif // GREENLET_PY311
}

#if GREENLET_PY312
void GREENLET_NOINLINE(PythonState::unexpose_frames)()
{
    if (!this->top_frame()) {
        return;
    }

    // See GreenletState::expose_frames() and the comment on frames_were_exposed
    // for more information about this logic.
    _PyInterpreterFrame *iframe = this->_top_frame->f_frame;
    while (iframe != nullptr) {
        _PyInterpreterFrame *prev_exposed = iframe->previous;
        assert(iframe->frame_obj);
        memcpy(&iframe->previous, &iframe->frame_obj->_f_frame_data[0],
               sizeof(void *));
        iframe = prev_exposed;
    }
}
#else
void PythonState::unexpose_frames()
{}
#endif

void PythonState::operator>>(PyThreadState *const tstate) noexcept
{
    tstate->context = this->_context.relinquish_ownership();
    /* Incrementing this value invalidates the contextvars cache,
       which would otherwise remain valid across switches */
    tstate->context_ver++;
#if GREENLET_USE_CFRAME
    tstate->cframe = this->cframe;
    /*
      If we were tracing, we need to keep tracing.
      There should never be the possibility of hitting the
      root_cframe here. See note above about why we can't
      just copy this from ``origin->cframe->use_tracing``.
    */
  #if !GREENLET_PY312
    tstate->cframe->use_tracing = this->use_tracing;
  #endif
#endif // GREENLET_USE_CFRAME
#if GREENLET_PY311
  #if GREENLET_PY314
    tstate->py_recursion_remaining = tstate->py_recursion_limit - this->py_recursion_depth;
    tstate->current_executor = this->current_executor;
    #ifdef Py_GIL_DISABLED
    ((_PyThreadStateImpl*)tstate)->c_stack_refs = this->c_stack_refs;
    // We're the running greenlet again: our C-stack refs live in the thread
    // state now and gc_visit_thread_stacks() covers them, so drop the strong
    // references tp_traverse held on our behalf while we were suspended.
    this->c_stack_ref_snapshot.clear();
    #endif
    this->unexpose_frames();
  #elif GREENLET_PY312
    tstate->py_recursion_remaining = tstate->py_recursion_limit - this->py_recursion_depth;
    tstate->c_recursion_remaining = Py_C_RECURSION_LIMIT - this->c_recursion_depth;
    this->unexpose_frames();
  #else // \/ 3.11
    tstate->recursion_remaining = tstate->recursion_limit - this->recursion_depth;
  #endif // GREENLET_PY312
  #if GREENLET_PY313
    tstate->current_frame = this->current_frame;
  #elif GREENLET_USE_CFRAME
    tstate->cframe->current_frame = this->current_frame;
  #endif
    tstate->datastack_chunk = this->datastack_chunk;
    tstate->datastack_top = this->datastack_top;
    tstate->datastack_limit = this->datastack_limit;
#if GREENLET_PY314 && defined(Py_GIL_DISABLED)
    if (this->top_frame()) {
        this->_top_frame->f_frame->stackpointer = this->stackpointer;
    }
#endif
    this->_top_frame.relinquish_ownership();
  #if GREENLET_PY313
    // See comments in operator<<. We own a strong reference to
    // this->delete_later, which may or may not be the same object as
    // tstate->delete_later (depending if something pushed an object
    // onto the trashcan). Again, because ``delete_later`` is managed
    // as a linked list, it's not clear that saving and restoring the
    // value, especially without ever setting it to NULL, accomplishes
    // much...but the code was added by a core dev, so assume correct.
    //
    // Recall that tstate->delete_later is supposed to have a refcount
    // of 0, because objects are added there from their ``tp_dealloc``
    // method. So we should only need to DECREF it if we're the ones
    // that INCREF'd it in operator<<. (This is different than the
    // core dev's original code which always did this.)
    if (this->delete_later == tstate->delete_later) {
        Py_XDECREF(tstate->delete_later);
        tstate->delete_later = this->delete_later;
        this->delete_later = nullptr;
    }
    else {
        // it got switched behind our back. So the reference we own
        // needs to be explicitly cleared.
        tstate->delete_later = this->delete_later;
        Py_CLEAR(this->delete_later);
    }
    tstate->critical_section = this->critical_section;
#ifdef Py_GIL_DISABLED
    // Re-acquire whatever operator<< suspended when this greenlet last yielded.
    // A no-op for a greenlet that held no locks, and for a brand-new one whose
    // chain starts empty. Mirrors the resume in _PyThreadState_Attach(); note
    // _PyCriticalSection_Resume() dereferences the head, so the != 0 guard is
    // load-bearing, not just a fast path.
    if (tstate->critical_section != 0) {
        _PyCriticalSection_Resume(tstate);
    }
#endif

  #elif GREENLET_PY312
    tstate->trash.delete_nesting = this->trash_delete_nesting;
  #else // not 3.12
    tstate->trash_delete_nesting = this->trash_delete_nesting;
  #endif // GREENLET_PY312
#else // not 3.11
    tstate->frame = this->_top_frame.relinquish_ownership();
    tstate->recursion_depth = this->recursion_depth;
    tstate->trash_delete_nesting = this->trash_delete_nesting;
#endif // GREENLET_PY311
}

inline void PythonState::will_switch_from(PyThreadState *const origin_tstate) noexcept
{
#if GREENLET_USE_CFRAME && !GREENLET_PY312
    // The weird thing is, we don't actually save this for an
    // effect on the current greenlet, it's saved for an
    // effect on the target greenlet. That is, we want
    // continuity of this setting across the greenlet switch.
    this->use_tracing = origin_tstate->cframe->use_tracing;
#endif
}

void PythonState::set_initial_state(const PyThreadState* const tstate) noexcept
{
    this->_top_frame = nullptr;
#if GREENLET_PY314
    this->py_recursion_depth = tstate->py_recursion_limit - tstate->py_recursion_remaining;
    this->current_executor = tstate->current_executor;
    #ifdef Py_GIL_DISABLED
    // Start with an empty C-stack-ref list, the way a brand-new thread does;
    // do NOT copy the parent thread state's head. Those _PyCStackRef nodes sit
    // on the parent greenlet's C stack, so once we start running on our own
    // stack and overwrite that region, following them reads garbage. The
    // free-threaded collector walks c_stack_refs for every thread in
    // gc_visit_thread_stacks(), so leaving the stale head here crashed it.
    // See https://github.com/python-greenlet/greenlet/issues/515.
    this->c_stack_refs = nullptr;
    #endif
    // this->stackpointer is left null because this->_top_frame is
    // null so there is no value to copy.
#elif GREENLET_PY312
    this->py_recursion_depth = tstate->py_recursion_limit - tstate->py_recursion_remaining;
#if GREENLET_314
    this->c_recursion_depth = 0; // unused on 3.14
#else
    this->c_recursion_depth = Py_C_RECURSION_LIMIT - tstate->c_recursion_remaining;
#endif
#elif GREENLET_PY311
    this->recursion_depth = tstate->recursion_limit - tstate->recursion_remaining;
#else
    this->recursion_depth = tstate->recursion_depth;
#endif
}
// TODO: Better state management about when we own the top frame.
int PythonState::tp_traverse(visitproc visit, void* arg, bool visit_top_frame) noexcept
{
    Py_VISIT(this->_context.borrow());
    if (visit_top_frame) {
        Py_VISIT(this->_top_frame.borrow());
    }
#if GREENLET_PY315
    // Visit the references held by our suspended frames.
    // This is important specially on free-threading where the
    // the suspended frames may contain deferred references to
    // objects, and if they are not traversed then the interpreter
    // can free objects early causing a use-after-free crash
    // at runtime exit.
    if (this->_top_frame) {
        for (_PyInterpreterFrame* iframe = this->_top_frame->f_frame;
             iframe != nullptr; iframe = iframe->previous) {
            // Skip generator/coroutine frames; their object's traverse
            // already visits them (gen_traverse), so we'd double-count.
            // expose_frames leaves them in the ->previous chain.
            if (iframe->owner != FRAME_OWNED_BY_THREAD) {
                continue;
            }
            Py_VISIT(iframe->frame_obj);
            Py_VISIT(iframe->f_locals);
            _Py_VISIT_STACKREF(iframe->f_funcobj);
            _Py_VISIT_STACKREF(iframe->f_executable);
            int frame_result = _PyGC_VisitFrameStack(iframe, visit, arg);
            if (frame_result) {
                return frame_result;
            }
        }
    }
#endif
#if GREENLET_PY314 && defined(Py_GIL_DISABLED)
    // Visit the objects this greenlet's C-stack refs were holding when it
    // suspended (captured by capture_c_stack_refs). The free-threaded collector
    // only walks the running thread's _PyCStackRef list in
    // gc_visit_thread_stacks(), so without this a collection could free an
    // object reachable only through a suspended greenlet's C-stack ref and we'd
    // use it after free once the greenlet resumed. The snapshot is empty while
    // we're the running greenlet, so this is a no-op there.
    for (const OwnedObject& ref : this->c_stack_ref_snapshot) {
        Py_VISIT(ref.borrow());
    }
#endif
    // Note that we DO NOT visit ``delete_later``. Even if it's
    // non-null and we technically own a reference to it, its
    // reference count already went to 0 once and it was in the
    // process of being deallocated. The trash can mechanism linked it
    // into a list that will be cleaned at some later time, and it has
    // become untracked by the GC.
    return 0;
}

void PythonState::tp_clear(bool own_top_frame) noexcept
{
    PythonStateContext::tp_clear();
#if GREENLET_PY314 && defined(Py_GIL_DISABLED)
    this->c_stack_ref_snapshot.clear();
#endif
    // If we get here owning a frame,
    // we got dealloc'd without being finished. We may or may not be
    // in the same thread.
    if (own_top_frame) {
#if GREENLET_PY315
        // Release the references held by our suspended frames.
        // this->top_frame gets implicitly cleared by the Py_CLEAR(iframe->frame_obj)
        // of the first complete frame, so in the end we relinquish ownership of it.
        if (this->_top_frame) {
            for (_PyInterpreterFrame* iframe = this->_top_frame->f_frame;
                 iframe != nullptr; iframe = iframe->previous) {
                if (iframe->owner != FRAME_OWNED_BY_THREAD) {
                    continue;
                }
                // Clear the references held by this frame's evaluation stack.
                _PyStackRef* locals = iframe->localsplus;
                _PyStackRef* sp = iframe->stackpointer;
                if (sp) {
                    while (sp > locals) {
                        sp--;
                        PyStackRef_CLEAR(*sp);
                    }
                    iframe->stackpointer = locals;
                }
                Py_CLEAR(iframe->f_locals);
                Py_CLEAR(iframe->frame_obj);
                PyStackRef_CLEAR(iframe->f_funcobj);
                PyStackRef_CLEAR(iframe->f_executable);
            }
        }
        this->_top_frame.relinquish_ownership();
#else
        this->_top_frame.CLEAR();
#endif
    }
}

#if GREENLET_USE_CFRAME
void PythonState::set_new_cframe(_PyCFrame& frame) noexcept
{
    frame = *PyThreadState_GET()->cframe;
    /* Make the target greenlet refer to the stack value. */
    this->cframe = &frame;
    /*
      And restore the link to the previous frame so this one gets
      unliked appropriately.
    */
    this->cframe->previous = &PyThreadState_GET()->root_cframe;
}
#endif

const PythonState::OwnedFrame& PythonState::top_frame() const noexcept
{
    return this->_top_frame;
}

void PythonState::did_finish(PyThreadState* tstate) noexcept
{
#if GREENLET_PY311
    // See https://github.com/gevent/gevent/issues/1924 and
    // https://github.com/python-greenlet/greenlet/issues/328. In
    // short, Python 3.11 allocates memory for frames as a sort of
    // linked list that's kept as part of PyThreadState in the
    // ``datastack_chunk`` member and friends. These are saved and
    // restored as part of switching greenlets.
    //
    // When we initially switch to a greenlet, we set those to NULL.
    // That causes the frame management code to treat this like a
    // brand new thread and start a fresh list of chunks, beginning
    // with a new "root" chunk. As we make calls in this greenlet,
    // those chunks get added, and as calls return, they get popped.
    // But the frame code (pystate.c) is careful to make sure that the
    // root chunk never gets popped.
    //
    // Thus, when a greenlet exits for the last time, there will be at
    // least a single root chunk that we must be responsible for
    // deallocating.
    //
    // The complex part is that these chunks are allocated and freed
    // using ``_PyObject_VirtualAlloc``/``Free``. Those aren't public
    // functions, and they aren't exported for linking. It so happens
    // that we know they are just thin wrappers around the Arena
    // allocator, so we can use that directly to deallocate in a
    // compatible way.
    //
    // CAUTION: Check this implementation detail on every major version.
    //
    // It might be nice to be able to do this in our destructor, but
    // can we be sure that no one else is using that memory? Plus, as
    // described below, our pointers may not even be valid anymore. As
    // a special case, there is one time that we know we can do this,
    // and that's from the destructor of the associated UserGreenlet
    // (NOT main greenlet)
    PyObjectArenaAllocator alloc;
    _PyStackChunk* chunk = nullptr;
    if (tstate) {
        // We really did finish, we can never be switched to again.
        chunk = tstate->datastack_chunk;
        // Unfortunately, we can't do much sanity checking. Our
        // this->datastack_chunk pointer is out of date (evaluation may
        // have popped down through it already) so we can't verify that
        // we deallocate it. I don't think we can even check datastack_top
        // for the same reason.

        PyObject_GetArenaAllocator(&alloc);
        tstate->datastack_chunk = nullptr;
        tstate->datastack_limit = nullptr;
        tstate->datastack_top = nullptr;

    }
    else if (this->datastack_chunk) {
        // The UserGreenlet (NOT the main greenlet!) is being deallocated. If we're
        // still holding a stack chunk, it's garbage because we know
        // we can never switch back to let cPython clean it up.
        // Because the last time we got switched away from, and we
        // haven't run since then, we know our chain is valid and can
        // be dealloced.
        chunk = this->datastack_chunk;
        PyObject_GetArenaAllocator(&alloc);
    }

    if (alloc.free && chunk) {
        // In case the arena mechanism has been torn down already.
        while (chunk) {
            _PyStackChunk *prev = chunk->previous;
            chunk->previous = nullptr;
            alloc.free(alloc.ctx, chunk, chunk->size);
            chunk = prev;
        }
    }

    this->datastack_chunk = nullptr;
    this->datastack_limit = nullptr;
    this->datastack_top = nullptr;
#endif
}


}; // namespace greenlet

#endif // GREENLET_PYTHON_STATE_CPP
