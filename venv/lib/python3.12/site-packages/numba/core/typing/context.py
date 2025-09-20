from collections import defaultdict
from collections.abc import Sequence
import typing as _tp
import types as pytypes
import weakref
import threading
import contextlib
import operator

from numba.core import types, errors, config
from numba.core.typeconv import Conversion, rules
from numba.core.typing import templates
from numba.core.utils import order_by_target_specificity
from .typeof import typeof, Purpose

from numba.core import utils


class Rating(object):
    __slots__ = 'promote', 'safe_convert', "unsafe_convert"

    def __init__(self):
        self.promote = 0
        self.safe_convert = 0
        self.unsafe_convert = 0

    def astuple(self):
        """Returns a tuple suitable for comparing with the worse situation
        start first.
        """
        return (self.unsafe_convert, self.safe_convert, self.promote)

    def __add__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        rsum = Rating()
        rsum.promote = self.promote + other.promote
        rsum.safe_convert = self.safe_convert + other.safe_convert
        rsum.unsafe_convert = self.unsafe_convert + other.unsafe_convert
        return rsum


class CallStack(Sequence):
    """
    A compile-time call stack
    """

    def __init__(self):
        self._stack = []
        self._lock = threading.RLock()
        # fail_cache only last for the current compilation session
        self._fail_cache = {}

    def __getitem__(self, index):
        """
        Returns item in the stack where index=0 is the top and index=1 is
        the second item from the top.
        """
        return self._stack[len(self) - index - 1]

    def __len__(self):
        return len(self._stack)

    @contextlib.contextmanager
    def register(self, target, typeinfer, func_id, args):
        with contextlib.ExitStack() as undo:
            # guard compiling the same function with the same signature
            if self.match(func_id.func, args):
                msg = "compiler re-entrant to the same function signature"
                raise errors.NumbaRuntimeError(msg)

            # Acquire lock
            undo.enter_context(self._lock)

            # Clear fail_cache at the start and end of a compilation session
            def clear_fail_cache(*exc):
                if config.DISABLE_TYPEINFER_FAIL_CACHE:
                    return   # bypass
                # Clear cache if stack is empty?
                if not self._stack:
                    self._fail_cache.clear()

            clear_fail_cache()
            undo.push(clear_fail_cache)

            # Setup callframe
            self._stack.append(CallFrame(target, typeinfer, func_id, args))

            @undo.push
            def undo_stack(*exc):
                self._stack.pop()

            yield

    def finditer(self, py_func):
        """
        Yields frame that matches the function object starting from the top
        of stack.
        """
        for frame in self:
            if frame.func_id.func is py_func:
                yield frame

    def findfirst(self, py_func):
        """
        Returns the first result from `.finditer(py_func)`; or None if no match.
        """
        try:
            return next(self.finditer(py_func))
        except StopIteration:
            return

    def match(self, py_func, args):
        """
        Returns first function that matches *py_func* and the arguments types in
        *args*; or, None if no match.
        """
        for frame in self.finditer(py_func):
            if frame.args == args:
                return frame

    def lookup_resolve_cache(self, func, args, kws) -> "_ResolveCache":
        """Lookup resolution cache for the given function type and argument
        types.
        """
        if not self._stack or config.DISABLE_TYPEINFER_FAIL_CACHE:
            # if callstack is empty, bypass fail_cache
            return _ResolveCache()

        def normalize_dict(obj):
            if isinstance(obj, dict):
                return tuple(sorted(kws.items()))
            return kws

        def hashable(obj):
            try:
                hash(obj)
            except TypeError:
                return False
            else:
                return True

        key = func, args, normalize_dict(kws)
        if not hashable(key):
            return _ResolveCache()
        return self._fail_cache.setdefault(key, _ResolveCache())


class _ResolveCache(object):
    """
    A cache for function resolution result.
    Currently only remember failed attempts.
    """
    _status: str
    _exc: _tp.Optional[BaseException]

    def __init__(self):
        self._status = "unmarked"
        self._exc = None

    def mark_error(self, exc) -> None:
        """Mark the function resolution as failed with an exception."""
        self._status = "error"
        self._exc = exc

    def mark_failed(self) -> None:
        """Mark the function resolution as failed."""
        self._status = "failed"

    def replay_failure(self) -> None:
        """Replay the failure if it has been marked as failed or error."""
        if self._status == "error":
            raise self._exc
        else:
            assert self._status == "failed"
            return None

    def has_failed_previously(self) -> bool:
        """Return True if the function resolution has failed previously."""
        return self._status in {"failed", "error"}


class CallFrame(object):
    """
    A compile-time call frame
    """
    def __init__(self, target, typeinfer, func_id, args):
        self.typeinfer = typeinfer
        self.func_id = func_id
        self.args = args
        self.target = target
        self._inferred_retty = set()

    def __repr__(self):
        return "CallFrame({}, {})".format(self.func_id, self.args)

    def add_return_type(self, return_type):
        """Add *return_type* to the list of inferred return-types.
        If there are too many, raise `TypingError`.
        """
        # The maximum limit is picked arbitrarily.
        # Don't think that this needs to be user configurable.
        RETTY_LIMIT = 16
        self._inferred_retty.add(return_type)
        if len(self._inferred_retty) >= RETTY_LIMIT:
            m = "Return type of recursive function does not converge"
            raise errors.TypingError(m)


class BaseContext(object):
    """A typing context for storing function typing constrain template.
    """

    def __init__(self):
        # A list of installed registries
        self._registries = {}
        # Typing declarations extracted from the registries or other sources
        self._functions = defaultdict(list)
        self._attributes = defaultdict(list)
        self._globals = utils.UniqueDict()
        self.tm = rules.default_type_manager
        self.callstack = CallStack()

        # Initialize
        self.init()

    def init(self):
        """
        Initialize the typing context.  Can be overridden by subclasses.
        """

    def refresh(self):
        """
        Refresh context with new declarations from known registries.
        Useful for third-party extensions.
        """
        self.load_additional_registries()
        # Some extensions may have augmented the builtin registry
        self._load_builtins()

    def explain_function_type(self, func):
        """
        Returns a string description of the type of a function
        """
        desc = []
        defns = []
        param = False
        if isinstance(func, types.Callable):
            sigs, param = func.get_call_signatures()
            defns.extend(sigs)

        elif func in self._functions:
            for tpl in self._functions[func]:
                param = param or hasattr(tpl, 'generic')
                defns.extend(getattr(tpl, 'cases', []))

        else:
            msg = "No type info available for {func!r} as a callable."
            desc.append(msg.format(func=func))

        if defns:
            desc = ['Known signatures:']
            for sig in defns:
                desc.append(' * {0}'.format(sig))

        return '\n'.join(desc)

    def resolve_function_type(self, func, args, kws):
        """
        Resolve function type *func* for argument types *args* and *kws*.
        A signature is returned.
        """
        cache = self.callstack.lookup_resolve_cache(func, args, kws)
        if cache.has_failed_previously():
            return cache.replay_failure()

        # Prefer user definition first
        try:
            res = self._resolve_user_function_type(func, args, kws)
        except errors.TypingError as e:
            # Capture any typing error
            last_exception = e
            res = None
        else:
            last_exception = None

        # Return early we know there's a working user function
        if res is not None:
            return res

        # Check builtin functions
        res = self._resolve_builtin_function_type(func, args, kws)

        # Re-raise last_exception if no function type has been found
        if res is None and last_exception is not None:
            cache.mark_error(last_exception)
            raise last_exception

        if res is None:
            cache.mark_failed()

        return res

    def _resolve_builtin_function_type(self, func, args, kws):
        # NOTE: we should reduce usage of this
        if func in self._functions:
            # Note: Duplicating code with types.Function.get_call_type().
            #       *defns* are CallTemplates.
            defns = self._functions[func]
            for defn in defns:
                for support_literals in [True, False]:
                    if support_literals:
                        res = defn.apply(args, kws)
                    else:
                        fixedargs = [types.unliteral(a) for a in args]
                        res = defn.apply(fixedargs, kws)
                    if res is not None:
                        return res

    def _resolve_user_function_type(self, func, args, kws, literals=None):
        # It's not a known function type, perhaps it's a global?
        functy = self._lookup_global(func)
        if functy is not None:
            func = functy

        if isinstance(func, types.Type):
            # If it's a type, it may support a __call__ method
            func_type = self.resolve_getattr(func, "__call__")
            if func_type is not None:
                # The function has a __call__ method, type its call.
                return self.resolve_function_type(func_type, args, kws)

        if isinstance(func, types.Callable):
            # XXX fold this into the __call__ attribute logic?
            return func.get_call_type(self, args, kws)

    def _get_attribute_templates(self, typ):
        """
        Get matching AttributeTemplates for the Numba type.
        """
        if typ in self._attributes:
            for attrinfo in self._attributes[typ]:
                yield attrinfo
        else:
            for cls in type(typ).__mro__:
                if cls in self._attributes:
                    for attrinfo in self._attributes[cls]:
                        yield attrinfo

    def resolve_getattr(self, typ, attr):
        """
        Resolve getting the attribute *attr* (a string) on the Numba type.
        The attribute's type is returned, or None if resolution failed.
        """
        def core(typ):
            out = self.find_matching_getattr_template(typ, attr)
            if out:
                return out['return_type']

        out = core(typ)
        if out is not None:
            return out

        # Try again without literals
        out = core(types.unliteral(typ))
        if out is not None:
            return out

        if isinstance(typ, types.Module):
            attrty = self.resolve_module_constants(typ, attr)
            if attrty is not None:
                return attrty

    def find_matching_getattr_template(self, typ, attr):

        templates = list(self._get_attribute_templates(typ))

        # get the order in which to try templates
        from numba.core.target_extension import get_local_target # circular
        target_hw = get_local_target(self)
        order = order_by_target_specificity(target_hw, templates, fnkey=attr)

        for template in order:
            return_type = template.resolve(typ, attr)
            if return_type is not None:
                return {
                    'template': template,
                    'return_type': return_type,
                }

    def resolve_setattr(self, target, attr, value):
        """
        Resolve setting the attribute *attr* (a string) on the *target* type
        to the given *value* type.
        A function signature is returned, or None if resolution failed.
        """
        for attrinfo in self._get_attribute_templates(target):
            expectedty = attrinfo.resolve(target, attr)
            # NOTE: convertibility from *value* to *expectedty* is left to
            # the caller.
            if expectedty is not None:
                return templates.signature(types.void, target, expectedty)

    def resolve_static_getitem(self, value, index):
        assert not isinstance(index, types.Type), index
        args = value, index
        kws = ()
        return self.resolve_function_type("static_getitem", args, kws)

    def resolve_static_setitem(self, target, index, value):
        assert not isinstance(index, types.Type), index
        args = target, index, value
        kws = {}
        return self.resolve_function_type("static_setitem", args, kws)

    def resolve_setitem(self, target, index, value):
        assert isinstance(index, types.Type), index
        fnty = self.resolve_value_type(operator.setitem)
        sig = fnty.get_call_type(self, (target, index, value), {})
        return sig

    def resolve_delitem(self, target, index):
        args = target, index
        kws = {}
        fnty = self.resolve_value_type(operator.delitem)
        sig = fnty.get_call_type(self, args, kws)
        return sig

    def resolve_module_constants(self, typ, attr):
        """
        Resolve module-level global constants.
        Return None or the attribute type
        """
        assert isinstance(typ, types.Module)
        attrval = getattr(typ.pymod, attr)
        try:
            return self.resolve_value_type(attrval)
        except ValueError:
            pass

    def resolve_value_type(self, val):
        """
        Return the numba type of a Python value that is being used
        as a runtime constant.
        ValueError is raised for unsupported types.
        """
        try:
            ty = typeof(val, Purpose.constant)
        except ValueError as e:
            # Make sure the exception doesn't hold a reference to the user
            # value.
            typeof_exc = utils.erase_traceback(e)
        else:
            return ty

        if isinstance(val, types.ExternalFunction):
            return val

        # Try to look up target specific typing information
        ty = self._get_global_type(val)
        if ty is not None:
            return ty

        raise typeof_exc

    def resolve_value_type_prefer_literal(self, value):
        """Resolve value type and prefer Literal types whenever possible.
        """
        lit = types.maybe_literal(value)
        if lit is None:
            return self.resolve_value_type(value)
        else:
            return lit

    def _get_global_type(self, gv):
        ty = self._lookup_global(gv)
        if ty is not None:
            return ty
        if isinstance(gv, pytypes.ModuleType):
            return types.Module(gv)

    def _load_builtins(self):
        # Initialize declarations
        from numba.core.typing import builtins, arraydecl, npdatetime  # noqa: F401, E501
        from numba.core.typing import ctypes_utils, bufproto           # noqa: F401, E501
        from numba.core.unsafe import eh                    # noqa: F401

        self.install_registry(templates.builtin_registry)

    def load_additional_registries(self):
        """
        Load target-specific registries.  Can be overridden by subclasses.
        """

    def install_registry(self, registry):
        """
        Install a *registry* (a templates.Registry instance) of function,
        attribute and global declarations.
        """
        try:
            loader = self._registries[registry]
        except KeyError:
            loader = templates.RegistryLoader(registry)
            self._registries[registry] = loader

        from numba.core.target_extension import (get_local_target,
                                                 resolve_target_str)
        current_target = get_local_target(self)

        def is_for_this_target(ftcls):
            metadata = getattr(ftcls, 'metadata', None)
            if metadata is None:
                return True

            target_str = metadata.get('target')
            if target_str is None:
                return True

            # There may be pending registrations for nonexistent targets.
            # Ideally it would be impossible to leave a registration pending
            # for an invalid target, but in practice this is exceedingly
            # difficult to guard against - many things are registered at import
            # time, and eagerly reporting an error when registering for invalid
            # targets would require that all target registration code is
            # executed prior to all typing registrations during the import
            # process; attempting to enforce this would impose constraints on
            # execution order during import that would be very difficult to
            # resolve and maintain in the presence of typical code maintenance.
            # Furthermore, these constraints would be imposed not only on
            # Numba internals, but also on its dependents.
            #
            # Instead of that enforcement, we simply catch any occurrences of
            # registrations for targets that don't exist, and report that
            # they're not for this target. They will then not be encountered
            # again during future typing context refreshes (because the
            # loader's new registrations are a stream_list that doesn't yield
            # previously-yielded items).
            try:
                ft_target = resolve_target_str(target_str)
            except errors.NonexistentTargetError:
                return False

            return current_target.inherits_from(ft_target)

        for ftcls in loader.new_registrations('functions'):
            if not is_for_this_target(ftcls):
                continue
            self.insert_function(ftcls(self))
        for ftcls in loader.new_registrations('attributes'):
            if not is_for_this_target(ftcls):
                continue
            self.insert_attributes(ftcls(self))
        for gv, gty in loader.new_registrations('globals'):
            existing = self._lookup_global(gv)
            if existing is None:
                self.insert_global(gv, gty)
            else:
                # A type was already inserted, see if we can add to it
                newty = existing.augment(gty)
                if newty is None:
                    raise TypeError("cannot augment %s with %s"
                                    % (existing, gty))
                self._remove_global(gv)
                self._insert_global(gv, newty)

    def _lookup_global(self, gv):
        """
        Look up the registered type for global value *gv*.
        """
        try:
            gv = weakref.ref(gv)
        except TypeError:
            pass
        try:
            return self._globals.get(gv, None)
        except TypeError:
            # Unhashable type
            return None

    def _insert_global(self, gv, gty):
        """
        Register type *gty* for value *gv*.  Only a weak reference
        to *gv* is kept, if possible.
        """
        def on_disposal(wr, pop=self._globals.pop):
            # pop() is pre-looked up to avoid a crash late at shutdown on 3.5
            # (https://bugs.python.org/issue25217)
            pop(wr)
        try:
            gv = weakref.ref(gv, on_disposal)
        except TypeError:
            pass
        self._globals[gv] = gty

    def _remove_global(self, gv):
        """
        Remove the registered type for global value *gv*.
        """
        try:
            gv = weakref.ref(gv)
        except TypeError:
            pass
        del self._globals[gv]

    def insert_global(self, gv, gty):
        self._insert_global(gv, gty)

    def insert_attributes(self, at):
        key = at.key
        self._attributes[key].append(at)

    def insert_function(self, ft):
        key = ft.key
        self._functions[key].append(ft)

    def insert_user_function(self, fn, ft):
        """Insert a user function.

        Args
        ----
        - fn:
            object used as callee
        - ft:
            function template
        """
        self._insert_global(fn, types.Function(ft))

    def can_convert(self, fromty, toty):
        """
        Check whether conversion is possible from *fromty* to *toty*.
        If successful, return a numba.typeconv.Conversion instance;
        otherwise None is returned.
        """
        if fromty == toty:
            return Conversion.exact
        else:
            # First check with the type manager (some rules are registered
            # at startup there, see numba.typeconv.rules)
            conv = self.tm.check_compatible(fromty, toty)
            if conv is not None:
                return conv

            # Fall back on type-specific rules
            forward = fromty.can_convert_to(self, toty)
            backward = toty.can_convert_from(self, fromty)
            if backward is None:
                return forward
            elif forward is None:
                return backward
            else:
                return min(forward, backward)

    def _rate_arguments(self, actualargs, formalargs, unsafe_casting=True,
                        exact_match_required=False):
        """
        Rate the actual arguments for compatibility against the formal
        arguments.  A Rating instance is returned, or None if incompatible.
        """
        if len(actualargs) != len(formalargs):
            return None
        rate = Rating()
        for actual, formal in zip(actualargs, formalargs):
            conv = self.can_convert(actual, formal)
            if conv is None:
                return None
            elif not unsafe_casting and conv >= Conversion.unsafe:
                return None
            elif exact_match_required and conv != Conversion.exact:
                return None

            if conv == Conversion.promote:
                rate.promote += 1
            elif conv == Conversion.safe:
                rate.safe_convert += 1
            elif conv == Conversion.unsafe:
                rate.unsafe_convert += 1
            elif conv == Conversion.exact:
                pass
            else:
                raise AssertionError("unreachable", conv)

        return rate

    def install_possible_conversions(self, actualargs, formalargs):
        """
        Install possible conversions from the actual argument types to
        the formal argument types in the C++ type manager.
        Return True if all arguments can be converted.
        """
        if len(actualargs) != len(formalargs):
            return False
        for actual, formal in zip(actualargs, formalargs):
            if self.tm.check_compatible(actual, formal) is not None:
                # This conversion is already known
                continue
            conv = self.can_convert(actual, formal)
            if conv is None:
                return False
            assert conv is not Conversion.exact
            self.tm.set_compatible(actual, formal, conv)
        return True

    def resolve_overload(self, key, cases, args, kws,
                         allow_ambiguous=True, unsafe_casting=True,
                         exact_match_required=False):
        """
        Given actual *args* and *kws*, find the best matching
        signature in *cases*, or None if none matches.
        *key* is used for error reporting purposes.
        If *allow_ambiguous* is False, a tie in the best matches
        will raise an error.
        If *unsafe_casting* is False, unsafe casting is forbidden.
        """
        assert not kws, "Keyword arguments are not supported, yet"
        options = {
            'unsafe_casting': unsafe_casting,
            'exact_match_required': exact_match_required,
        }
        # Rate each case
        candidates = []
        for case in cases:
            if len(args) == len(case.args):
                rating = self._rate_arguments(args, case.args, **options)
                if rating is not None:
                    candidates.append((rating.astuple(), case))

        # Find the best case
        candidates.sort(key=lambda i: i[0])
        if candidates:
            best_rate, best = candidates[0]
            if not allow_ambiguous:
                # Find whether there is a tie and if so, raise an error
                tied = []
                for rate, case in candidates:
                    if rate != best_rate:
                        break
                    tied.append(case)
                if len(tied) > 1:
                    args = (key, args, '\n'.join(map(str, tied)))
                    msg = "Ambiguous overloading for %s %s:\n%s" % args
                    raise TypeError(msg)
            # Simply return the best matching candidate in order.
            # If there is a tie, since list.sort() is stable, the first case
            # in the original order is returned.
            # (this can happen if e.g. a function template exposes
            #  (int32, int32) -> int32 and (int64, int64) -> int64,
            #  and you call it with (int16, int16) arguments)
            return best

    def unify_types(self, *typelist):
        # Sort the type list according to bit width before doing
        # pairwise unification (with thanks to aterrel).
        def keyfunc(obj):
            """Uses bitwidth to order numeric-types.
            Fallback to stable, deterministic sort.
            """
            return getattr(obj, 'bitwidth', 0)
        typelist = sorted(typelist, key=keyfunc)
        unified = typelist[0]
        for tp in typelist[1:]:
            unified = self.unify_pairs(unified, tp)
            if unified is None:
                break
        return unified

    def unify_pairs(self, first, second):
        """
        Try to unify the two given types.  A third type is returned,
        or None in case of failure.
        """
        if first == second:
            return first

        if first is types.undefined:
            return second
        elif second is types.undefined:
            return first

        # Types with special unification rules
        unified = first.unify(self, second)
        if unified is not None:
            return unified

        unified = second.unify(self, first)
        if unified is not None:
            return unified

        # Other types with simple conversion rules
        conv = self.can_convert(fromty=first, toty=second)
        if conv is not None and conv <= Conversion.safe:
            # Can convert from first to second
            return second

        conv = self.can_convert(fromty=second, toty=first)
        if conv is not None and conv <= Conversion.safe:
            # Can convert from second to first
            return first

        if isinstance(first, types.Literal) or \
           isinstance(second, types.Literal):
            first = types.unliteral(first)
            second = types.unliteral(second)
            return self.unify_pairs(first, second)

        # Cannot unify
        return None


class Context(BaseContext):

    def load_additional_registries(self):
        from . import (
            cffi_utils,
            cmathdecl,
            enumdecl,
            listdecl,
            mathdecl,
            npydecl,
            setdecl,
            dictdecl,
        )
        self.install_registry(cffi_utils.registry)
        self.install_registry(cmathdecl.registry)
        self.install_registry(enumdecl.registry)
        self.install_registry(listdecl.registry)
        self.install_registry(mathdecl.registry)
        self.install_registry(npydecl.registry)
        self.install_registry(setdecl.registry)
        self.install_registry(dictdecl.registry)
