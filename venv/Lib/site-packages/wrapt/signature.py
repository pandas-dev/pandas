"""The `with_signature` decorator: override the signature (and related
introspection attributes) of a wrapped callable without mutating the wrapped
function itself. Accepts a prototype callable, a prebuilt `inspect.Signature`,
or a factory that derives a signature from the wrapped function at decoration
time.

This replaces the need for the `adapter=` argument of `wrapt.decorator`,
which is planned for deprecation.
"""

from inspect import CO_VARARGS, CO_VARKEYWORDS, Parameter, Signature
from inspect import signature as _inspect_signature

from .__wrapt__ import (
    BaseObjectProxy,
    BoundFunctionWrapper,
    CallableObjectProxy,
    FunctionWrapper,
)

_MISSING = object()


def _cached(attr):
    """Read-only property that memoizes its derivation in a `_self_*` slot.

    `functools.cached_property` is unusable on ObjectProxy subclasses: wrapt
    overrides `__dict__` with a property returning the wrapped object's dict,
    so cached_property's `instance.__dict__[name] = value` would mutate the
    wrapped function. Storing under a `_self_*` name routes through wrapt's
    own setattr handling and lands in the proxy's real instance dict. Using
    a `@property` (data descriptor) also ensures we win over any value
    copied into the instance dict by `FunctionWrapper.__init__`.
    """
    slot = f"_self_cached_{attr.lstrip('_')}"

    def getter(self):
        value = getattr(self, slot, _MISSING)
        if value is _MISSING:
            value = self._derive(attr)
            self.__self_setattr__(slot, value)
        return value

    return property(getter)


def _derive_annotations(sig):
    ann = {
        p.name: p.annotation
        for p in sig.parameters.values()
        if p.annotation is not Parameter.empty
    }
    if sig.return_annotation is not Signature.empty:
        ann["return"] = sig.return_annotation
    return ann


def _derive_defaults(sig):
    defaults = tuple(
        p.default
        for p in sig.parameters.values()
        if p.kind in (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
        and p.default is not Parameter.empty
    )
    return defaults or None


def _derive_kwdefaults(sig):
    kwdefaults = {
        p.name: p.default
        for p in sig.parameters.values()
        if p.kind is Parameter.KEYWORD_ONLY and p.default is not Parameter.empty
    }
    return kwdefaults or None


def _derive_varnames(sig):
    pos_only, pos_or_kw, kw_only = [], [], []
    var_pos = var_kw = None
    for p in sig.parameters.values():
        if p.kind is Parameter.POSITIONAL_ONLY:
            pos_only.append(p.name)
        elif p.kind is Parameter.POSITIONAL_OR_KEYWORD:
            pos_or_kw.append(p.name)
        elif p.kind is Parameter.VAR_POSITIONAL:
            var_pos = p.name
        elif p.kind is Parameter.KEYWORD_ONLY:
            kw_only.append(p.name)
        elif p.kind is Parameter.VAR_KEYWORD:
            var_kw = p.name
    names = pos_only + pos_or_kw + kw_only
    if var_pos:
        names.append(var_pos)
    if var_kw:
        names.append(var_kw)
    return tuple(names)


class _SignatureCode(BaseObjectProxy):
    """Code-object proxy deriving argument-related co_* attrs from a Signature.

    Non-argument bits (co_flags for coroutine/generator, co_filename, etc.)
    fall through to the wrapped function's real __code__.
    """

    def __init__(self, wrapped, signature):
        super().__init__(wrapped)
        self._self_signature = signature

    def _derive(self, attr):
        sig = self._self_signature
        if attr == "co_argcount":
            kinds = (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
            return sum(p.kind in kinds for p in sig.parameters.values())
        if attr == "co_posonlyargcount":
            return sum(
                p.kind is Parameter.POSITIONAL_ONLY for p in sig.parameters.values()
            )
        if attr == "co_kwonlyargcount":
            return sum(
                p.kind is Parameter.KEYWORD_ONLY for p in sig.parameters.values()
            )
        if attr == "co_varnames":
            return _derive_varnames(sig)
        if attr == "co_flags":
            kinds = {p.kind for p in sig.parameters.values()}
            flags = self.__wrapped__.co_flags & ~(CO_VARARGS | CO_VARKEYWORDS)
            if Parameter.VAR_POSITIONAL in kinds:
                flags |= CO_VARARGS
            if Parameter.VAR_KEYWORD in kinds:
                flags |= CO_VARKEYWORDS
            return flags
        raise AttributeError(attr)

    co_argcount = _cached("co_argcount")
    co_posonlyargcount = _cached("co_posonlyargcount")
    co_kwonlyargcount = _cached("co_kwonlyargcount")
    co_varnames = _cached("co_varnames")
    co_flags = _cached("co_flags")


class _SignatureMixin:
    """Shared Signature-derived introspection attrs for wrapper + surrogate."""

    def _derive(self, attr):
        sig = self._self_signature
        if attr == "__annotations__":
            return _derive_annotations(sig)
        if attr == "__defaults__":
            return _derive_defaults(sig)
        if attr == "__kwdefaults__":
            return _derive_kwdefaults(sig)
        if attr == "__code__":
            return _SignatureCode(self.__wrapped__.__code__, sig)
        raise AttributeError(attr)

    @property
    def __signature__(self):
        return self._self_signature

    __annotations__ = _cached("__annotations__")
    __defaults__ = _cached("__defaults__")
    __kwdefaults__ = _cached("__kwdefaults__")
    __code__ = _cached("__code__")


class _SignatureFunctionSurrogate(_SignatureMixin, CallableObjectProxy):
    """Function surrogate exposing Signature-derived introspection attrs.

    Used as __func__ of a bound wrapper so that inspect.signature -- which
    treats the bound wrapper as a MethodType and consults __func__ -- sees
    our override and strips self/cls correctly.
    """

    def __init__(self, wrapped, signature):
        super().__init__(wrapped)
        self._self_signature = signature


class _BoundSignatureFunctionWrapper(BoundFunctionWrapper):
    @property
    def __func__(self):
        return _SignatureFunctionSurrogate(
            self.__wrapped__.__func__, self._self_parent._self_signature
        )

    @property
    def __signature__(self):
        return self._self_parent._self_signature


class _SignatureFunctionWrapper(_SignatureMixin, FunctionWrapper):
    __bound_function_wrapper__ = _BoundSignatureFunctionWrapper

    def __init__(self, wrapped, wrapper, signature):
        super().__init__(wrapped, wrapper)
        self._self_signature = signature

    @property
    def __func__(self):
        return _SignatureFunctionSurrogate(self.__wrapped__, self._self_signature)


def with_signature(wrapped=None, /, *, prototype=None, signature=None, factory=None):
    """Override the signature of a wrapped callable.

    Exactly one of `prototype`, `signature`, or `factory` must be supplied:

    - `prototype`: a callable whose signature will be used.
    - `signature`: a prebuilt `inspect.Signature` object.
    - `factory`: a callable `factory(wrapped)` invoked at decoration time
      that returns either a `Signature` or a prototype callable.

    The resulting wrapper exposes the override via `__signature__`, and
    derives `__annotations__`, `__defaults__`, `__kwdefaults__`, and the
    argument-related attributes of `__code__` from the same source. The
    wrapped function is not mutated. Calling behaviour is unchanged.
    """

    specified = sum(x is not None for x in (prototype, signature, factory))
    if specified == 0:
        raise TypeError(
            "with_signature requires one of prototype=, signature=, or factory="
        )
    if specified > 1:
        raise TypeError(
            "with_signature accepts only one of prototype=, signature=, or factory="
        )

    def _decorator(wrapped):
        def _wrapper(wrapped, instance, args, kwargs):
            return wrapped(*args, **kwargs)

        if signature is not None:
            resolved = signature
        elif prototype is not None:
            resolved = _inspect_signature(prototype)
        else:
            produced = factory(wrapped)
            resolved = (
                produced
                if isinstance(produced, Signature)
                else _inspect_signature(produced)
            )

        return _SignatureFunctionWrapper(wrapped, _wrapper, resolved)

    if wrapped is None:
        return _decorator
    return _decorator(wrapped)
