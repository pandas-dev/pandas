from __future__ import annotations

import contextlib
import inspect
import logging
import math
import os
import secrets
import sys
import time
import warnings
from abc import ABCMeta, abstractmethod
from collections.abc import Callable, Hashable
from contextlib import contextmanager
from dataclasses import dataclass
from itertools import count, starmap
from threading import Condition, RLock, get_ident, local
from typing import TYPE_CHECKING, Final, Literal, NoReturn, TypedDict, TypeVar, cast
from weakref import WeakKeyDictionary, WeakValueDictionary

from ._error import SoftFileLockLifetimeWarning, Timeout
from ._util import break_lock_file

#: No explicit file permission mode was passed. Lock files then open with 0o666 so umask and default ACLs pick
#: the final permissions, and fchmod is skipped to preserve POSIX default ACL inheritance.
_UNSET_FILE_MODE: Final[int] = -1

#: Ceiling on the retry counter used as a power of two, so a long contended wait cannot overflow the backoff multiply.
_MAX_BACKOFF_EXPONENT: Final[int] = 20

#: How a context manager reconciles a body failure with a release failure on exit (see the property of this name).
ContextErrorPolicy = Literal["chain", "group"]
_CONTEXT_ERROR_POLICIES: Final[frozenset[str]] = frozenset({"chain", "group"})

#: What a descriptor-owning backend does with an ``os.close`` failure after relinquishing ownership (see the property).
CloseErrorPolicy = Literal["default", "raise", "suppress"]
_CLOSE_ERROR_POLICIES: Final[frozenset[str]] = frozenset({"default", "raise", "suppress"})

if TYPE_CHECKING:
    from collections.abc import Generator
    from types import TracebackType
    from typing import Protocol

    from ._read_write import ReadWriteLock
    from ._soft_rw import SoftReadWriteLock

    class _ForkResettable(Protocol):
        def _reset_after_fork_in_child(self) -> None: ...

    class _ForkDescriptorOwner(Protocol):
        def _descriptors_for_fork(self) -> tuple[tuple[int, tuple[int, int] | None], ...]: ...

    # Matched against the class object itself rather than `type[...]` of it. A metaclass supplies this method to the
    # class while leaving instances without it, so a `type[_ForkResettableClass]` bound rejects `ReadWriteLock`.
    class _ForkResettableClass(Protocol):
        def _reset_class_after_fork(self) -> None: ...

    class _RegisterAtFork(Protocol):
        def __call__(
            self,
            *,
            before: Callable[[], None] | None = None,
            after_in_parent: Callable[[], None] | None = None,
            after_in_child: Callable[[], None] | None = None,
        ) -> None: ...

    if sys.version_info >= (3, 11):  # pragma: no cover (py311+)
        from typing import Self
    else:  # pragma: no cover (<py311)
        from typing_extensions import Self

_LOGGER: Final[logging.Logger] = logging.getLogger("filelock")
_REGISTER_AT_FORK: Final[_RegisterAtFork | None] = cast("_RegisterAtFork | None", getattr(os, "register_at_fork", None))
_HAS_REGISTER_AT_FORK: Final[bool] = _REGISTER_AT_FORK is not None

_ExtraValue = TypeVar("_ExtraValue")
_MarkerValue = TypeVar("_MarkerValue")
_SubclassValue = TypeVar("_SubclassValue")
_LockInitValue = float | int | bool | str | None | Callable[[int], None]


class LockOptions(TypedDict, total=False):
    """Every option the metaclass forwards, so a subclass adding its own can still type what it passes through."""

    timeout: float
    mode: int
    thread_local: bool
    blocking: bool
    is_singleton: bool
    poll_interval: float
    lifetime: float | None
    context_error_policy: ContextErrorPolicy
    close_error_policy: CloseErrorPolicy
    fallback_to_soft: bool
    preserve_lock_file: bool
    on_acquired: Callable[[int], None] | None


def _exception_group_cls() -> type[BaseException]:
    # BaseExceptionGroup is a builtin on 3.11+; on 3.10 it needs the exceptiongroup backport. filelock keeps zero
    # runtime dependencies, so the backport is imported lazily rather than required, and only group mode needs it.
    if sys.version_info >= (3, 11):  # pragma: no cover (py311+)
        return BaseExceptionGroup  # ruff:ignore[undefined-name]  # builtin on 3.11+
    # Alias the import so BaseExceptionGroup above stays the builtin rather than an unbound local of this function.
    from exceptiongroup import (  # ruff:ignore[import-outside-top-level]  # pragma: no cover (<py311)
        BaseExceptionGroup as _Backport,
    )

    return _Backport  # pragma: no cover (<py311)


def _raise_grouped_errors(
    message: str,
    first_error: BaseException,
    second_error: BaseException,
    *additional_errors: BaseException,
    marker: tuple[str, _MarkerValue] | None = None,
) -> NoReturn:
    errors = (first_error, second_error, *additional_errors)
    _detach_grouped_contexts(errors)
    group = _exception_group_cls()(message, errors)
    if marker is not None:
        setattr(group, marker[0], marker[1])
    raise group from None


def _detach_grouped_contexts(errors: tuple[BaseException, ...]) -> None:
    seen: set[int] = set()
    pending = list(errors)
    while pending:
        error = pending.pop()
        if id(error) in seen:
            continue
        seen.add(id(error))
        if (context := error.__context__) is not None and (
            context is error
            or _same_exception_tree(error, context)
            or any(context is root or _contains_exception(root, context) for root in errors)
        ):
            error.__context__ = None
        elif context is not None:
            pending.append(context)
        if error.__cause__ is not None:
            pending.append(error.__cause__)
        if isinstance(error, _exception_group_cls()):
            pending.extend(cast("_ExceptionGroupProtocol", error).exceptions)


def _same_exception_tree(first: BaseException, second: BaseException) -> bool:
    pending = [(first, second)]
    seen: set[tuple[int, int]] = set()
    while pending:
        first_error, second_error = pending.pop()
        if first_error is second_error:
            continue
        if (pair := (id(first_error), id(second_error))) in seen:
            continue
        seen.add(pair)
        if (
            type(first_error) is not type(second_error)
            or not isinstance(first_error, _exception_group_cls())
            or not isinstance(second_error, _exception_group_cls())
        ):
            return False
        first_group = cast("_ExceptionGroupProtocol", first_error)
        second_group = cast("_ExceptionGroupProtocol", second_error)
        if first_group.message != second_group.message or len(first_group.exceptions) != len(second_group.exceptions):
            return False
        pending.extend(zip(first_group.exceptions, second_group.exceptions, strict=True))
    return True


def _contains_exception(error: BaseException, target: BaseException | None) -> bool:
    if target is None or not isinstance(error, _exception_group_cls()):
        return False
    pending = list(cast("_ExceptionGroupProtocol", error).exceptions)
    seen: set[int] = set()
    while pending:
        child = pending.pop()
        if child is target:
            return True
        if id(child) in seen:
            continue
        seen.add(id(child))
        if isinstance(child, _exception_group_cls()):
            pending.extend(cast("_ExceptionGroupProtocol", child).exceptions)
    return False


def _append_exception_context(error: BaseException, context: BaseException) -> None:
    if _exception_graph_contains(error, context) or _exception_graph_contains(context, error):
        return
    if error.__context__ is None:
        error.__context__ = context
        return
    tail = error
    seen: set[int] = set()
    while id(tail) not in seen:
        seen.add(id(tail))
        if (next_error := tail.__cause__ if tail.__cause__ is not None else tail.__context__) is None:
            tail.__context__ = context
            return
        tail = next_error


def _exception_graph_contains(error: BaseException, target: BaseException) -> bool:
    pending = [error]
    seen: set[int] = set()
    while pending:
        current = pending.pop()
        if current is target:
            return True
        if id(current) in seen:  # pragma: no cover - arbitrary caller exceptions can contain cycles
            continue
        seen.add(id(current))
        if current.__cause__ is not None:
            pending.append(current.__cause__)
        if current.__context__ is not None:
            pending.append(current.__context__)
        if isinstance(current, _exception_group_cls()):
            pending.extend(cast("_ExceptionGroupProtocol", current).exceptions)
    return False


def _grouped_errors(
    error: BaseException, message: str, marker: tuple[str, _MarkerValue]
) -> tuple[BaseException, ...] | None:
    if not isinstance(error, _exception_group_cls()):
        return None
    group = cast("_ExceptionGroupProtocol", error)
    return group.exceptions if group.message == message and getattr(group, marker[0], None) is marker[1] else None


if TYPE_CHECKING:

    class _ExceptionGroupProtocol(Protocol):
        @property
        def message(self) -> str: ...

        @property
        def exceptions(self) -> tuple[BaseException, ...]: ...


def _raise_chained_errors(first_error: BaseException, second_error: BaseException | None = None) -> NoReturn:
    if second_error is None:
        first_context = first_error.__context__
        try:
            raise first_error  # ruff:ignore[raise-within-try]  # the handler restores caller-supplied context before propagation
        except BaseException:
            first_error.__context__ = first_context
            raise
    if (second_context := second_error.__context__) is not None and second_context is not first_error:
        _detach_exception_context(second_context, first_error)
        _append_exception_context(first_error, second_context)
    first_context = first_error.__context__
    try:
        raise first_error  # ruff:ignore[raise-within-try]  # the second raise needs this error as implicit context
    except BaseException:  # ruff:ignore[blind-except]  # first_error may be a control-flow exception
        first_error.__context__ = first_context
        try:
            raise second_error  # ruff:ignore[raise-within-try]  # the handler makes the chain interpreter-independent
        except BaseException:
            second_error.__context__ = first_error
            first_error.__context__ = first_context
            raise


def _detach_exception_context(error: BaseException, target: BaseException) -> None:
    pending = [error]
    seen: set[int] = set()
    while pending:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        if current.__context__ is target:
            current.__context__ = None
        elif current.__context__ is not None:
            pending.append(current.__context__)
        if current.__cause__ is not None:
            pending.append(current.__cause__)
        if isinstance(current, _exception_group_cls()):
            pending.extend(cast("_ExceptionGroupProtocol", current).exceptions)


def _raise_body_and_release(body_error: BaseException, release_error: BaseException) -> NoReturn:
    # Group mode: surface the body failure and the release failure as sibling leaves instead of letting one hide in the
    # other's __context__. BaseExceptionGroup returns a plain ExceptionGroup when both leaves subclass Exception, so
    # ``except*`` and ``except Exception`` still catch them; a BaseException leaf (KeyboardInterrupt, CancelledError)
    # keeps the group outside ordinary handlers. ``from None`` stops the group itself gaining a redundant __context__.
    _raise_grouped_errors("lock body and release both failed", body_error, release_error)


def _raise_cleanup_errors(
    message: str,
    primary_error: BaseException,
    *cleanup_errors: BaseException | None,
) -> NoReturn:
    _raise_grouped_errors(
        message,
        primary_error,
        *(error for error in cleanup_errors if error is not None),
    )


# On Windows os.path.realpath calls CreateFileW with share_mode=0, which blocks concurrent DeleteFileW and causes
# livelocks under threaded contention with SoftFileLock. os.path.abspath is purely string-based and avoids this.
_resolve_dir: Final[Callable[[str], str]] = os.path.abspath if sys.platform == "win32" else os.path.realpath


def _canonical(path: str | os.PathLike[str]) -> str:
    """
    Return one stable key for *path*, collapsing equivalent spellings without following a final symlink.

    Relative, absolute, and ``./`` spellings of one lock file must map to a single singleton instance, deadlock-registry
    entry, and removal key. Resolving the whole path with ``realpath`` would follow a final symlink and alias a lock
    target the backend deliberately rejects, so the registry identity would differ from the backend's. Resolving only
    the parent directory and re-appending the literal final component collapses the equivalent spellings while keeping a
    final symlink a distinct key. On Windows the parent is resolved with ``abspath`` so junctions and reparse points are
    not followed either.
    """
    parent, name = os.path.split(os.fspath(path))
    return os.path.join(_resolve_dir(parent or os.curdir), name)  # ruff:ignore[os-path-join]  # string join matches abspath/realpath


class _ThreadLocalRegistry(local):
    def __init__(self) -> None:
        super().__init__()
        self.held: dict[Hashable, int] = {}


_registry: Final[_ThreadLocalRegistry] = _ThreadLocalRegistry()


_T = TypeVar("_T", bound="BaseFileLock")


class FileLockMeta(ABCMeta):
    _instances: WeakValueDictionary[str, BaseFileLock]
    _instances_lock: RLock
    _instances_under_construction: set[str]

    def __call__(  # ruff:ignore[too-many-arguments]  # forwards the public constructor's documented parameters
        cls: type[_T],
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        mode: int = _UNSET_FILE_MODE,
        thread_local: bool = True,  # ruff:ignore[boolean-type-hint-positional-argument, boolean-default-value-positional-argument]  # public API: positional bool kept for backwards compatibility
        *,
        blocking: bool = True,
        is_singleton: bool = False,
        poll_interval: float = 0.05,
        lifetime: float | None = None,
        context_error_policy: ContextErrorPolicy = "chain",
        close_error_policy: CloseErrorPolicy = "default",
        fallback_to_soft: bool = True,
        preserve_lock_file: bool = False,
        on_acquired: Callable[[int], None] | None = None,
        **kwargs: _ExtraValue,
    ) -> _T:
        _ensure_current_process()
        lifetime = _resolve_lifetime(lifetime, cls, stacklevel=cls._constructor_lifetime_warning_stacklevel)
        # Validate before building the instance: a raise inside __init__ would leave a half-constructed object whose
        # __del__ then trips over the missing context.
        context_error_policy = _resolve_context_error_policy(context_error_policy)
        close_error_policy = _resolve_close_error_policy(close_error_policy)
        preserve_lock_file = _resolve_preserve_lock_file(
            preserve=preserve_lock_file, supported=cls._preserve_lock_file_supported, cls_name=cls.__name__
        )
        on_acquired = _resolve_on_acquired(on_acquired, supported=cls._on_acquired_supported, cls_name=cls.__name__)
        params: dict[str, _LockInitValue | _ExtraValue] = {
            "timeout": timeout,
            "mode": mode,
            "thread_local": thread_local,
            "blocking": blocking,
            "is_singleton": is_singleton,
            "poll_interval": poll_interval,
            "lifetime": lifetime,
            "context_error_policy": context_error_policy,
            "close_error_policy": close_error_policy,
            "fallback_to_soft": fallback_to_soft,
            "preserve_lock_file": preserve_lock_file,
            "on_acquired": on_acquired,
            **kwargs,
        }
        if not is_singleton:
            return cls._create_instance(lock_file, params)

        # Look up, build and store under one lock. Without it two threads racing the first construction for a
        # path both miss the cache and each build their own instance, so callers relying on is_singleton for
        # reentrant locking across instances end up with two "singletons" and acquire()'s deadlock check then
        # rejects a legitimate reentrant acquire; the unguarded writes to the WeakValueDictionary are a data
        # race besides. ReadWriteLock and SoftReadWriteLock already guard their singleton caches this way.
        # Key the cache on the canonical form so equivalent spellings of one path share a singleton, and it matches the
        # deadlock-registry key acquire() uses.
        singleton_key = _canonical(lock_file)
        with cls._instances_lock:
            if (instance := cls._instances.get(singleton_key)) is None:
                if singleton_key in cls._instances_under_construction:  # pragma: needs fork
                    msg = f"Singleton lock construction is already active for {lock_file!s}"
                    raise RuntimeError(msg)
                construction_registry = cls._instances_under_construction
                construction_pid = os.getpid()
                construction_registry.add(singleton_key)
                try:
                    instance = cls._create_instance(lock_file, params)
                finally:
                    construction_registry.discard(singleton_key)
                if os.getpid() != construction_pid:  # pragma: needs fork
                    msg = "Lock construction cannot continue after fork; construct a new lock in the child"
                    raise RuntimeError(msg)
                cls._instances[singleton_key] = instance
                return instance

        params_to_check = {
            "thread_local": (thread_local, instance.is_thread_local()),
            "timeout": (timeout, instance.timeout),
            "mode": (mode, instance._context.mode),  # ruff:ignore[private-member-access]  # compares against the managed instance's own context
            "blocking": (blocking, instance.blocking),
            "poll_interval": (poll_interval, instance.poll_interval),
            "lifetime": (lifetime, instance.lifetime),
            "context_error_policy": (context_error_policy, instance.context_error_policy),
            "close_error_policy": (close_error_policy, instance.close_error_policy),
            "fallback_to_soft": (fallback_to_soft, instance.fallback_to_soft),
            "preserve_lock_file": (preserve_lock_file, instance.preserve_lock_file),
        }
        non_matching_params = {
            name: (passed_param, set_param)
            for name, (passed_param, set_param) in params_to_check.items()
            if passed_param != set_param
        }
        # Callables compare by identity, not equality: two equal callables can close over different state, so a
        # singleton must reject a different hook object even if it compares equal. Keep it out of the scalar dict above.
        hook_mismatch = on_acquired is not instance.on_acquired
        if not non_matching_params and not hook_mismatch:
            return instance  # ty: ignore[invalid-return-type]  # https://github.com/astral-sh/ty/issues/3231

        msg = "Singleton lock instances cannot be initialized with differing arguments"
        msg += "\nNon-matching arguments: "
        for param_name, (passed_param, set_param) in non_matching_params.items():
            msg += f"\n\t{param_name} (existing lock has {set_param} but {passed_param} was passed)"
        if hook_mismatch:
            msg += f"\n\ton_acquired (existing lock has {instance.on_acquired} but {on_acquired} was passed)"
        raise ValueError(msg)

    def _create_instance(
        cls: type[_T], lock_file: str | os.PathLike[str], params: dict[str, _LockInitValue | _ExtraValue]
    ) -> _T:
        model = _init_parameter_model(cls)
        if model.accepts_kwargs:
            return super().__call__(lock_file, **params)

        unsupported = sorted(
            name
            for name, value in params.items()
            if name not in model.accepted_params
            and ((parameter := model.default_params.get(name)) is None or value != parameter.default)
        )
        if unsupported:
            msg = f"{cls.__name__} does not support non-default lock options: {', '.join(unsupported)}"
            raise TypeError(msg)
        # virtualenv narrows a BaseFileLock descendant's signature; omit base defaults it does not accept (#340).
        return super().__call__(
            lock_file,
            **{name: value for name, value in params.items() if name in model.accepted_params},
        )


_INIT_PARAMETER_MODELS: Final[WeakKeyDictionary[type[BaseFileLock], _InitParameterModel]] = WeakKeyDictionary()


def _init_parameter_model(cls: type[BaseFileLock]) -> _InitParameterModel:
    # A strong cache would keep dynamically created subclasses alive for the process lifetime.
    with _fork_transition(), _FORK_STATE.parameter_models_lock:
        if (model := _INIT_PARAMETER_MODELS.get(cls)) is None:
            parameters = inspect.signature(cls.__init__).parameters.values()
            model = _InitParameterModel(
                accepted_params=frozenset(
                    parameter.name
                    for parameter in parameters
                    if parameter.kind in {inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY}
                ),
                accepts_kwargs=any(parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in parameters),
                default_params={
                    name: parameter
                    for name, parameter in inspect.signature(type(cls).__call__).parameters.items()
                    if parameter.default is not inspect.Parameter.empty
                },
            )
            _INIT_PARAMETER_MODELS[cls] = model
        return model


@dataclass(frozen=True)
class _InitParameterModel:
    accepted_params: frozenset[str]
    accepts_kwargs: bool
    default_params: dict[str, inspect.Parameter]


def _resolve_lifetime(lifetime: float | None, cls: type[BaseFileLock], *, stacklevel: int) -> float | None:
    """
    Validate ``lifetime`` and drop a value the backend cannot honor.

    ``lifetime`` is a deliberate age-based lease: a lock file older than ``lifetime`` is broken even while its holder is
    still alive. Existence locks (:class:`SoftFileLock`) implement that behavior by unlinking a reclaimable pathname,
    which can overlap a live holder. A native OS lock lives on the inode, so unlinking the pathname by age cannot revoke
    the kernel lock; a contender would lock a fresh inode and overlap the live holder (#590). Ignore the request with a
    warning rather than accept a setting that breaks mutual exclusion.
    """
    if lifetime is not None:
        if isinstance(lifetime, bool) or not isinstance(lifetime, (int, float)):
            msg = f"lifetime must be a finite non-negative number or None, not {type(lifetime).__name__}"
            raise TypeError(msg)
        if lifetime < 0 or (isinstance(lifetime, float) and not math.isfinite(lifetime)):
            msg = f"lifetime must be finite and non-negative, not {lifetime!r}"
            raise ValueError(msg)
    if lifetime is not None and not cls._lifetime_supported:
        warnings.warn(
            f"lifetime is ignored for {cls.__name__}: {cls._lifetime_unsupported_reason}; "
            f"only SoftFileLock supports lifetime-based expiry",
            stacklevel=stacklevel,
        )
        return None
    if lifetime is not None and cls._lifetime_replacements is not None:
        strict_lock, lease = cls._lifetime_replacements
        warnings.warn(
            f"{cls.__name__}(lifetime=...) uses age-based expiry and can overlap a live holder; "
            f"use {lease} for expiry or {strict_lock} for fail-closed locking",
            SoftFileLockLifetimeWarning,
            stacklevel=stacklevel,
        )
    return lifetime


def _resolve_context_error_policy(policy: str) -> ContextErrorPolicy:
    if policy not in _CONTEXT_ERROR_POLICIES:
        msg = f"context_error_policy must be 'chain' or 'group', got {policy!r}"
        raise ValueError(msg)
    if policy == "group":  # fail fast at construction rather than only when a dual failure happens to occur
        try:
            _exception_group_cls()
        except ImportError as exc:  # pragma: no cover  # only on 3.10 without the exceptiongroup backport
            msg = "context_error_policy='group' requires Python 3.11+ or the 'exceptiongroup' backport installed"
            raise ValueError(msg) from exc
    return cast("ContextErrorPolicy", policy)


def _resolve_close_error_policy(policy: str) -> CloseErrorPolicy:
    if policy not in _CLOSE_ERROR_POLICIES:
        msg = f"close_error_policy must be 'default', 'raise', or 'suppress', got {policy!r}"
        raise ValueError(msg)
    return cast("CloseErrorPolicy", policy)


def _resolve_preserve_lock_file(*, preserve: bool, supported: bool, cls_name: str) -> bool:
    # An existence lock unlinks its marker to release, so preserving the pathname would defeat unlocking. Reject the
    # request rather than silently ignore it, since a caller asking for a stable identity must know it cannot be kept.
    if preserve and not supported:
        msg = f"preserve_lock_file=True is not supported by {cls_name}: unlinking its marker is how it releases"
        raise ValueError(msg)
    return preserve


def _resolve_on_acquired(
    on_acquired: Callable[[int], None] | None, *, supported: bool, cls_name: str
) -> Callable[[int], None] | None:
    if on_acquired is None:
        return None
    # An existence lock stores protocol state in its marker, so a caller writing through the descriptor would corrupt
    # stale detection and ownership metadata; only native locks lend out the descriptor.
    if not supported:
        msg = f"on_acquired is not supported by {cls_name}: only native locks expose the lock descriptor"
        raise ValueError(msg)
    # A hook that fails and then also fails to release surfaces both errors as a BaseExceptionGroup. Require that class
    # at construction rather than at the rare moment both fail, matching how context_error_policy='group' validates.
    try:
        _exception_group_cls()
    except ImportError as exc:  # pragma: no cover  # only on 3.10 without the exceptiongroup backport
        msg = "on_acquired requires Python 3.11+ or the 'exceptiongroup' backport for its rollback error path"
        raise ValueError(msg) from exc
    return on_acquired


class BaseFileLock(contextlib.ContextDecorator, metaclass=FileLockMeta):  # ruff:ignore[too-many-public-methods]  # public config properties
    """
    Abstract base class for a file lock object.

    Provides the common reentrant API and state management. Subclasses implement the locking mechanism
    (:class:`UnixFileLock <filelock.UnixFileLock>`, :class:`WindowsFileLock <filelock.WindowsFileLock>`,
    :class:`SoftFileLock <filelock.SoftFileLock>`).

    """

    _instances: WeakValueDictionary[str, BaseFileLock]
    _instances_lock: RLock
    _instances_under_construction: set[str]

    #: How the cross-instance deadlock message names the conflicting holder; the async subclass says "task".
    _deadlock_holder_desc: str = "FileLock instance in this thread"

    #: Whether an age-based :attr:`lifetime` lease may break this lock. Only existence locks set it (they reclaim by
    #: unlinking a pathname); native OS locks leave it ``False`` since a kernel lock cannot be revoked by file age.
    _lifetime_supported: bool = False

    #: Strict-lock and lease replacements for a backend with legacy age-based expiry.
    _lifetime_replacements: tuple[str, str] | None = None

    #: Why a backend that refuses ``lifetime`` cannot honor it, named in the warning that drops the value.
    _lifetime_unsupported_reason: str = "a native OS lock cannot be broken safely by file age"

    #: Async construction adds one metaclass frame before lifetime validation.
    _constructor_lifetime_warning_stacklevel: int = 3

    #: Whether :attr:`preserve_lock_file` may be ``True``. Native locks keep the pathname on release, so they support
    #: it; existence locks unlink their marker to release and reject it.
    _preserve_lock_file_supported: bool = True

    #: Whether an :attr:`on_acquired` hook may be set. Native locks lend the descriptor out; existence locks keep
    #: protocol state in the marker and reject it.
    _on_acquired_supported: bool = True

    #: Whether a shared instance serializes its physical acquire and release behind one gate. A backend that publishes
    #: several files per owner needs it; a single-file backend is atomic and leaves it off to skip the gate entirely.
    _serialize_transitions: bool = False

    #: Ceiling in seconds on the jittered backoff between contended acquisition retries. ``0`` keeps the fixed
    #: poll cadence; a multi-file backend sets it so contending processes desynchronize instead of livelocking.
    _poll_backoff_cap: float = 0.0

    def __init_subclass__(cls, **kwargs: _SubclassValue) -> None:
        """Give each lock subclass its own singleton registry and lock."""
        super().__init_subclass__(**kwargs)
        cls._instances = WeakValueDictionary()
        cls._instances_lock = RLock()
        cls._instances_under_construction = set()
        _register_fork_class(cls)

    @classmethod
    def _reset_class_after_fork(cls) -> None:  # pragma: forked child
        cls._instances = WeakValueDictionary()
        cls._instances_lock = RLock()
        cls._instances_under_construction = set()

    def __init__(  # ruff:ignore[too-many-arguments]  # public constructor: one parameter per documented lock option
        self,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        mode: int = _UNSET_FILE_MODE,
        thread_local: bool = True,  # ruff:ignore[boolean-type-hint-positional-argument, boolean-default-value-positional-argument]  # public API: positional bool kept for backwards compatibility
        *,
        blocking: bool = True,
        is_singleton: bool = False,
        poll_interval: float = 0.05,
        lifetime: float | None = None,
        context_error_policy: ContextErrorPolicy = "chain",
        close_error_policy: CloseErrorPolicy = "default",
        fallback_to_soft: bool = True,
        preserve_lock_file: bool = False,
        on_acquired: Callable[[int], None] | None = None,
    ) -> None:
        """
        Create a new lock object.

        :param lock_file: path to the file
        :param timeout: default timeout when acquiring the lock, in seconds. It will be used as fallback value in the
            acquire method, if no timeout value (``None``) is given. If you want to disable the timeout, set it to a
            negative value. A timeout of 0 means that there is exactly one attempt to acquire the file lock.
        :param mode: file permissions for the lockfile. When not specified, the OS controls permissions via umask and
            default ACLs, preserving POSIX default ACL inheritance in shared directories.
        :param thread_local: Whether this object's internal context should be thread local or not. If this is set to
            ``False`` then the lock will be reentrant across threads. When ``True`` (the default), **all fields of the
            lock's internal context are per-thread**, including the configuration values ``poll_interval``, ``timeout``,
            ``blocking``, ``mode``, and ``lifetime``. Setting one of these properties from one thread does not change
            the value seen by another thread; threads that did not perform the write continue to see the value supplied
            at construction time. If you need configuration values to be visible across threads, construct the lock
            with ``thread_local=False``.
        :param blocking: whether the lock should be blocking or not
        :param is_singleton: If this is set to ``True`` then only one instance of this class will be created per lock
            file. This is useful if you want to use the lock object for reentrant locking without needing to pass the
            same object around.
        :param poll_interval: default interval for polling the lock file, in seconds. It will be used as fallback value
            in the acquire method, if no poll_interval value (``None``) is given.
        :param lifetime: for :class:`SoftFileLock`, the age in seconds after which a waiting process may delete the
            marker, even while its holder remains alive. This legacy expiry mode does not provide strict mutual
            exclusion. ``None`` (the default) disables age-based expiry. Native OS locks (:class:`FileLock`) cannot be
            revoked by file age and ignore a non-``None`` ``lifetime`` with a warning.
        :param context_error_policy: how a context manager reconciles a failure in its body with a failure while
            releasing on exit. ``"chain"`` (the default) keeps Python's behavior: the release error propagates with the
            body error in its ``__context__``. ``"group"`` raises a :class:`BaseExceptionGroup` holding the body error
            first and the release error second, so neither hides the other.
        :param close_error_policy: what to do with an ``os.close`` failure after relinquishing descriptor ownership.
            ``"default"`` keeps each backend's historical behavior (Unix native locks drop a FUSE/Docker ``EIO``;
            Windows native locks and :class:`SoftFileLock` propagate); ``"raise"`` always propagates the ``OSError``;
            ``"suppress"`` always ignores it. Held state is released either way. It does not affect unlock failures or
            lock-file deletion.
        :param fallback_to_soft: for :class:`UnixFileLock`, whether to switch to :class:`SoftFileLock` when the
            filesystem's ``flock`` returns ``ENOSYS``. ``True`` (the default) keeps the historical fallback;
            ``False`` fails closed, letting the ``ENOSYS`` propagate so a caller that needs kernel-enforced
            locking is never silently downgraded. It has no effect on Windows or :class:`SoftFileLock`.
        :param preserve_lock_file: for native locks (:class:`FileLock`), whether filelock promises not to unlink the
            lock pathname on release. ``False`` (the default) keeps each backend's cleanup: Windows removes the lock
            file, Unix already leaves it. ``True`` keeps a stable file identity for ACLs, auditing, and holder metadata:
            Windows skips its post-release unlink and Unix refuses to enter the ``ENOSYS`` soft fallback (which releases
            by unlinking). :class:`SoftFileLock` rejects ``True``. The promise covers filelock's own release path only;
            it cannot stop another process or the filesystem from removing the pathname.
        :param on_acquired: for native locks (:class:`FileLock`), a callable invoked with the borrowed lock descriptor
            once per physical acquisition, after filelock holds the native lock and finished backend initialization but
            before :meth:`~BaseFileLock.acquire` returns. Recursive acquisitions do not call it again. The callback may
            read, write, seek, truncate, or set metadata through ``os`` on the descriptor, but must not close, unlock,
            or take ownership of it, and filelock does not fsync its writes. If it raises, filelock releases the lock
            and re-raises. :class:`SoftFileLock` rejects the hook.

        """
        self._creator_pid = os.getpid()
        self._transition_lock = RLock()
        self._is_thread_local = thread_local
        self._is_singleton = is_singleton
        self._context_error_policy = context_error_policy  # already validated by the metaclass
        self._close_error_policy = close_error_policy  # already validated by the metaclass
        self._fallback_to_soft = fallback_to_soft
        self._preserve_lock_file = preserve_lock_file  # already validated by the metaclass
        self._on_acquired = on_acquired  # already validated by the metaclass

        self._context: FileLockContext = (ThreadLocalFileContext if thread_local else FileLockContext)(
            lock_file=os.fspath(lock_file),
            timeout=timeout,
            mode=mode,
            blocking=blocking,
            poll_interval=poll_interval,
            lifetime=lifetime,
        )
        _register_fork_object(self)

    def is_thread_local(self) -> bool:
        """:returns: a flag indicating if this lock is thread local or not"""
        return self._is_thread_local

    @property
    def is_singleton(self) -> bool:
        """
        A flag indicating if this lock is singleton or not.

        .. versionadded:: 3.13.0

        """
        return self._is_singleton

    @property
    def context_error_policy(self) -> ContextErrorPolicy:
        """
        How a context manager reconciles a body failure with a release failure on exit.

        .. versionadded:: 3.30.0

        """
        return self._context_error_policy

    @property
    def close_error_policy(self) -> CloseErrorPolicy:
        """
        What a lock does with an ``os.close`` failure after relinquishing descriptor ownership.

        .. versionadded:: 3.30.0

        """
        return self._close_error_policy

    def _close_released_fd(self, fd: int, *, default_suppresses: bool) -> None:
        # CPython never retries close() after EINTR because the descriptor number may already be reused, so neither does
        # this. close_error_policy decides the error's fate after the backend relinquishes descriptor ownership.
        try:
            os.close(fd)
        except OSError:
            if self._close_error_policy == "suppress" or (self._close_error_policy == "default" and default_suppresses):
                return
            raise

    @property
    def fallback_to_soft(self) -> bool:
        """
        Whether a :class:`FileLock` falls back to :class:`SoftFileLock` when the filesystem lacks ``flock``.

        Only :class:`UnixFileLock` acts on it: when ``False`` an ``ENOSYS`` from ``flock`` propagates instead of
        switching to existence-lock semantics.

        .. versionadded:: 3.30.0

        """
        return self._fallback_to_soft

    @property
    def preserve_lock_file(self) -> bool:
        """
        Whether filelock promises not to unlink the lock pathname on release.

        When ``True``, Windows skips its post-release unlink and Unix refuses the ``ENOSYS`` soft fallback.
        :class:`SoftFileLock` rejects ``True`` because unlinking its marker is how it releases.

        .. versionadded:: 3.30.0

        """
        return self._preserve_lock_file

    @property
    def on_acquired(self) -> Callable[[int], None] | None:
        """
        The callback run with the borrowed lock descriptor once per physical acquisition, or ``None``.

        Native locks only. It runs after the native lock is held and backend initialization finished, before
        :meth:`~BaseFileLock.acquire` returns; a raise rolls back the acquisition. :class:`SoftFileLock` rejects it.

        .. versionadded:: 3.30.0

        """
        return self._on_acquired

    @property
    def lock_file(self) -> str:
        """Path to the lock file."""
        return self._context.lock_file

    @property
    def timeout(self) -> float:
        """
        The default timeout value, in seconds.

        .. versionadded:: 2.0.0

        """
        return self._context.timeout

    @timeout.setter
    def timeout(self, value: float | str) -> None:
        """
        Change the default timeout value.

        :param value: the new value, in seconds

        """
        self._context.timeout = float(value)

    @property
    def blocking(self) -> bool:
        """
        Whether the locking is blocking or not.

        .. versionadded:: 3.14.0

        """
        return self._context.blocking

    @blocking.setter
    def blocking(self, value: bool) -> None:
        """
        Change the default blocking value.

        :param value: the new value as bool

        """
        self._context.blocking = value

    @property
    def poll_interval(self) -> float:
        """
        The default polling interval, in seconds.

        .. versionadded:: 3.24.0

        """
        return self._context.poll_interval

    @poll_interval.setter
    def poll_interval(self, value: float) -> None:
        """
        Change the default polling interval.

        :param value: the new value, in seconds

        """
        self._context.poll_interval = value

    @property
    def lifetime(self) -> float | None:
        """
        The soft marker age in seconds that permits expiry, or ``None`` to disable age-based expiry.

        A non-``None`` value permits a waiter to enter while the previous holder remains active, so it does not provide
        strict mutual exclusion. Native locks ignore the value with a warning.

        .. versionadded:: 3.24.0

        """
        return self._context.lifetime

    @lifetime.setter
    def lifetime(self, value: float | None) -> None:
        """
        Change the legacy age-based expiry threshold.

        :param value: the new value in seconds, or ``None`` to disable expiration

        :raises ValueError: if *value* is negative or not finite
        :raises TypeError: if *value* is not ``None`` and not a real number

        """
        self._context.lifetime = _resolve_lifetime(value, type(self), stacklevel=3)

    @property
    def mode(self) -> int:
        """The file permissions for the lockfile."""
        return 0o644 if self._context.mode == _UNSET_FILE_MODE else self._context.mode

    @property
    def has_explicit_mode(self) -> bool:
        """Whether the file permissions were explicitly set."""
        return self._context.mode != _UNSET_FILE_MODE

    def _open_mode(self) -> int:
        """Mode for ``os.open``: 0o666 when unset so umask and ACLs decide, otherwise the explicit mode."""
        return 0o666 if self._context.mode == _UNSET_FILE_MODE else self._context.mode

    @property
    def is_locked(self) -> bool:
        """
        A boolean indicating if the lock file is holding the lock currently.

        .. versionchanged:: 2.0.0

            This was previously a method and is now a property.

        """
        _ensure_current_process()
        return self._context.lock_file_fd is not None

    @property
    def lock_counter(self) -> int:
        """The number of times this lock has been acquired (but not yet released)."""
        _ensure_current_process()
        return self._context.lock_counter

    def __enter__(self) -> Self:
        """
        Acquire the lock.

        :returns: the lock object

        """
        self.acquire()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Release the lock, reconciling a release failure with any body failure per :attr:`context_error_policy`."""
        self._release_in_context(exc_value)

    def _release_in_context(self, body_error: BaseException | None) -> None:
        # Release from a context-manager exit. "chain" lets a release failure propagate with the body error already in
        # its __context__ (Python's default); "group" raises both as sibling leaves so neither one hides the other.
        try:
            self.release()
        except BaseException as release_error:
            if body_error is None or self._context_error_policy == "chain":
                raise
            _raise_body_and_release(body_error, release_error)

    def __del__(self) -> None:
        """Force-release so a dropped reference never leaks a held lock."""
        if vars(self).get("_creator_pid") != os.getpid():
            return  # pragma: forked child
        # A finalizer must not raise. A release error during garbage collection would otherwise surface as an
        # unraisable-exception warning, attributed to whichever code triggered collection. The dropped lock still gets
        # best-effort cleanup; an explicit release() reports the same error to a caller who can act on it.
        with contextlib.suppress(Exception):
            self.release(force=True)

    def acquire(
        self,
        timeout: float | None = None,
        poll_interval: float | None = None,
        *,
        poll_intervall: float | None = None,
        blocking: bool | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> AcquireReturnProxy:
        """
        Try to acquire the file lock.

        :param timeout: maximum wait time for acquiring the lock, ``None`` means use the default :attr:`~timeout` is and
            if ``timeout < 0``, there is no timeout and this method will block until the lock could be acquired
        :param poll_interval: interval of trying to acquire the lock file, ``None`` means use the default
            :attr:`~poll_interval`
        :param poll_intervall: deprecated, kept for backwards compatibility, use ``poll_interval`` instead
        :param blocking: defaults to True. If False, function will return immediately if it cannot obtain a lock on the
            first attempt. Otherwise, this method will block until the timeout expires or the lock is acquired.
        :param cancel_check: a callable returning ``True`` when the acquisition should be canceled. Checked on each poll
            iteration. When triggered, raises :class:`~Timeout` just like an expired timeout.

        :returns: a context object that will unlock the file when the context is exited

        :raises Timeout: if fails to acquire lock within the timeout period

        .. code-block:: python

            # You can use this method in the context manager (recommended)
            with lock.acquire():
                pass

            # Or use an equivalent try-finally construct:
            lock.acquire()
            try:
                pass
            finally:
                lock.release()

        .. versionchanged:: 2.0.0

            This method returns now a *proxy* object instead of *self*, so that it can be used in a with statement
            without side effects.

        """
        self._raise_if_inherited()
        if timeout is None:
            timeout = self._context.timeout

        if blocking is None:
            blocking = self._context.blocking

        if poll_intervall is not None:
            msg = "use poll_interval instead of poll_intervall"
            warnings.warn(msg, DeprecationWarning, stacklevel=2)
            poll_interval = poll_intervall

        poll_interval = poll_interval if poll_interval is not None else self._context.poll_interval

        start_time = time.perf_counter()
        # Wait for admission before touching any state: a caller refused entry must leave the counter, the registry and
        # the descriptor exactly as it found them.
        with self._transition_admission(
            blocking=blocking,
            cancel_check=cancel_check,
            timeout=timeout,
            poll_interval=poll_interval,
            start_time=start_time,
        ):
            # Bump the counter up front; _undo_acquire rolls it back if acquisition fails.
            self._context.lock_counter += 1

            canonical = _canonical(self.lock_file)
            self._raise_if_would_deadlock(canonical, timeout=timeout, blocking=blocking)

            try:
                self._poll_until_acquired(
                    blocking=blocking,
                    cancel_check=cancel_check,
                    timeout=timeout,
                    poll_interval=poll_interval,
                    start_time=start_time,
                )
            except BaseException:
                self._reconcile_failed_acquire(canonical)
                raise
            self._commit_acquire(canonical)
            return AcquireReturnProxy(lock=self)

    @contextlib.contextmanager
    def _transition_admission(
        self,
        *,
        blocking: bool,
        cancel_check: Callable[[], bool] | None,
        timeout: float,
        poll_interval: float,
        start_time: float,
    ) -> Generator[None]:
        # One thread at a time drives the physical transition of a shared instance. A protocol that publishes several
        # files per owner leaves a half-built claim visible otherwise, and a second thread would read it as a holder.
        # Only such a backend opts in; a single-file backend and a thread-local context each need no gate.
        if not self._serialize_transitions or self._is_thread_local:
            yield
            return
        while not self._transition_lock.acquire(blocking=False):  # pragma: needs hard-link
            if not blocking or (cancel_check is not None and cancel_check()):
                raise Timeout(self.lock_file)
            if timeout >= 0 and time.perf_counter() - start_time >= timeout:
                raise Timeout(self.lock_file)
            time.sleep(poll_interval)
        try:  # pragma: needs hard-link
            yield
        finally:  # pragma: needs hard-link
            self._transition_lock.release()

    def release(self, force: bool = False) -> None:  # ruff:ignore[boolean-type-hint-positional-argument, boolean-default-value-positional-argument]  # public API: positional bool kept for backwards compatibility
        """
        Release the file lock. The lock is only completely released when the lock counter reaches 0. The lock file
        itself may be deleted automatically, the behavior is platform-specific.

        :param force: If true, the lock counter is ignored and the lock is released in every case.

        """
        # A shared instance releases under the same gate its acquisition ran through, so a thread entering the lock
        # never observes a partially torn-down owner.
        serialize = self._serialize_transitions and not self._is_thread_local
        with self._transition_lock if serialize else contextlib.nullcontext():
            if self._creator_pid != os.getpid() or not self.is_locked:
                return
            if not force and self._context.lock_counter > 1:
                self._context.lock_counter -= 1
                return

            lock_id, lock_filename = id(self), self.lock_file
            _LOGGER.debug("Attempting to release lock %s on %s", lock_id, lock_filename)
            try:
                self._release_with_fork_tracking()
            except BaseException:
                # A failure after the OS unlock (during close or unlink) still released the lock: the backend cleared
                # its descriptor, so commit the counter and registry to released even as the cleanup error propagates.
                # A failure that left the lock held keeps the counter so a later release can retry the OS unlock.
                if not self.is_locked:
                    self._commit_release()
                raise
            self._commit_release()
            _LOGGER.debug("Lock %s released on %s", lock_id, lock_filename)

    def _raise_if_inherited(self) -> None:
        if self._creator_pid != os.getpid():  # pragma: forked child
            msg = f"{type(self).__name__} on {self.lock_file} was inherited across fork; construct a new instance"
            raise RuntimeError(msg)

    def _mark_descriptor_owned(self, fd: int, identity: tuple[int, int] | None = None) -> None:
        self._context.pending_lock_file_fd = None
        self._context.pending_lock_file_fd_identity = None
        self._context.lock_file_fd = fd
        self._context.lock_file_fd_identity = identity

    def _mark_descriptor_pending(self, fd: int, identity: tuple[int, int] | None = None) -> None:
        self._context.pending_lock_file_fd = fd
        self._context.pending_lock_file_fd_identity = identity

    def _mark_descriptor_released(self) -> None:
        self._context.pending_lock_file_fd = None
        self._context.pending_lock_file_fd_identity = None
        self._context.lock_file_fd = None
        self._context.lock_file_fd_identity = None

    def _reset_after_fork_in_child(self) -> None:  # pragma: forked child
        # fork copies the lock in whatever state the parent's threads left it, so give the child an unheld one.
        self._transition_lock = RLock()
        self._context.owner_claim_paths = ()
        self._context.claim_root = None
        self._context.lock_file_fd = None
        self._context.lock_file_fd_token = None
        self._context.lock_file_fd_identity = None
        self._context.pending_lock_file_fd = None
        self._context.pending_lock_file_fd_identity = None
        self._context.lock_counter = 0
        self._context.lock_file_key = None

    def _descriptors_for_fork(self) -> tuple[tuple[int, tuple[int, int] | None], ...]:  # pragma: needs fork
        descriptors: list[tuple[int, tuple[int, int] | None]] = []
        if self._context.lock_file_fd is not None and self._context.lock_file_fd_token is None:
            descriptors.append((self._context.lock_file_fd, self._context.lock_file_fd_identity))
        if self._context.pending_lock_file_fd is not None:
            descriptors.append((self._context.pending_lock_file_fd, self._context.pending_lock_file_fd_identity))
        return tuple(descriptors)

    def _raise_if_would_deadlock(self, canonical: str, *, timeout: float, blocking: bool) -> None:
        """
        Fail fast when a *different* live instance already holds this path in the current deadlock scope.

        Only the first, indefinitely-blocking acquire can self-deadlock this way: waiting in the OS primitive would
        block on a lock this flow already owns. A finite timeout or ``blocking=False`` keeps the normal Timeout path.
        """
        would_block = self._context.lock_counter == 1 and not self.is_locked and timeout < 0 and blocking
        if would_block and _registry.held.get(self._registry_key(canonical)) not in {None, id(self)}:
            self._context.lock_counter -= 1
            msg = (
                f"Deadlock: lock '{self.lock_file}' is already held by a different {self._deadlock_holder_desc}. "
                f"Use is_singleton=True to enable reentrant locking across instances."
            )
            raise RuntimeError(msg)

    def _registry_key(self, canonical: str) -> Hashable:
        return canonical if (scope := self._deadlock_scope()) is None else (scope, canonical)

    @staticmethod
    def _deadlock_scope() -> Hashable | None:
        """
        Execution unit whose own hold would deadlock a new acquire.

        ``None`` scopes holders to the thread, which the thread-local registry already separates. Async locks
        override it with the running task, since one event loop thread runs many tasks and only a reacquire from
        the *same* task can self-deadlock.
        """
        return None

    def _poll_until_acquired(
        self,
        *,
        blocking: bool,
        cancel_check: Callable[[], bool] | None,
        timeout: float,
        poll_interval: float,
        start_time: float,
    ) -> None:
        lock_id = id(self)
        lock_filename = self.lock_file
        attempt = 0
        while True:
            self._raise_if_inherited()
            if not self.is_locked:
                self._try_break_expired_lock()
                _LOGGER.debug("Attempting to acquire lock %s on %s", lock_id, lock_filename)
                self._acquire_with_fork_tracking()
                self._raise_if_inherited()
            if self.is_locked:
                _LOGGER.debug("Lock %s acquired on %s", lock_id, lock_filename)
                return
            if self._check_give_up(
                blocking=blocking,
                cancel_check=cancel_check,
                timeout=timeout,
                start_time=start_time,
            ):
                raise Timeout(lock_filename)
            attempt += 1
            delay = self._poll_delay(poll_interval, attempt)
            msg = "Lock %s not acquired on %s, waiting %s seconds ..."
            _LOGGER.debug(msg, lock_id, lock_filename, delay)
            time.sleep(delay)

    def _poll_delay(self, poll_interval: float, attempt: int) -> float:
        # A single-file lock retries on a fixed cadence. A backend that publishes several files per acquisition sets a
        # cap, and then contending processes back off across a jittered, exponentially widening window instead of
        # colliding on every poll; poll_interval stays the floor so a lone waiter is still responsive.
        if not self._poll_backoff_cap:
            return poll_interval
        # Cap the exponent before doubling: under heavy contention attempt reaches the thousands, and 2**attempt would
        # overflow the float multiply long before the window itself stops growing past the cap.
        window = min(
            self._poll_backoff_cap, poll_interval * 2 ** min(attempt, _MAX_BACKOFF_EXPONENT)
        )  # pragma: needs hard-link
        return max(poll_interval, secrets.randbelow(int(window * 1_000_000) + 1) / 1_000_000)  # pragma: needs hard-link

    def _reconcile_failed_acquire(self, canonical: str) -> None:
        # An acquire that raised while still holding the native lock (a hook that failed and whose rollback could not
        # release) must keep the registry entry so a later release can retry the OS unlock; otherwise roll the counter
        # back. is_locked was already reconciled by whichever release ran.
        if self.is_locked:
            self._commit_acquire(canonical)
        else:
            self._undo_acquire()

    def _invoke_on_acquired(self) -> None:
        # The wrapper runs in the backend executor for async locks, preserving the callback's documented thread.
        if self._on_acquired is None or self._context.lock_counter != 1:
            return
        try:
            self._on_acquired(cast("int", self._context.lock_file_fd))
        except BaseException as callback_error:  # arbitrary caller code; roll back on any failure
            callback_context = callback_error.__context__
            try:
                self._release_with_fork_tracking()
            except BaseException as release_error:  # ruff:ignore[blind-except]  # both errors surface via the group below
                _raise_body_and_release(callback_error, release_error)
            callback_error.__context__ = callback_context
            raise

    def _acquire_with_fork_tracking(self) -> None:
        with _fork_transition(self):
            try:
                self._acquire()
            except BaseException as acquisition_error:
                self._rollback_failed_acquire(acquisition_error)
                raise
            try:
                self._register_context_descriptor()
            except BaseException as registration_error:  # pragma: needs fork
                self._rollback_failed_registration(registration_error)
                raise
        if self.is_locked:
            self._invoke_on_acquired()

    def _rollback_failed_acquire(self, acquisition_error: BaseException) -> None:
        if not self.is_locked:
            return
        registration_error: BaseException | None = None
        tracking_error: BaseException | None = None
        try:
            self._register_context_descriptor()
        except BaseException as error:  # ruff:ignore[blind-except]  # preserve registration and acquisition failures
            registration_error = error
            try:
                # Rollback may fail too; retain the fd so a child can close it without another identity probe.
                self._register_unverified_context_descriptor()
            except BaseException as error:  # ruff:ignore[blind-except]  # pragma: no cover - allocation/control-flow during fallback
                tracking_error = error
        try:
            self._release_with_fork_tracking()
        except BaseException as rollback_error:  # ruff:ignore[blind-except]  # preserve rollback and acquisition failures
            _raise_cleanup_errors(
                "lock acquisition cleanup failed",
                acquisition_error,
                registration_error,
                tracking_error,
                rollback_error,
            )
        if registration_error is not None:  # pragma: needs fork
            _raise_cleanup_errors(
                "lock acquisition cleanup failed", acquisition_error, registration_error, tracking_error
            )

    def _rollback_failed_registration(self, registration_error: BaseException) -> None:  # pragma: needs fork
        tracking_error: BaseException | None = None
        try:
            # Rollback may fail too; retain the fd so a child can close it without another identity probe.
            self._register_unverified_context_descriptor()
        except BaseException as error:  # ruff:ignore[blind-except]  # pragma: no cover - allocation/control-flow during fallback
            tracking_error = error
        try:
            self._release_with_fork_tracking()
        except BaseException as rollback_error:  # ruff:ignore[blind-except]  # preserve rollback and registration failures
            _raise_cleanup_errors(
                "descriptor registration cleanup failed", registration_error, tracking_error, rollback_error
            )
        if tracking_error is not None:  # pragma: no cover - requires failed in-memory fallback
            _raise_cleanup_errors("descriptor registration cleanup failed", registration_error, tracking_error)

    def _release_with_fork_tracking(self) -> None:
        with _fork_transition(self):
            try:
                self._release()
            finally:
                self._unregister_released_descriptor()

    def _register_context_descriptor(self) -> None:
        if self._context.lock_file_fd is not None and self._context.lock_file_fd_token is None:
            self._context.lock_file_fd_token = _register_owned_descriptor(
                self._context.lock_file_fd,
                self._context.lock_file_fd_identity,
            )

    def _register_unverified_context_descriptor(self) -> None:
        # The rollback only reaches here still holding a descriptor it never managed to register.
        if self._context.lock_file_fd is not None and self._context.lock_file_fd_token is None:  # pragma: no branch
            self._context.lock_file_fd_token = _register_unverified_owned_descriptor(self._context.lock_file_fd)

    def _unregister_released_descriptor(self) -> None:
        if self._context.lock_file_fd is None:
            if (token := self._context.lock_file_fd_token) is not None:  # pragma: needs fork
                _unregister_owned_descriptor(token)
            self._context.lock_file_fd_token = None
            self._context.lock_file_fd_identity = None

    def _undo_acquire(self) -> None:
        """Roll back the counter after a failed acquire, dropping the registry entry once nothing holds the path."""
        self._context.lock_counter = max(0, self._context.lock_counter - 1)
        if self._context.lock_counter == 0:
            self._drop_registry_entry()

    def _commit_acquire(self, canonical: str) -> None:
        """Record this instance as the holder once the first acquire succeeds, so peers can detect the deadlock."""
        if self._context.lock_counter == 1:
            key = self._registry_key(canonical)
            # The holder scope is resolved once at commit so a later release from another flow drops the right entry.
            self._context.lock_file_key = key
            _registry.held[key] = id(self)

    def _drop_registry_entry(self) -> None:
        """Forget the key owned by this hold without resolving a mutable path again."""
        key = self._context.lock_file_key
        self._context.lock_file_key = None
        if key is not None:
            _registry.held.pop(key, None)

    def _commit_release(self) -> None:
        """Record the lock as fully released: reset the recursion counter and drop the deadlock-registry entry."""
        self._context.lock_counter = 0
        self._drop_registry_entry()

    def _try_break_expired_lock(self) -> None:
        """Remove the lock file if its modification time exceeds the configured :attr:`lifetime`."""
        if (lifetime := self._context.lifetime) is None:
            return
        with contextlib.suppress(OSError):
            # lstat, not stat: an attacker with write access to the lock directory can replace a held
            # lock file with a symlink pointing at an old file, making stat() report the target's stale
            # mtime so a waiter breaks a live lock and two processes hold it at once. lstat reads the
            # symlink's own mtime, matching the O_NOFOLLOW reads elsewhere.
            st = os.lstat(self.lock_file)
            if time.time() - st.st_mtime < lifetime:
                return
            break_lock_file(self.lock_file, st.st_mtime, st.st_ino)

    def _check_give_up(
        self,
        *,
        blocking: bool,
        cancel_check: Callable[[], bool] | None,
        timeout: float,
        start_time: float,
    ) -> bool:
        lock_id, lock_filename = id(self), self.lock_file
        if blocking is False:
            _LOGGER.debug("Failed to immediately acquire lock %s on %s", lock_id, lock_filename)
            return True
        if cancel_check is not None and cancel_check():
            _LOGGER.debug("Cancellation requested for lock %s on %s", lock_id, lock_filename)
            return True
        if 0 <= timeout < time.perf_counter() - start_time:
            _LOGGER.debug("Timeout on acquiring lock %s on %s", lock_id, lock_filename)
            return True
        return False

    @abstractmethod
    def _acquire(self) -> None:
        """If the file lock could be acquired, self._context.lock_file_fd holds the file descriptor of the lock file."""
        raise NotImplementedError

    @abstractmethod
    def _release(self) -> None:
        """Releases the lock and sets self._context.lock_file_fd to None."""
        raise NotImplementedError


# acquire() returns this wrapper instead of self so entering the with-statement does not call __enter__ a second
# time; returning self would re-acquire the lock in BaseFileLock.__enter__ without a matching release (issue #37).
class AcquireReturnProxy:
    """A context-aware object that will release the lock file when exiting."""

    def __init__(self, lock: BaseFileLock | ReadWriteLock | SoftReadWriteLock) -> None:
        self.lock: BaseFileLock | ReadWriteLock | SoftReadWriteLock = lock

    def __enter__(self) -> BaseFileLock | ReadWriteLock | SoftReadWriteLock:
        return self.lock

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        if isinstance(self.lock, BaseFileLock):
            self.lock._release_in_context(exc_value)  # ruff:ignore[private-member-access]  # forwards __exit__ to the owned lock's context release
        else:  # a reader/writer lock does not carry a context_error_policy
            self.lock.release()


@dataclass
class FileLockContext:
    """Holds the context for a ``BaseFileLock`` object."""

    # A separate class so ThreadLocalFileContext can make the whole context thread-local.

    lock_file: str
    timeout: float
    mode: int
    blocking: bool
    poll_interval: float

    #: The lock lifetime in seconds; ``None`` means the lock never expires.
    lifetime: float | None = None

    #: File descriptor from os.open for the lock file; not None while the lock is held.
    lock_file_fd: int | None = None

    #: Registry token for the descriptor owned by this thread's context.
    lock_file_fd_token: int | None = None

    #: Identity captured by a backend that already inspected the descriptor.
    lock_file_fd_identity: tuple[int, int] | None = None

    #: Descriptor opened by a backend but not yet committed as the held lock.
    pending_lock_file_fd: int | None = None

    #: Identity captured for a descriptor whose acquisition has not committed.
    pending_lock_file_fd_identity: tuple[int, int] | None = None

    #: Depth of nested acquisitions; the lock is released only when it returns to 0.
    lock_counter: int = 0

    #: Canonical registry key captured when the first physical acquisition commits.
    lock_file_key: Hashable | None = None

    #: Claim pathnames this owner published, removed by name on release so no holder ever unlinks a peer's claim.
    owner_claim_paths: tuple[str, ...] = ()

    #: Canonical lock path resolved when an acquisition starts. A waiter polling a relative path must keep publishing
    #: into the directory it started in, even when another thread changes the working directory mid-wait.
    claim_root: str | None = None


class ThreadLocalFileContext(FileLockContext, local):
    """A thread local version of the ``FileLockContext`` class."""


@dataclass(frozen=True)
class _OwnedDescriptor:
    fd: int
    creator_pid: int
    device: int | None
    inode: int | None


class _ForkTransitionContext(local):
    depth: int = 0


class _ForkState:
    def __init__(self) -> None:
        self.gate = Condition(RLock())
        self.registry_lock = RLock()
        self.parameter_models_lock = RLock()
        self.transition_context = _ForkTransitionContext()
        self.transitions: dict[int, dict[int, _ForkDescriptorOwner | None]] = {}
        self.active_transitions = 0
        self.admission_closed = False
        self.fork_owner_depths: dict[int, int] = {}
        self.pinned_objects: dict[int, list[tuple[_ForkResettable, ...]]] = {}
        self.pinned_classes: dict[int, list[tuple[_ForkResettableClass, ...]]] = {}
        self.provisional_descriptor_tokens: dict[int, list[tuple[int, ...]]] = {}
        self.pid = os.getpid()

    def reset_synchronization(self) -> None:  # pragma: forked child
        self.gate = Condition(RLock())
        self.registry_lock = RLock()
        self.parameter_models_lock = RLock()
        self.transition_context = _ForkTransitionContext()
        self.transitions = {}
        self.active_transitions = 0
        self.admission_closed = False
        self.fork_owner_depths = {}
        self.pinned_objects = {}
        self.pinned_classes = {}
        self.provisional_descriptor_tokens = {}


_FORK_OBJECTS: Final[WeakValueDictionary[int, _ForkResettable]] = WeakValueDictionary()
_FORK_CLASSES: Final[WeakValueDictionary[int, _ForkResettableClass]] = WeakValueDictionary()
_OWNED_DESCRIPTORS: Final[dict[int, _OwnedDescriptor]] = {}
_DESCRIPTOR_TOKENS: Final[count[int]] = count()
_TRANSITION_TOKENS: Final[count[int]] = count()
_FORK_STATE: Final = _ForkState()
_FORK_AUDIT_EVENTS: Final[frozenset[str]] = frozenset({"os.fork", "os.forkpty"})


def _register_fork_hooks() -> None:
    if _REGISTER_AT_FORK is None:
        return  # pragma: lacks fork
    sys.addaudithook(_audit_fork_safety)  # pragma: needs fork
    _REGISTER_AT_FORK(  # pragma: needs fork
        before=_pin_fork_objects,
        after_in_parent=_resume_parent_after_fork,
        after_in_child=_reset_child_after_fork,
    )


@contextmanager
def _fork_transition(descriptor_owner: _ForkDescriptorOwner | None = None) -> Generator[None]:
    if not _HAS_REGISTER_AT_FORK:
        yield  # pragma: lacks fork
        return  # pragma: lacks fork
    _ensure_current_process()  # pragma: needs fork
    creator_pid = os.getpid()  # pragma: needs fork
    token = _enter_fork_transition(descriptor_owner)  # pragma: needs fork
    try:  # pragma: needs fork
        yield
    finally:
        if os.getpid() == creator_pid:  # pragma: needs fork
            _leave_fork_transition(token)


def _register_fork_object(instance: _ForkResettable) -> None:
    if not _HAS_REGISTER_AT_FORK:
        return  # pragma: lacks fork
    with _fork_transition(), _FORK_STATE.registry_lock:  # pragma: needs fork
        _FORK_OBJECTS[id(instance)] = instance
        _refresh_owner_pins()


def _register_fork_class(cls: _ForkResettableClass) -> None:
    if not _HAS_REGISTER_AT_FORK:
        return  # pragma: lacks fork
    with _fork_transition(), _FORK_STATE.registry_lock:  # pragma: needs fork
        _FORK_CLASSES[id(cls)] = cls
        _refresh_owner_pins()


def _register_owned_descriptor(fd: int, identity: tuple[int, int] | None = None) -> int | None:
    if not _HAS_REGISTER_AT_FORK:
        return None  # pragma: lacks fork
    with _fork_transition():  # pragma: needs fork
        if identity is None:
            stat_result = os.fstat(fd)
            identity = stat_result.st_dev, stat_result.st_ino
        return _record_owned_descriptor(fd, identity)


def _register_unverified_owned_descriptor(fd: int) -> int | None:
    if not _HAS_REGISTER_AT_FORK:
        return None  # pragma: lacks fork
    with _fork_transition():  # pragma: needs fork
        return _record_owned_descriptor(fd, None)


def _record_owned_descriptor(fd: int, identity: tuple[int, int] | None) -> int:  # pragma: needs fork
    with _FORK_STATE.registry_lock:
        token = next(_DESCRIPTOR_TOKENS)
        _OWNED_DESCRIPTORS[token] = _OwnedDescriptor(
            fd=fd,
            creator_pid=os.getpid(),
            device=None if identity is None else identity[0],
            inode=None if identity is None else identity[1],
        )
    return token


def _unregister_owned_descriptor(token: int) -> None:  # pragma: needs fork
    with _fork_transition(), _FORK_STATE.registry_lock:
        _OWNED_DESCRIPTORS.pop(token, None)


def _pin_fork_objects() -> None:  # pragma: needs fork
    _ensure_current_process()
    thread_id = get_ident()
    with _FORK_STATE.gate:
        _FORK_STATE.fork_owner_depths[thread_id] = _FORK_STATE.fork_owner_depths.get(thread_id, 0) + 1
        _FORK_STATE.admission_closed = True
        while _FORK_STATE.active_transitions > len(_FORK_STATE.transitions.get(thread_id, ())):
            _FORK_STATE.gate.wait()
        transition_owners = tuple(_FORK_STATE.transitions.get(thread_id, {}).values())
    _verify_unverified_descriptors()
    owners = {id(owner): owner for owner in transition_owners if owner is not None}
    provisional_descriptor_tokens = tuple(
        starmap(
            _snapshot_descriptor_for_fork,
            (descriptor for owner in owners.values() for descriptor in owner._descriptors_for_fork()),  # ruff:ignore[private-member-access]  # snapshots each owner's own fork descriptors
        )
    )
    with _FORK_STATE.registry_lock:
        _FORK_STATE.provisional_descriptor_tokens.setdefault(thread_id, []).append(provisional_descriptor_tokens)
        _FORK_STATE.pinned_objects.setdefault(thread_id, []).append(tuple(_FORK_OBJECTS.values()))
        _FORK_STATE.pinned_classes.setdefault(thread_id, []).append(tuple(_FORK_CLASSES.values()))


def _resume_parent_after_fork() -> None:  # pragma: needs fork
    thread_id = get_ident()
    with _FORK_STATE.registry_lock:
        for token in _FORK_STATE.provisional_descriptor_tokens[thread_id].pop():
            _OWNED_DESCRIPTORS.pop(token, None)
        _FORK_STATE.pinned_objects[thread_id].pop()
        _FORK_STATE.pinned_classes[thread_id].pop()
        if not _FORK_STATE.provisional_descriptor_tokens[thread_id]:  # pragma: needs fork
            del _FORK_STATE.provisional_descriptor_tokens[thread_id]
            del _FORK_STATE.pinned_objects[thread_id]
            del _FORK_STATE.pinned_classes[thread_id]
    with _FORK_STATE.gate:
        if _FORK_STATE.fork_owner_depths[thread_id] == 1:
            del _FORK_STATE.fork_owner_depths[thread_id]
        else:  # pragma: no cover - earlier at-fork callbacks may deadlock first
            _FORK_STATE.fork_owner_depths[thread_id] -= 1
        _FORK_STATE.admission_closed = bool(_FORK_STATE.fork_owner_depths)
        _FORK_STATE.gate.notify_all()


def _ensure_current_process() -> None:
    _reset_child_after_fork()


def _reset_child_after_fork() -> None:  # pragma: forked child
    if (pid := os.getpid()) == _FORK_STATE.pid:
        return
    thread_id = get_ident()
    pinned_objects = _FORK_STATE.pinned_objects.get(thread_id, [()])[-1]
    pinned_classes = _FORK_STATE.pinned_classes.get(thread_id, [()])[-1]
    descriptors: list[_OwnedDescriptor] = []
    for token, descriptor in tuple(_OWNED_DESCRIPTORS.items()):
        if descriptor.creator_pid != pid:
            descriptors.append(_OWNED_DESCRIPTORS.pop(token))
    _FORK_STATE.reset_synchronization()
    _FORK_STATE.pid = pid
    _INIT_PARAMETER_MODELS.clear()
    _detach_child_state(descriptors, pinned_objects, pinned_classes)


def _enter_fork_transition(descriptor_owner: _ForkDescriptorOwner | None) -> int:  # pragma: needs fork
    thread_id = get_ident()
    with _FORK_STATE.gate:
        while (
            _FORK_STATE.admission_closed
            and thread_id not in _FORK_STATE.fork_owner_depths
            and thread_id not in _FORK_STATE.transitions
        ):
            _FORK_STATE.gate.wait()  # pragma: no cover - exercised in isolated interpreter
        token = next(_TRANSITION_TOKENS)
        _FORK_STATE.transitions.setdefault(thread_id, {})[token] = descriptor_owner
        _FORK_STATE.active_transitions += 1
    _FORK_STATE.transition_context.depth += 1
    return token


def _leave_fork_transition(token: int) -> None:  # pragma: needs fork
    _FORK_STATE.transition_context.depth -= 1
    with _FORK_STATE.gate:
        thread_id = get_ident()
        del _FORK_STATE.transitions[thread_id][token]
        if not _FORK_STATE.transitions[thread_id]:
            del _FORK_STATE.transitions[thread_id]
        _FORK_STATE.active_transitions -= 1
        _FORK_STATE.gate.notify_all()


def _verify_unverified_descriptors() -> None:  # pragma: needs fork
    with _FORK_STATE.registry_lock:
        for token, descriptor in tuple(_OWNED_DESCRIPTORS.items()):
            if descriptor.creator_pid != os.getpid() or descriptor.device is not None:
                continue
            try:
                stat_result = os.fstat(descriptor.fd)
            except OSError:
                continue
            _OWNED_DESCRIPTORS[token] = _OwnedDescriptor(
                fd=descriptor.fd,
                creator_pid=descriptor.creator_pid,
                device=stat_result.st_dev,
                inode=stat_result.st_ino,
            )


def _snapshot_descriptor_for_fork(fd: int, identity: tuple[int, int] | None) -> int:  # pragma: needs fork
    if identity is None:  # pragma: needs fork
        try:
            stat_result = os.fstat(fd)
        except OSError:
            pass
        else:
            identity = stat_result.st_dev, stat_result.st_ino
    return _record_owned_descriptor(fd, identity)


def _detach_child_state(  # pragma: forked child
    descriptors: list[_OwnedDescriptor],
    pinned_objects: tuple[_ForkResettable, ...],
    pinned_classes: tuple[_ForkResettableClass, ...],
) -> None:
    for descriptor in descriptors:
        if descriptor.device is None:
            continue
        try:
            stat_result = os.fstat(descriptor.fd)
        except OSError:
            continue
        if (stat_result.st_dev, stat_result.st_ino) != (descriptor.device, descriptor.inode):
            continue
        with contextlib.suppress(OSError):
            os.close(descriptor.fd)
    for instance in pinned_objects:
        instance._reset_after_fork_in_child()  # ruff:ignore[private-member-access]  # resets each pinned instance in the fork child
    for cls in pinned_classes:
        cls._reset_class_after_fork()
    _registry.held.clear()


def _refresh_owner_pins() -> None:  # pragma: needs fork
    with _FORK_STATE.gate:
        fork_owner_thread_ids = tuple(_FORK_STATE.fork_owner_depths)
    if get_ident() not in fork_owner_thread_ids:  # pragma: no cover - at-fork callbacks disable tracing
        return
    objects = tuple(_FORK_OBJECTS.values())
    classes = tuple(_FORK_CLASSES.values())
    for thread_id in fork_owner_thread_ids:
        if object_snapshots := _FORK_STATE.pinned_objects.get(thread_id):  # pragma: needs fork
            object_snapshots[:] = [objects] * len(object_snapshots)
            _FORK_STATE.pinned_classes[thread_id][:] = [classes] * len(object_snapshots)


# Audit events carry unrelated interpreter-owned values; object is the sole accurate common element type.
def _audit_fork_safety(  # pragma: no cover - CPython disables tracing while Python audit hooks run
    event: str, _args: tuple[object, ...]
) -> None:
    if event in _FORK_AUDIT_EVENTS:
        if _FORK_STATE.transition_context.depth or _FORK_STATE.fork_owner_depths:
            msg = f"{event} is unsafe while filelock is changing descriptor ownership"
            raise RuntimeError(msg)
    elif event == "_posixsubprocess.fork_exec" and _FORK_STATE.transition_context.depth:
        msg = "fork_exec is unsafe while filelock is changing descriptor ownership"
        raise RuntimeError(msg)


_register_fork_hooks()


__all__ = [
    "_UNSET_FILE_MODE",
    "AcquireReturnProxy",
    "BaseFileLock",
    "CloseErrorPolicy",
    "ContextErrorPolicy",
    "FileLockContext",
    "FileLockMeta",
    "LockOptions",
    "_append_exception_context",
    "_canonical",
    "_ensure_current_process",
    "_fork_transition",
    "_grouped_errors",
    "_raise_body_and_release",
    "_raise_chained_errors",
    "_raise_cleanup_errors",
    "_raise_grouped_errors",
    "_register_fork_class",
    "_register_fork_object",
    "_register_owned_descriptor",
    "_unregister_owned_descriptor",
]
