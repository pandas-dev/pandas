import sys

if sys.version_info >= (3, 10):
    from inspect import FullArgSpec, Signature
    from types import ModuleType, TracebackType
    from typing import (
        Any,
        AsyncIterator,
        Callable,
        Concatenate,
        Generator,
        Generic,
        Iterator,
        ParamSpec,
        Protocol,
        TypeVar,
        overload,
    )

    # Mirrors wrapt.__all__ at runtime (src/wrapt/__init__.py). Anything not
    # listed here is a private stub-internal helper (TypeVar, Protocol,
    # type alias) and is prefixed with a leading underscore.
    __all__ = (
        "AutoObjectProxy",
        "BaseObjectProxy",
        "BoundFunctionWrapper",
        "CallableObjectProxy",
        "FunctionWrapper",
        "LazyObjectProxy",
        "ObjectProxy",
        "PartialCallableObjectProxy",
        "partial",
        "AdapterFactory",
        "adapter_factory",
        "bind_state_to_wrapper",
        "async_to_sync",
        "decorator",
        "lru_cache",
        "mark_as_async",
        "mark_as_sync",
        "sync_to_async",
        "synchronized",
        "with_signature",
        "discover_post_import_hooks",
        "notify_module_loaded",
        "register_post_import_hook",
        "when_imported",
        "apply_patch",
        "function_wrapper",
        "lazy_import",
        "patch_function_wrapper",
        "resolve_path",
        "transient_function_wrapper",
        "wrap_function_wrapper",
        "wrap_object",
        "wrap_object_attribute",
        "WeakFunctionProxy",
    )

    _P = ParamSpec("_P")
    _R = TypeVar("_R", covariant=True)

    _T = TypeVar("_T", bound=Any)

    # Need two sets of ParamSpec/TypeVar for generics in some cases to ensure
    # that mypy, pyrefly and ty all work correctly. More specifically need a
    # separate set for cases where extracting the type or instance from the
    # first argument of a callable where binding is involved.

    _P1 = ParamSpec("_P1")
    _R1 = TypeVar("_R1", covariant=True)

    _T1 = TypeVar("_T1", bound=Any)

    _P2 = ParamSpec("_P2")
    _R2 = TypeVar("_R2", covariant=True)

    _T2 = TypeVar("_T2", bound=Any)

    class _Boolean(Protocol):
        def __bool__(self) -> bool: ...

    # ObjectProxy
    #
    # BaseObjectProxy mirrors the dunder surface of the runtime class
    # (wrappers.ObjectProxy / _wrappers.ObjectProxy). Every protocol/operator
    # dunder that the runtime defines unconditionally is declared here so
    # static type checkers accept `len(proxy)`, `with proxy: ...`,
    # `proxy + x`, etc. without complaint. Dunders that the runtime only
    # attaches dynamically per wrapped object (__iter__ lives on
    # ObjectProxy for backward compatibility; __next__, __aiter__,
    # __anext__, __await__ and __length_hint__ are added per-instance by
    # AutoObjectProxy based on the wrapped interface) are intentionally
    # NOT claimed statically on BaseObjectProxy.
    #
    # `__enter__`/`__aenter__` return `Any` because the runtime forwards to
    # the wrapped object's `__enter__`, which may return a value of any type
    # (for example, `threading.Lock.__enter__` returns `bool`, not the lock
    # itself). Returning `_T` would be a false claim that the value bound by
    # `with proxy as x` is the wrapped type, which is only sometimes true.
    # Other operator dunders return `Any` for the same reason: the result
    # depends on the wrapped type.

    class BaseObjectProxy(Generic[_T]):
        __wrapped__: _T

        __name__: str
        __qualname__: str

        def __init__(self, wrapped: _T) -> None: ...
        def __getattr__(self, name: str) -> Any: ...
        def __mro_entries__(self, bases: tuple[type, ...]) -> tuple[type, ...]: ...

        # Context managers.
        def __enter__(self) -> Any: ...
        def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: BaseException | None,
            traceback: TracebackType | None,
        ) -> bool | None: ...
        async def __aenter__(self) -> Any: ...
        async def __aexit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: BaseException | None,
            traceback: TracebackType | None,
        ) -> bool | None: ...

        # Container protocol.
        def __len__(self) -> int: ...
        def __contains__(self, key: object, /) -> bool: ...
        def __getitem__(self, key: Any, /) -> Any: ...
        def __setitem__(self, key: Any, value: Any, /) -> None: ...
        def __delitem__(self, key: Any, /) -> None: ...
        def __reversed__(self) -> Iterator[Any]: ...

        # Rich comparisons.
        def __lt__(self, other: object, /) -> bool: ...
        def __le__(self, other: object, /) -> bool: ...
        def __gt__(self, other: object, /) -> bool: ...
        def __ge__(self, other: object, /) -> bool: ...

        # Numeric conversions.
        def __bool__(self) -> bool: ...
        def __int__(self) -> int: ...
        def __float__(self) -> float: ...
        def __complex__(self) -> complex: ...
        def __bytes__(self) -> bytes: ...
        def __index__(self) -> int: ...
        def __round__(self, ndigits: int | None = ..., /) -> Any: ...

        # Unary arithmetic.
        def __neg__(self) -> Any: ...
        def __pos__(self) -> Any: ...
        def __abs__(self) -> Any: ...
        def __invert__(self) -> Any: ...

        # Binary arithmetic.
        def __add__(self, other: Any, /) -> Any: ...
        def __sub__(self, other: Any, /) -> Any: ...
        def __mul__(self, other: Any, /) -> Any: ...
        def __matmul__(self, other: Any, /) -> Any: ...
        def __truediv__(self, other: Any, /) -> Any: ...
        def __floordiv__(self, other: Any, /) -> Any: ...
        def __mod__(self, other: Any, /) -> Any: ...
        def __divmod__(self, other: Any, /) -> Any: ...
        def __pow__(self, other: Any, modulo: Any = ..., /) -> Any: ...
        def __lshift__(self, other: Any, /) -> Any: ...
        def __rshift__(self, other: Any, /) -> Any: ...
        def __and__(self, other: Any, /) -> Any: ...
        def __or__(self, other: Any, /) -> Any: ...
        def __xor__(self, other: Any, /) -> Any: ...

        # Reflected binary arithmetic.
        def __radd__(self, other: Any, /) -> Any: ...
        def __rsub__(self, other: Any, /) -> Any: ...
        def __rmul__(self, other: Any, /) -> Any: ...
        def __rmatmul__(self, other: Any, /) -> Any: ...
        def __rtruediv__(self, other: Any, /) -> Any: ...
        def __rfloordiv__(self, other: Any, /) -> Any: ...
        def __rmod__(self, other: Any, /) -> Any: ...
        def __rdivmod__(self, other: Any, /) -> Any: ...
        def __rpow__(self, other: Any, modulo: Any = ..., /) -> Any: ...
        def __rlshift__(self, other: Any, /) -> Any: ...
        def __rrshift__(self, other: Any, /) -> Any: ...
        def __rand__(self, other: Any, /) -> Any: ...
        def __ror__(self, other: Any, /) -> Any: ...
        def __rxor__(self, other: Any, /) -> Any: ...

        # In-place arithmetic.
        def __iadd__(self, other: Any, /) -> Any: ...
        def __isub__(self, other: Any, /) -> Any: ...
        def __imul__(self, other: Any, /) -> Any: ...
        def __imatmul__(self, other: Any, /) -> Any: ...
        def __itruediv__(self, other: Any, /) -> Any: ...
        def __ifloordiv__(self, other: Any, /) -> Any: ...
        def __imod__(self, other: Any, /) -> Any: ...
        def __ipow__(self, other: Any, /) -> Any: ...  # type: ignore[misc]
        def __ilshift__(self, other: Any, /) -> Any: ...
        def __irshift__(self, other: Any, /) -> Any: ...
        def __iand__(self, other: Any, /) -> Any: ...
        def __ior__(self, other: Any, /) -> Any: ...
        def __ixor__(self, other: Any, /) -> Any: ...

        # Copy / pickle.
        def __copy__(self) -> Any: ...
        def __deepcopy__(self, memo: dict[int, Any], /) -> Any: ...
        def __reduce__(self) -> Any: ...

        # wrapt-specific escape hatches (not Python-language dunders despite
        # the dunder-shaped names). Ordinary attribute access, __setattr__
        # and arithmetic operations all forward to the wrapped object; these
        # hooks let subclasses and wrapt-aware callers reach around the
        # forwarding layer when they need to.
        #
        # __self_dict__: the proxy's own instance dict (since __dict__ is
        # overridden to delegate to the wrapped object).
        #
        # __self_setattr__: set an attribute directly on the proxy, bypassing
        # the forwarding __setattr__. Used for stashing state onto a wrapper.
        #
        # __object_proxy__: the class used to re-wrap results of arithmetic
        # and bitwise operations (e.g. __add__ returns
        # ``self.__object_proxy__(self.__wrapped__ + other)``). Subclasses
        # override it to control the type of proxy produced from operations.
        __self_dict__: dict[str, Any]
        @property
        def __object_proxy__(self) -> type[BaseObjectProxy[Any]]: ...
        def __self_setattr__(self, name: str, value: Any) -> None: ...

    class ObjectProxy(BaseObjectProxy[_T]):
        def __new__(cls, *args: Any, **kwargs: Any) -> ObjectProxy[_T]: ...
        def __init__(self, wrapped: _T) -> None: ...
        def __iter__(self) -> Iterator[Any]: ...

    class AutoObjectProxy(BaseObjectProxy[_T]):
        # AutoObjectProxy attaches the dunders below to a per-instance
        # subclass at construction time, based on what the wrapped object
        # exposes (see proxies.AutoObjectProxy.__new__). We statically
        # claim all of them here so type checkers accept code that uses
        # them -- AutoObjectProxy's contract is "take on the interface
        # of the wrapped object", so the stub reflects that intent at
        # the cost of being permissive for proxies whose wrapped value
        # doesn't actually support a given dunder (a runtime
        # AttributeError, which mirrors the wrapt design).

        def __new__(cls, wrapped: _T) -> AutoObjectProxy[_T]: ...
        def __init__(self, wrapped: _T) -> None: ...

        # Hook called by BaseObjectProxy.__setattr__ whenever __wrapped__
        # is reassigned. AutoObjectProxy's implementation re-wires the
        # interface-conditional dunders below to match the new wrapped
        # object; subclasses may override it to run additional fixup logic
        # when the wrapped object changes.
        def __wrapped_setattr_fixups__(self) -> None: ...
        def __call__(self, *args: Any, **kwargs: Any) -> Any: ...
        def __iter__(self) -> Iterator[Any]: ...
        def __next__(self) -> Any: ...
        def __aiter__(self) -> AsyncIterator[Any]: ...
        async def __anext__(self) -> Any: ...
        def __length_hint__(self) -> int: ...
        def __await__(self) -> Generator[Any, Any, Any]: ...

    # LazyObjectProxy

    class LazyObjectProxy(AutoObjectProxy[_T]):
        def __new__(
            cls,
            callback: Callable[[], _T] | None = None,
            *,
            interface: Any = ...,
        ) -> LazyObjectProxy[_T]: ...
        def __init__(
            self,
            callback: Callable[[], _T] | None = None,
            *,
            interface: Any = ...,
        ) -> None: ...

    @overload
    def lazy_import(name: str) -> LazyObjectProxy[ModuleType]: ...
    @overload
    def lazy_import(
        name: str, attribute: str, *, interface: Any = ...
    ) -> LazyObjectProxy[Any]: ...

    # CallableObjectProxy

    class CallableObjectProxy(BaseObjectProxy[_T]):
        def __call__(self, *args: Any, **kwargs: Any) -> Any: ...

    # PartialCallableObjectProxy

    class PartialCallableObjectProxy:
        def __init__(
            self, func: Callable[..., Any], *args: Any, **kwargs: Any
        ) -> None: ...
        def __call__(self, *args: Any, **kwargs: Any) -> Any: ...

    def partial(
        func: Callable[..., Any], /, *args: Any, **kwargs: Any
    ) -> Callable[..., Any]: ...

    # WeakFunctionProxy

    class WeakFunctionProxy:
        def __init__(
            self,
            wrapped: Callable[..., Any],
            callback: Callable[..., Any] | None = None,
        ) -> None: ...
        def __call__(self, *args: Any, **kwargs: Any) -> Any: ...

    # FunctionWrapper

    _WrappedFunction = Callable[_P, _R]

    _GenericCallableWrapperFunction = Callable[
        [_WrappedFunction[_P, _R], Any, tuple[Any, ...], dict[str, Any]], _R
    ]

    _ClassMethodWrapperFunction = Callable[
        [type[Any], _WrappedFunction[_P, _R], Any, tuple[Any, ...], dict[str, Any]], _R
    ]

    _InstanceMethodWrapperFunction = Callable[
        [Any, _WrappedFunction[_P, _R], Any, tuple[Any, ...], dict[str, Any]], _R
    ]

    _WrapperFunction = (
        _GenericCallableWrapperFunction[_P, _R]
        | _ClassMethodWrapperFunction[_P, _R]
        | _InstanceMethodWrapperFunction[_P, _R]
    )

    class _FunctionWrapperBase(ObjectProxy[_WrappedFunction[_P, _R]]):
        _self_instance: Any
        _self_wrapper: _WrapperFunction[_P, _R]
        _self_enabled: bool | _Boolean | Callable[[], bool] | None
        _self_binding: str
        _self_parent: Any
        _self_owner: Any

    class BoundFunctionWrapper(_FunctionWrapperBase[_P1, _R1]):
        def __call__(self, *args: _P1.args, **kwargs: _P1.kwargs) -> _R1: ...

        # Note that for following overloads, testing with mypy and ty they still do
        # not handle static methods being decorated but to best knowledge this is
        # a limitation in those type checkers. Testing with pyrefly fails on any
        # type of bound method. Testing with pyright handles case correctly.
        #
        # Also, note that use of _T2, _P2 and _R2 in first two cases is also required
        # to ensure correct handling by mypy and ty, so do not change to use of _T1,
        # _P1 and _R1.

        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: BoundFunctionWrapper[Concatenate[_T2, _P2], _R2],
            instance: _T2,
            owner: type[_T2] | None = None,
            /,
        ) -> BoundFunctionWrapper[_P2, _R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: BoundFunctionWrapper[Concatenate[_T2, _P2], _R2],
            instance: _T2,
            owner: type[Any] | None = None,
            /,
        ) -> BoundFunctionWrapper[_P2, _R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self, instance: None, owner: type[_T1] | None = None, /
        ) -> BoundFunctionWrapper[_P1, _R1]: ...
        @overload
        def __get__(  # Required to ensure pyright works correctly
            self, instance: _T1, owner: type[_T1] | None = None, /
        ) -> BoundFunctionWrapper[_P1, _R1]: ...

    class FunctionWrapper(_FunctionWrapperBase[_P1, _R1]):
        def __init__(
            self,
            wrapped: _WrappedFunction[_P1, _R1],
            wrapper: _WrapperFunction[_P1, _R1],
            enabled: bool | _Boolean | Callable[[], bool] | None = None,
        ) -> None: ...
        def __call__(self, *args: _P1.args, **kwargs: _P1.kwargs) -> _R1: ...

        # Note that for following overloads, testing with mypy and ty they still do
        # not handle static methods being decorated but to best knowledge this is
        # a limitation in those type checkers. Testing with pyrefly fails on any
        # type of bound method. Testing with pyright handles case correctly.
        #
        # Also, note that use of _T2, _P2 and _R2 in first two cases is also required
        # to ensure correct handling by mypy and ty, so do not change to use of _T1,
        # _P1 and _R1.

        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: FunctionWrapper[Concatenate[_T2, _P2], _R2],
            instance: _T2,
            owner: type[Any] | None = None,
            /,
        ) -> BoundFunctionWrapper[_P2, _R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: FunctionWrapper[Concatenate[_T2, _P2], _R2],
            instance: _T2,
            owner: type[_T2] | None = None,
            /,
        ) -> BoundFunctionWrapper[_P2, _R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self, instance: None, owner: type[_T1] | None = None, /
        ) -> BoundFunctionWrapper[_P1, _R1]: ...
        @overload
        def __get__(  # Required to ensure pyright works correctly
            self, instance: _T1, owner: type[_T1] | None = None, /
        ) -> BoundFunctionWrapper[_P1, _R1]: ...

    # AdapterFactory/adapter_factory()

    class AdapterFactory:
        def __call__(
            self, wrapped: Callable[..., Any]
        ) -> str | FullArgSpec | Callable[..., Any]: ...

    class _DelegatedAdapterFactory(AdapterFactory):
        def __init__(self, factory: Callable[..., Any]) -> None: ...

    adapter_factory = _DelegatedAdapterFactory

    # decorator()

    class _Descriptor(Protocol):
        def __get__(self, instance: Any, owner: type[Any] | None = None) -> Any: ...

    class _FunctionDecorator:
        @overload
        def __call__(
            self,
            callable: (
                Callable[_P, _R]
                | Callable[Concatenate[type[_T], _P], _R]  # Required for pylance
                # | Callable[Concatenate[Any, _P], _R]  # Breaks mypy, pyrefly and ty
                | Callable[[type[_T]], _R]  # Required for pylance
            ),
        ) -> FunctionWrapper[_P, _R]: ...
        @overload
        def __call__(self, callable: _Descriptor) -> FunctionWrapper[_P, Any]: ...

    class _PartialFunctionDecorator:
        @overload
        def __call__(
            self, wrapper: _GenericCallableWrapperFunction[_P, _R], /
        ) -> _FunctionDecorator: ...
        @overload
        def __call__(
            self, wrapper: _ClassMethodWrapperFunction[_P, _R], /
        ) -> _FunctionDecorator: ...
        @overload
        def __call__(
            self, wrapper: _InstanceMethodWrapperFunction[_P, _R], /
        ) -> _FunctionDecorator: ...

    # ... Decorator applied to class type.

    @overload
    def decorator(wrapper: type[_T], /) -> _FunctionDecorator: ...

    # ... Decorator applied to function or method.

    @overload
    def decorator(
        wrapper: _GenericCallableWrapperFunction[_P, _R], /
    ) -> _FunctionDecorator: ...
    @overload
    def decorator(
        wrapper: _ClassMethodWrapperFunction[_P, _R], /
    ) -> _FunctionDecorator: ...
    @overload
    def decorator(
        wrapper: _InstanceMethodWrapperFunction[_P, _R], /
    ) -> _FunctionDecorator: ...

    # ... Positional arguments.

    @overload
    def decorator(
        *,
        enabled: bool | _Boolean | Callable[[], bool] | None = None,
        adapter: str | FullArgSpec | AdapterFactory | Callable[..., Any] | None = None,
        proxy: type[FunctionWrapper[Any, Any]] | None = None,
    ) -> _PartialFunctionDecorator: ...

    # function_wrapper()

    @overload
    def function_wrapper(wrapper: type[Any]) -> _FunctionDecorator: ...
    @overload
    def function_wrapper(
        wrapper: _GenericCallableWrapperFunction[_P, _R],
    ) -> _FunctionDecorator: ...
    @overload
    def function_wrapper(
        wrapper: _ClassMethodWrapperFunction[_P, _R],
    ) -> _FunctionDecorator: ...
    @overload
    def function_wrapper(
        wrapper: _InstanceMethodWrapperFunction[_P, _R],
    ) -> _FunctionDecorator: ...
    # @overload
    # def function_wrapper(wrapper: Any) -> _FunctionDecorator: ... # Don't use, breaks stuff.

    # wrap_function_wrapper()

    def wrap_function_wrapper(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        wrapper: _WrapperFunction[_P, _R],
    ) -> FunctionWrapper[_P, _R]: ...

    # patch_function_wrapper()

    class _WrapperDecorator:
        def __call__(
            self, wrapper: _WrapperFunction[_P, _R]
        ) -> FunctionWrapper[_P, _R]: ...

    def patch_function_wrapper(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        enabled: bool | _Boolean | Callable[[], bool] | None = None,
    ) -> _WrapperDecorator: ...

    # transient_function_wrapper()

    class _TransientDecorator:
        def __call__(self, wrapper: _WrapperFunction[_P, _R]) -> _FunctionDecorator: ...

    def transient_function_wrapper(
        target: ModuleType | type[Any] | Any | str, name: str
    ) -> _TransientDecorator: ...

    # resolve_path()

    def resolve_path(
        target: ModuleType | type[Any] | Any | str, name: str
    ) -> tuple[ModuleType | type[Any] | Any, str, Callable[..., Any]]: ...

    # apply_patch()

    def apply_patch(
        parent: ModuleType | type[Any] | Any,
        attribute: str,
        replacement: Any,
    ) -> None: ...

    # wrap_object()

    _WrapperFactory = Callable[
        [Callable[..., Any], tuple[Any, ...], dict[str, Any]], type[ObjectProxy[Any]]
    ]

    def wrap_object(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        factory: _WrapperFactory | type[ObjectProxy[Any]],
        args: tuple[Any, ...] = (),
        kwargs: dict[str, Any] | None = None,
    ) -> Any: ...

    # wrap_object_attribute()

    def wrap_object_attribute(
        module: ModuleType | type[Any] | Any | str,
        name: str,
        factory: _WrapperFactory | type[ObjectProxy[Any]],
        args: tuple[Any, ...] = (),
        kwargs: dict[str, Any] | None = None,
    ) -> Any: ...

    # register_post_import_hook()

    def register_post_import_hook(
        hook: Callable[[ModuleType], Any] | str, name: str
    ) -> None: ...

    # discover_post_import_hooks()

    def discover_post_import_hooks(group: str) -> None: ...

    # notify_module_loaded()

    def notify_module_loaded(module: ModuleType) -> None: ...

    # when_imported()

    class _ImportHookDecorator:
        def __call__(self, hook: Callable[[ModuleType], Any]) -> Callable[..., Any]: ...

    def when_imported(name: str) -> _ImportHookDecorator: ...

    # synchronized()

    class _SynchronizedObject:
        def __call__(self, wrapped: Callable[_P, _R]) -> Callable[_P, _R]: ...
        def __enter__(self) -> Any: ...
        def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: BaseException | None,
            traceback: TracebackType | None,
        ) -> bool | None: ...
        async def __aenter__(self) -> Any: ...
        async def __aexit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: BaseException | None,
            traceback: TracebackType | None,
        ) -> bool | None: ...

    @overload
    def synchronized(wrapped: Callable[_P, _R]) -> Callable[_P, _R]: ...
    @overload
    def synchronized(wrapped: Any) -> _SynchronizedObject: ...

    # mark_as_sync(), mark_as_async(), async_to_sync(), sync_to_async()

    @overload
    def mark_as_sync(wrapped: Callable[_P, _R], /) -> Callable[_P, _R]: ...
    @overload
    def mark_as_sync(
        *, generator: bool | None = None
    ) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]: ...
    @overload
    def mark_as_async(wrapped: Callable[_P, _R], /) -> Callable[_P, _R]: ...
    @overload
    def mark_as_async(
        *, generator: bool | None = None
    ) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]: ...
    def async_to_sync(wrapped: Callable[_P, _R]) -> Callable[_P, _R]: ...
    def sync_to_async(wrapped: Callable[_P, _R]) -> Callable[_P, _R]: ...

    # bind_state_to_wrapper()

    class _StateBindingWrapper:
        name: str
        wrapper_factory: _Descriptor | None
        def __init__(self, *, name: str = "state") -> None: ...
        def __call__(self, wrapper_factory: _Descriptor) -> _StateBindingWrapper: ...
        def __get__(
            self, instance: Any, owner: type[Any]
        ) -> (
            _StateBindingWrapper
            | Callable[[Callable[..., Any]], FunctionWrapper[..., Any]]
        ): ...

    bind_state_to_wrapper = _StateBindingWrapper

    # lru_cache()

    class _BoundLRUCacheFunctionWrapper(BoundFunctionWrapper[_P1, _R1]):
        def cache_info(self) -> Any | None: ...
        def cache_clear(self) -> None: ...
        def cache_parameters(self) -> dict[str, Any] | None: ...

    class _LRUCacheFunctionWrapper(FunctionWrapper[_P1, _R1]):
        __bound_function_wrapper__: type[_BoundLRUCacheFunctionWrapper[_P1, _R1]]
        def cache_info(self) -> Any | None: ...
        def cache_clear(self) -> None: ...
        def cache_parameters(self) -> dict[str, Any] | None: ...

    @overload
    def lru_cache(func: Callable[_P, _R], /) -> _LRUCacheFunctionWrapper[_P, _R]: ...
    @overload
    def lru_cache(
        func: None = None, /, **kwargs: Any
    ) -> Callable[[Callable[_P, _R]], _LRUCacheFunctionWrapper[_P, _R]]: ...

    # with_signature()

    def with_signature(
        *,
        prototype: Callable[..., Any] | None = None,
        signature: Signature | None = None,
        factory: (
            Callable[[Callable[..., Any]], Signature | Callable[..., Any]] | None
        ) = None,
    ) -> Callable[[Callable[_P, _R]], FunctionWrapper[_P, _R]]: ...
