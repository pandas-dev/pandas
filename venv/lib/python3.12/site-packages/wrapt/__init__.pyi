import sys

if sys.version_info >= (3, 10):
    from inspect import FullArgSpec
    from types import ModuleType, TracebackType
    from typing import (
        Any,
        Callable,
        Concatenate,
        Generic,
        ParamSpec,
        Protocol,
        TypeVar,
        overload,
    )

    P = ParamSpec("P")
    R = TypeVar("R", covariant=True)

    T = TypeVar("T", bound=Any)

    # Need two sets of ParamSpec/TypeVar for generics in some cases to ensure
    # that mypy, pyrefly and ty all work correctly. More specifically need a
    # separate set for cases where extracting the type or instance from the
    # first argument of a callable where binding is involved.

    P1 = ParamSpec("P1")
    R1 = TypeVar("R1", covariant=True)

    T1 = TypeVar("T1", bound=Any)

    P2 = ParamSpec("P2")
    R2 = TypeVar("R2", covariant=True)

    T2 = TypeVar("T2", bound=Any)

    class Boolean(Protocol):
        def __bool__(self) -> bool: ...

    # ObjectProxy

    class BaseObjectProxy(Generic[T]):
        __wrapped__: T
        def __init__(self, wrapped: T) -> None: ...

    class ObjectProxy(BaseObjectProxy[T]):
        def __init__(self, wrapped: T) -> None: ...

    class AutoObjectProxy(BaseObjectProxy[T]):
        def __init__(self, wrapped: T) -> None: ...

    # LazyObjectProxy

    class LazyObjectProxy(AutoObjectProxy[T]):
        def __init__(
            self, callback: Callable[[], T] | None, *, interface: Any = ...
        ) -> None: ...

    @overload
    def lazy_import(name: str) -> LazyObjectProxy[ModuleType]: ...
    @overload
    def lazy_import(
        name: str, attribute: str, *, interface: Any = ...
    ) -> LazyObjectProxy[Any]: ...

    # CallableObjectProxy

    class CallableObjectProxy(BaseObjectProxy[T]):
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

    WrappedFunction = Callable[P, R]

    GenericCallableWrapperFunction = Callable[
        [WrappedFunction[P, R], Any, tuple[Any, ...], dict[str, Any]], R
    ]

    ClassMethodWrapperFunction = Callable[
        [type[Any], WrappedFunction[P, R], Any, tuple[Any, ...], dict[str, Any]], R
    ]

    InstanceMethodWrapperFunction = Callable[
        [Any, WrappedFunction[P, R], Any, tuple[Any, ...], dict[str, Any]], R
    ]

    WrapperFunction = (
        GenericCallableWrapperFunction[P, R]
        | ClassMethodWrapperFunction[P, R]
        | InstanceMethodWrapperFunction[P, R]
    )

    class _FunctionWrapperBase(ObjectProxy[WrappedFunction[P, R]]):
        _self_instance: Any
        _self_wrapper: WrapperFunction[P, R]
        _self_enabled: bool | Boolean | Callable[[], bool] | None
        _self_binding: str
        _self_parent: Any
        _self_owner: Any

    class BoundFunctionWrapper(_FunctionWrapperBase[P1, R1]):
        def __call__(self, *args: P1.args, **kwargs: P1.kwargs) -> R1: ...

        # Note that for following overloads, testing with mypy and ty they still do
        # not handle static methods being decorated but to best knowledge this is
        # a limitation in those type checkers. Testing with pyrefly fails on any
        # type of bound method. Testing with pyright handles case correctly.
        #
        # Also, note that use of T2, P2 and R2 in first two cases is also required
        # to ensure correct handling by mypy and ty, so do not change to use of T1,
        # P1 and R1.

        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: BoundFunctionWrapper[Concatenate[T2, P2], R2],
            instance: T2,
            owner: type[T2] | None = None,
        ) -> BoundFunctionWrapper[P2, R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: BoundFunctionWrapper[Concatenate[T2, P2], R2],
            instance: T2,
            owner: type[Any] | None = None,
        ) -> BoundFunctionWrapper[P2, R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self, instance: None, owner: type[T1] | None = None
        ) -> BoundFunctionWrapper[P1, R1]: ...
        @overload
        def __get__(  # Required to ensure pyright works correctly
            self, instance: T1, owner: type[T1] | None = None
        ) -> BoundFunctionWrapper[P1, R1]: ...

    class FunctionWrapper(_FunctionWrapperBase[P1, R1]):
        def __init__(
            self,
            wrapped: WrappedFunction[P1, R1],
            wrapper: WrapperFunction[P1, R1],
            enabled: bool | Boolean | Callable[[], bool] | None = None,
        ) -> None: ...
        def __call__(self, *args: P1.args, **kwargs: P1.kwargs) -> R1: ...

        # Note that for following overloads, testing with mypy and ty they still do
        # not handle static methods being decorated but to best knowledge this is
        # a limitation in those type checkers. Testing with pyrefly fails on any
        # type of bound method. Testing with pyright handles case correctly.
        #
        # Also, note that use of T2, P2 and R2 in first two cases is also required
        # to ensure correct handling by mypy and ty, so do not change to use of T1,
        # P1 and R1.

        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: FunctionWrapper[Concatenate[T2, P2], R2],
            instance: T2,
            owner: type[Any] | None = None,
        ) -> BoundFunctionWrapper[P2, R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self: FunctionWrapper[Concatenate[T2, P2], R2],
            instance: T2,
            owner: type[T2] | None = None,
        ) -> BoundFunctionWrapper[P2, R2]: ...
        @overload
        def __get__(  # Required to ensure mypy, pyrefly and ty works correctly
            self, instance: None, owner: type[T1] | None = None
        ) -> BoundFunctionWrapper[P1, R1]: ...
        @overload
        def __get__(  # Required to ensure pyright works correctly
            self, instance: T1, owner: type[T1] | None = None
        ) -> BoundFunctionWrapper[P1, R1]: ...

    # AdapterFactory/adapter_factory()

    class AdapterFactory(Protocol):
        def __call__(
            self, wrapped: Callable[..., Any]
        ) -> str | FullArgSpec | Callable[..., Any]: ...

    def adapter_factory(wrapped: Callable[..., Any]) -> AdapterFactory: ...

    # decorator()

    class Descriptor(Protocol):
        def __get__(self, instance: Any, owner: type[Any] | None = None) -> Any: ...

    class FunctionDecorator:
        @overload
        def __call__(
            self,
            callable: (
                Callable[P, R]
                | Callable[Concatenate[type[T], P], R]  # Required for pylance
                # | Callable[Concatenate[Any, P], R]  # Breaks mypy, pyrefly and ty
                | Callable[[type[T]], R]  # Required for pylance
            ),
        ) -> FunctionWrapper[P, R]: ...
        @overload
        def __call__(self, callable: Descriptor) -> FunctionWrapper[P, Any]: ...

    class PartialFunctionDecorator:
        @overload
        def __call__(
            self, wrapper: GenericCallableWrapperFunction[P, R], /
        ) -> FunctionDecorator: ...
        @overload
        def __call__(
            self, wrapper: ClassMethodWrapperFunction[P, R], /
        ) -> FunctionDecorator: ...
        @overload
        def __call__(
            self, wrapper: InstanceMethodWrapperFunction[P, R], /
        ) -> FunctionDecorator: ...

    # ... Decorator applied to class type.

    @overload
    def decorator(wrapper: type[T], /) -> FunctionDecorator: ...

    # ... Decorator applied to function or method.

    @overload
    def decorator(
        wrapper: GenericCallableWrapperFunction[P, R], /
    ) -> FunctionDecorator: ...
    @overload
    def decorator(
        wrapper: ClassMethodWrapperFunction[P, R], /
    ) -> FunctionDecorator: ...
    @overload
    def decorator(
        wrapper: InstanceMethodWrapperFunction[P, R], /
    ) -> FunctionDecorator: ...

    # ... Positional arguments.

    @overload
    def decorator(
        *,
        enabled: bool | Boolean | Callable[[], bool] | None = None,
        adapter: str | FullArgSpec | AdapterFactory | Callable[..., Any] | None = None,
        proxy: type[FunctionWrapper[Any, Any]] | None = None,
    ) -> PartialFunctionDecorator: ...

    # function_wrapper()

    @overload
    def function_wrapper(wrapper: type[Any]) -> FunctionDecorator: ...
    @overload
    def function_wrapper(
        wrapper: GenericCallableWrapperFunction[P, R],
    ) -> FunctionDecorator: ...
    @overload
    def function_wrapper(
        wrapper: ClassMethodWrapperFunction[P, R],
    ) -> FunctionDecorator: ...
    @overload
    def function_wrapper(
        wrapper: InstanceMethodWrapperFunction[P, R],
    ) -> FunctionDecorator: ...
    # @overload
    # def function_wrapper(wrapper: Any) -> FunctionDecorator: ... # Don't use, breaks stuff.

    # wrap_function_wrapper()

    def wrap_function_wrapper(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        wrapper: WrapperFunction[P, R],
    ) -> FunctionWrapper[P, R]: ...

    # patch_function_wrapper()

    class WrapperDecorator:
        def __call__(self, wrapper: WrapperFunction[P, R]) -> FunctionWrapper[P, R]: ...

    def patch_function_wrapper(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        enabled: bool | Boolean | Callable[[], bool] | None = None,
    ) -> WrapperDecorator: ...

    # transient_function_wrapper()

    class TransientDecorator:
        def __call__(self, wrapper: WrapperFunction[P, R]) -> FunctionDecorator: ...

    def transient_function_wrapper(
        target: ModuleType | type[Any] | Any | str, name: str
    ) -> TransientDecorator: ...

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

    WrapperFactory = Callable[
        [Callable[..., Any], tuple[Any, ...], dict[str, Any]], type[ObjectProxy[Any]]
    ]

    def wrap_object(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        factory: WrapperFactory | type[ObjectProxy[Any]],
        args: tuple[Any, ...],
        kwargs: dict[str, Any],
    ) -> Any: ...

    # wrap_object_attribute()

    def wrap_object_attribute(
        target: ModuleType | type[Any] | Any | str,
        name: str,
        factory: WrapperFactory | type[ObjectProxy[Any]],
        args: tuple[Any, ...] = (),
        kwargs: dict[str, Any] = {},
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

    class ImportHookDecorator:
        def __call__(self, hook: Callable[[ModuleType], Any]) -> Callable[..., Any]: ...

    def when_imported(name: str) -> ImportHookDecorator: ...

    # synchronized()

    class SynchronizedObject:
        def __call__(self, wrapped: Callable[P, R]) -> Callable[P, R]: ...
        def __enter__(self) -> Any: ...
        def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_value: BaseException | None,
            traceback: TracebackType | None,
        ) -> bool | None: ...

    @overload
    def synchronized(wrapped: Callable[P, R]) -> Callable[P, R]: ...
    @overload
    def synchronized(wrapped: Any) -> SynchronizedObject: ...
