from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache, wraps
from typing import TYPE_CHECKING, Generic, TypeVar

from dask import config
from dask._compatibility import importlib_metadata
from dask.utils import funcname

if TYPE_CHECKING:
    from typing_extensions import ParamSpec

    BackendFuncParams = ParamSpec("BackendFuncParams")
    BackendFuncReturn = TypeVar("BackendFuncReturn")


class DaskBackendEntrypoint:
    """Base Collection-Backend Entrypoint Class

    Most methods in this class correspond to collection-creation
    for a specific library backend. Once a collection is created,
    the existing data will be used to dispatch compute operations
    within individual tasks. The backend is responsible for
    ensuring that these data-directed dispatch functions are
    registered when ``__init__`` is called.
    """

    @classmethod
    def to_backend_dispatch(cls):
        """Return a dispatch function to move data to this backend"""
        raise NotImplementedError

    @staticmethod
    def to_backend(data):
        """Create a new collection with this backend"""
        raise NotImplementedError


@lru_cache(maxsize=1)
def detect_entrypoints(entry_point_name):
    return {
        ep.name: ep for ep in importlib_metadata.entry_points(group=entry_point_name)
    }


BackendEntrypointType = TypeVar(
    "BackendEntrypointType",
    bound="DaskBackendEntrypoint",
)


class CreationDispatch(Generic[BackendEntrypointType]):
    """Simple backend dispatch for collection-creation functions"""

    _lookup: dict[str, BackendEntrypointType]
    _module_name: str
    _config_field: str
    _default: str
    _entrypoint_class: type[BackendEntrypointType]
    _entrypoint_root: str

    def __init__(
        self,
        module_name: str,
        default: str,
        entrypoint_class: type[BackendEntrypointType],
        name: str | None = None,
        entrypoint_root: str = "dask",
    ):
        self._lookup = {}
        self._module_name = module_name
        self._config_field = f"{module_name}.backend"
        self._default = default
        self._entrypoint_class = entrypoint_class
        self._entrypoint_root = entrypoint_root
        if name:
            self.__name__ = name

    def register_backend(
        self, name: str, backend: BackendEntrypointType
    ) -> BackendEntrypointType:
        """Register a target class for a specific array-backend label"""
        if not isinstance(backend, self._entrypoint_class):
            raise ValueError(
                f"This CreationDispatch only supports "
                f"{self._entrypoint_class} registration. "
                f"Got {type(backend)}"
            )
        self._lookup[name] = backend
        return backend

    def dispatch(self, backend: str):
        """Return the desired backend entrypoint"""
        try:
            impl = self._lookup[backend]
        except KeyError:
            # Check entrypoints for the specified backend
            entrypoints = detect_entrypoints(
                f"{self._entrypoint_root}.{self._module_name}.backends"
            )
            if backend in entrypoints:
                return self.register_backend(backend, entrypoints[backend].load()())
        else:
            return impl
        raise ValueError(f"No backend dispatch registered for {backend}")

    @property
    def backend(self) -> str:
        """Return the desired collection backend"""
        return config.get(self._config_field, self._default) or self._default

    @backend.setter
    def backend(self, value: str):
        raise RuntimeError(
            f"Set the backend by configuring the {self._config_field} option"
        )

    def register_inplace(
        self,
        backend: str,
        name: str | None = None,
    ) -> Callable[
        [Callable[BackendFuncParams, BackendFuncReturn]],
        Callable[BackendFuncParams, BackendFuncReturn],
    ]:
        """Register dispatchable function"""

        def decorator(
            fn: Callable[BackendFuncParams, BackendFuncReturn]
        ) -> Callable[BackendFuncParams, BackendFuncReturn]:
            dispatch_name = name or fn.__name__
            dispatcher = self.dispatch(backend)
            dispatcher.__setattr__(dispatch_name, fn)

            @wraps(fn)
            def wrapper(*args, **kwargs):
                func = getattr(self, dispatch_name)
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    try:
                        exc = type(e)(
                            f"An error occurred while calling the {funcname(func)} "
                            f"method registered to the {self.backend} backend.\n"
                            f"Original Message: {e}"
                        )
                    except TypeError:
                        raise e
                    else:
                        raise exc from e

            wrapper.__name__ = dispatch_name
            return wrapper

        return decorator

    def __getattr__(self, item: str):
        """
        Return the appropriate attribute for the current backend
        """
        backend = self.dispatch(self.backend)
        return getattr(backend, item)
