"""
Wrapt is a library for decorators, wrappers and monkey patching.
"""


def _format_version(parts):
    base = ".".join(parts[:3])
    if len(parts) == 3:
        return base
    suffix = parts[3]
    return (
        f"{base}.{suffix}" if suffix.startswith(("dev", "post")) else f"{base}{suffix}"
    )


__version_info__ = ("2", "2", "2")
__version__ = _format_version(__version_info__)

from .__wrapt__ import (
    BaseObjectProxy,
    BoundFunctionWrapper,
    CallableObjectProxy,
    FunctionWrapper,
    PartialCallableObjectProxy,
    partial,
)
from .caching import lru_cache
from .decorators import (
    AdapterFactory,
    adapter_factory,
    bind_state_to_wrapper,
    decorator,
)
from .importer import (
    discover_post_import_hooks,
    notify_module_loaded,
    register_post_import_hook,
    when_imported,
)
from .patches import (
    apply_patch,
    function_wrapper,
    patch_function_wrapper,
    resolve_path,
    transient_function_wrapper,
    wrap_function_wrapper,
    wrap_object,
    wrap_object_attribute,
)
from .proxies import AutoObjectProxy, LazyObjectProxy, ObjectProxy, lazy_import
from .signature import with_signature
from .synchronization import (
    async_to_sync,
    mark_as_async,
    mark_as_sync,
    sync_to_async,
    synchronized,
)
from .weakrefs import WeakFunctionProxy

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
