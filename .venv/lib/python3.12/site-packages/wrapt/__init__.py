"""
Wrapt is a library for decorators, wrappers and monkey patching.
"""

__version_info__ = ("2", "1", "1")
__version__ = ".".join(__version_info__)

from .__wrapt__ import (
    BaseObjectProxy,
    BoundFunctionWrapper,
    CallableObjectProxy,
    FunctionWrapper,
    PartialCallableObjectProxy,
    partial,
)
from .decorators import AdapterFactory, adapter_factory, decorator, synchronized
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
    "decorator",
    "synchronized",
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
