from __future__ import annotations

import sys
from functools import cache
from typing import TYPE_CHECKING, Any, Protocol, cast

from narwhals._compliant import CompliantNamespace
from narwhals._typing_compat import TypeVar

if TYPE_CHECKING:
    from collections.abc import Iterator
    from importlib.metadata import EntryPoints

    from typing_extensions import LiteralString, TypeAlias

    from narwhals._compliant.typing import (
        CompliantDataFrameAny,
        CompliantFrameAny,
        CompliantLazyFrameAny,
        CompliantSeriesAny,
    )
    from narwhals.utils import Version


__all__ = ["Plugin", "from_native"]

CompliantAny: TypeAlias = (
    "CompliantDataFrameAny | CompliantLazyFrameAny | CompliantSeriesAny"
)
"""A statically-unknown, Compliant object originating from a plugin."""

FrameT = TypeVar(
    "FrameT",
    bound="CompliantFrameAny",
    default="CompliantDataFrameAny | CompliantLazyFrameAny",
)
FromNativeR_co = TypeVar(
    "FromNativeR_co", bound=CompliantAny, covariant=True, default=CompliantAny
)


@cache
def _discover_entrypoints() -> EntryPoints:
    from importlib.metadata import entry_points as eps

    group = "narwhals.plugins"
    if sys.version_info < (3, 10):
        return cast("EntryPoints", eps().get(group, ()))
    return eps(group=group)


class PluginNamespace(CompliantNamespace[FrameT, Any], Protocol[FrameT, FromNativeR_co]):
    def from_native(self, data: Any, /) -> FromNativeR_co: ...


class Plugin(Protocol[FrameT, FromNativeR_co]):
    @property
    def NATIVE_PACKAGE(self) -> LiteralString: ...  # noqa: N802

    def __narwhals_namespace__(
        self, version: Version
    ) -> PluginNamespace[FrameT, FromNativeR_co]: ...
    def is_native(self, native_object: object, /) -> bool: ...


@cache
def _might_be(cls: type, type_: str) -> bool:  # pragma: no cover
    try:
        return any(type_ in o.__module__.split(".") for o in cls.mro())
    except TypeError:
        return False


def _is_native_plugin(native_object: Any, plugin: Plugin) -> bool:
    pkg = plugin.NATIVE_PACKAGE
    return (
        sys.modules.get(pkg) is not None
        and _might_be(type(native_object), pkg)  # type: ignore[arg-type]
        and plugin.is_native(native_object)
    )


def _iter_from_native(native_object: Any, version: Version) -> Iterator[CompliantAny]:
    for entry_point in _discover_entrypoints():
        plugin: Plugin = entry_point.load()
        if _is_native_plugin(native_object, plugin):
            compliant_namespace = plugin.__narwhals_namespace__(version=version)
            yield compliant_namespace.from_native(native_object)


def from_native(native_object: Any, version: Version) -> CompliantAny | None:
    """Attempt to convert `native_object` to a Compliant object, using any available plugin(s).

    Arguments:
        native_object: Raw object from user.
        version: Narwhals API version.

    Returns:
        If the following conditions are met
            - at least 1 plugin is installed
            - at least 1 installed plugin supports `type(native_object)`

        Then for the **first matching plugin**, the result of the call below.
        This *should* be an object accepted by a Narwhals Dataframe, Lazyframe, or Series:

            plugin: Plugin
            plugin.__narwhals_namespace__(version).from_native(native_object)

        In all other cases, `None` is returned instead.
    """
    return next(_iter_from_native(native_object, version), None)


def _show_suggestions(native_object_type: type) -> str | None:
    if _might_be(native_object_type, "daft"):  # pragma: no cover
        return (
            "Hint: it looks like you passed a `daft.DataFrame` but don't have `narwhals-daft` installed.\n"
            "Please refer to https://github.com/narwhals-dev/narwhals-daft for installation instructions."
        )
    return None
