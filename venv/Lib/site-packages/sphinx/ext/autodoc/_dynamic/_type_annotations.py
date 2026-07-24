from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.util import inspect
from sphinx.util.typing import stringify_annotation

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from sphinx.util.typing import _StringifyMode


def _record_typehints(
    *,
    autodoc_annotations: dict[str, dict[str, str]],
    name: str,
    obj: Any,
    short_literals: bool,
    type_aliases: Mapping[str, str] | None,
    unqualified_typehints: bool,
) -> None:
    """Record type hints to env object."""
    mode: _StringifyMode
    if unqualified_typehints:
        mode = 'smart'
    else:
        mode = 'fully-qualified'

    try:
        if callable(obj):
            annotation = autodoc_annotations.setdefault(name, {})
            sig = inspect.signature(obj, type_aliases=type_aliases)
            for param in sig.parameters.values():
                if param.annotation is not param.empty:
                    annotation[param.name] = stringify_annotation(
                        param.annotation, mode, short_literals=short_literals
                    )
            if sig.return_annotation is not sig.empty:
                annotation['return'] = stringify_annotation(
                    sig.return_annotation, mode, short_literals=short_literals
                )
    except (TypeError, ValueError):
        pass
