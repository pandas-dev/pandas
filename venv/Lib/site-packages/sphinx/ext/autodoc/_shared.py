"""Shared utilities for autodoc that don't have a better home."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.util import logging
from sphinx.util.inspect import safe_getattr

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence
    from typing import Any, Final, Literal, NoReturn, Protocol

    from sphinx.config import Config
    from sphinx.util.typing import _RestifyMode

    class _AttrGetter(Protocol):  # NoQA: PYI046
        def __call__(self, obj: Any, name: str, default: Any = ..., /) -> Any: ...


LOGGER: Final[logging.SphinxLoggerAdapter] = logging.getLogger('sphinx.ext.autodoc')


class _AutodocConfig:
    __slots__ = (
        'autoclass_content',
        'autodoc_class_signature',
        'autodoc_default_options',
        'autodoc_docstring_signature',
        'autodoc_inherit_docstrings',
        'autodoc_member_order',
        'autodoc_mock_imports',
        'autodoc_preserve_defaults',
        'autodoc_type_aliases',
        'autodoc_typehints',
        'autodoc_typehints_description_target',
        'autodoc_typehints_format',
        'autodoc_use_type_comments',
        # non-autodoc config
        'python_display_short_literal_types',
        'strip_signature_backslash',
    )

    autoclass_content: Literal['both', 'class', 'init']
    autodoc_class_signature: Literal['mixed', 'separated']
    autodoc_default_options: Mapping[str, str | bool]
    autodoc_docstring_signature: bool
    autodoc_inherit_docstrings: bool
    autodoc_member_order: Literal['alphabetical', 'bysource', 'groupwise']
    autodoc_mock_imports: Sequence[str]
    autodoc_preserve_defaults: bool
    autodoc_type_aliases: Mapping[str, str]
    autodoc_typehints: Literal['signature', 'description', 'none', 'both']
    autodoc_typehints_description_target: Literal[
        'all', 'documented', 'documented_params'
    ]
    autodoc_typehints_format: Literal['fully-qualified', 'short']
    autodoc_use_type_comments: bool
    # non-autodoc config
    python_display_short_literal_types: bool
    strip_signature_backslash: bool

    @classmethod
    def from_config(cls, config: Config) -> _AutodocConfig:
        return cls(
            autoclass_content=config.autoclass_content,
            autodoc_class_signature=config.autodoc_class_signature,
            autodoc_default_options=config.autodoc_default_options,
            autodoc_docstring_signature=config.autodoc_docstring_signature,
            autodoc_inherit_docstrings=config.autodoc_inherit_docstrings,
            autodoc_member_order=config.autodoc_member_order,
            autodoc_mock_imports=config.autodoc_mock_imports,
            autodoc_preserve_defaults=config.autodoc_preserve_defaults,
            autodoc_type_aliases=config.autodoc_type_aliases,
            autodoc_typehints=config.autodoc_typehints,
            autodoc_typehints_description_target=config.autodoc_typehints_description_target,
            autodoc_typehints_format=config.autodoc_typehints_format,
            autodoc_use_type_comments=config.autodoc_use_type_comments,
            python_display_short_literal_types=config.python_display_short_literal_types,
            strip_signature_backslash=config.strip_signature_backslash,
        )

    def __init__(
        self,
        *,
        autoclass_content: Literal['both', 'class', 'init'] = 'class',
        autodoc_class_signature: Literal['mixed', 'separated'] = 'mixed',
        autodoc_default_options: Mapping[str, str | bool] = {}.keys().mapping,
        autodoc_docstring_signature: bool = True,
        autodoc_inherit_docstrings: bool = True,
        autodoc_member_order: Literal['alphabetical', 'bysource', 'groupwise'] = (
            'alphabetical'
        ),
        autodoc_mock_imports: Sequence[str] = (),
        autodoc_preserve_defaults: bool = False,
        autodoc_type_aliases: Mapping[str, str] = {}.keys().mapping,
        autodoc_typehints: Literal[
            'signature', 'description', 'none', 'both'
        ] = 'signature',
        autodoc_typehints_description_target: Literal[
            'all', 'documented', 'documented_params'
        ] = 'all',
        autodoc_typehints_format: Literal['fully-qualified', 'short'] = 'short',
        autodoc_use_type_comments: bool = True,
        python_display_short_literal_types: bool = False,
        strip_signature_backslash: bool = False,
    ) -> None:
        for name in self.__slots__:
            super().__setattr__(name, locals()[name])

    def __repr__(self) -> str:
        items = ((name, getattr(self, name)) for name in self.__slots__)
        args = ', '.join(f'{name}={value!r}' for name, value in items)
        return f'_AutodocConfig({args})'

    def __setattr__(self, key: str, value: Any) -> NoReturn:
        msg = f'{self.__class__.__name__} is immutable'
        raise AttributeError(msg)

    def __delattr__(self, key: str) -> NoReturn:
        msg = f'{self.__class__.__name__} is immutable'
        raise AttributeError(msg)


class _AutodocAttrGetter:
    """getattr() override for types such as Zope interfaces."""

    _attr_getters: Sequence[tuple[type, Callable[[Any, str, Any], Any]]]

    __slots__ = ('_attr_getters',)

    def __init__(
        self, attr_getters: dict[type, Callable[[Any, str, Any], Any]], /
    ) -> None:
        super().__setattr__('_attr_getters', tuple(attr_getters.items()))

    def __call__(self, obj: Any, name: str, *defargs: Any) -> Any:
        for typ, func in self._attr_getters:
            if isinstance(obj, typ):
                return func(obj, name, *defargs)

        return safe_getattr(obj, name, *defargs)

    def __repr__(self) -> str:
        return f'_AutodocAttrGetter({dict(self._attr_getters)!r})'

    def __setattr__(self, key: str, value: Any) -> NoReturn:
        msg = f'{self.__class__.__name__} is immutable'
        raise AttributeError(msg)

    def __delattr__(self, key: str) -> NoReturn:
        msg = f'{self.__class__.__name__} is immutable'
        raise AttributeError(msg)


def _get_render_mode(
    typehints_format: Literal['fully-qualified', 'short'],
    /,
) -> _RestifyMode:
    if typehints_format == 'short':
        return 'smart'
    return 'fully-qualified-except-typing'
