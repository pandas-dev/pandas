from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from docutils.utils import assemble_option_dict

from sphinx.deprecation import RemovedInSphinx11Warning
from sphinx.ext.autodoc._sentinels import ALL, EMPTY, SUPPRESS
from sphinx.locale import __

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Set
    from typing import Any, Final, Literal, Self

    from sphinx.ext.autodoc._property_types import _AutodocObjType
    from sphinx.ext.autodoc._sentinels import ALL_T, EMPTY_T, SUPPRESS_T
    from sphinx.util.typing import OptionSpec


# common option names for autodoc directives
AUTODOC_DEFAULT_OPTIONS = (
    'members',
    'undoc-members',
    'no-index',
    'no-index-entry',
    'inherited-members',
    'show-inheritance',
    'private-members',
    'special-members',
    'ignore-module-all',
    'exclude-members',
    'member-order',
    'imported-members',
    'class-doc-from',
    'no-value',
)

AUTODOC_EXTENDABLE_OPTIONS = frozenset({
    'members',
    'private-members',
    'special-members',
    'exclude-members',
})


class _AutoDocumenterOptions:
    # TODO: make immutable.

    no_index: Literal[True] | None = None
    no_index_entry: Literal[True] | None = None
    _tab_width: int = 8

    # module-like options
    members: ALL_T | list[str] | None = None
    undoc_members: Literal[True] | None = None
    inherited_members: Set[str] | None = None
    show_inheritance: Literal[True] | None = None
    synopsis: str | None = None
    platform: str | None = None
    deprecated: Literal[True] | None = None
    member_order: Literal['alphabetical', 'bysource', 'groupwise'] | None = None
    exclude_members: EMPTY_T | set[str] | None = None
    private_members: ALL_T | list[str] | None = None
    special_members: ALL_T | list[str] | None = None
    imported_members: Literal[True] | None = None
    ignore_module_all: Literal[True] | None = None
    no_value: Literal[True] | None = None

    # class-like options (class, exception)
    class_doc_from: Literal['both', 'class', 'init'] | None = None

    # assignment-like (data, attribute)
    annotation: SUPPRESS_T | str | None = None

    noindex: Literal[True] | None = None

    def __init__(self, **kwargs: Any) -> None:
        vars(self).update(kwargs)

    def __repr__(self) -> str:
        args = ', '.join(f'{k}={v!r}' for k, v in vars(self).items())
        return f'_AutoDocumenterOptions({args})'

    def __getattr__(self, name: str) -> object:
        return None  # return None for missing attributes

    def copy(self) -> Self:
        return self.__class__(**vars(self))

    @classmethod
    def from_directive_options(cls, opts: Mapping[str, Any], /) -> Self:
        return cls(**{k.replace('-', '_'): v for k, v in opts.items() if v is not None})

    # Mapping interface:

    def __getitem__(self, item: str) -> Any:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        try:
            return getattr(self, item)
        except AttributeError:
            raise KeyError(item) from None

    def __setitem__(self, key: str, value: Any) -> None:
        msg = f'{self.__class__.__name__!r} object does not support indexed assignment'
        raise TypeError(msg)

    def __delitem__(self, key: str) -> None:
        msg = f'{self.__class__.__name__!r} object does not support indexed deletion'
        raise TypeError(msg)

    def __contains__(self, item: str) -> bool:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        return hasattr(self, item)

    def __keys(self) -> list[str]:
        return [key for key in dir(self) if not key.startswith('_')]

    def __iter__(self) -> Iterator[str]:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        yield from self.__keys()

    def __len__(self) -> int:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        return len(self.__keys())

    def keys(self) -> Iterable[str]:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        yield from self.__keys()

    def items(self) -> Iterable[tuple[str, Any]]:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        for key in self.__keys():
            yield key, getattr(self, key)

    def values(self) -> Iterable[Any]:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        for key in self.__keys():
            yield getattr(self, key)

    def get(self, key: str, default: Any | None = None) -> Any | None:
        warnings.warn(
            'The mapping interface for autodoc options objects is deprecated, '
            'and will be removed in Sphinx 11. Use attribute access instead.',
            RemovedInSphinx11Warning,
            stacklevel=2,
        )
        try:
            return getattr(self, key)
        except AttributeError:
            return default


def identity(x: Any) -> Any:
    return x


def members_option(arg: str | None) -> ALL_T | list[str] | None:
    """Used to convert the :members: option to auto directives."""
    if arg is None or arg is True:
        return ALL
    if arg is False:
        return None
    return [stripped for x in arg.split(',') if (stripped := x.strip())]


def exclude_members_option(arg: str | None) -> EMPTY_T | set[str]:
    """Used to convert the :exclude-members: option."""
    if arg is None or arg is True:
        return EMPTY
    return {stripped for x in arg.split(',') if (stripped := x.strip())}


def inherited_members_option(arg: str | None) -> set[str]:
    """Used to convert the :inherited-members: option to auto directives."""
    if arg is None or arg is True:
        return {'object'}
    if arg:
        return {x.strip() for x in arg.split(',')}
    return set()


def member_order_option(
    arg: str | None,
) -> Literal['alphabetical', 'bysource', 'groupwise'] | None:
    """Used to convert the :member-order: option to auto directives."""
    if arg is None or arg is True:
        return None
    if arg in {'alphabetical', 'bysource', 'groupwise'}:
        return arg  # type: ignore[return-value]
    raise ValueError(__('invalid value for member-order option: %s') % arg)


def class_doc_from_option(arg: str | None) -> Literal['both', 'class', 'init']:
    """Used to convert the :class-doc-from: option to autoclass directives."""
    if arg in {'both', 'class', 'init'}:
        return arg  # type: ignore[return-value]
    raise ValueError(__('invalid value for class-doc-from option: %s') % arg)


def annotation_option(arg: str | None) -> SUPPRESS_T | str | Literal[False]:
    if arg is None or arg is True:
        # suppress showing the representation of the object
        return SUPPRESS
    return arg


def bool_option(arg: str | None) -> bool:
    """Used to convert flag options to auto directives.  (Instead of
    directives.flag(), which returns None).
    """
    return True


def merge_members_option(options: dict[str, Any]) -> None:
    """Merge :private-members: and :special-members: options to the
    :members: option.
    """
    if options.get('members') is ALL:
        # merging is not needed when members: ALL
        return

    members = options.setdefault('members', [])
    for key in ('private-members', 'special-members'):
        other_members = options.get(key)
        if other_members is not None and other_members is not ALL:
            for member in other_members:
                if member not in members:
                    members.append(member)


_OPTION_SPEC_COMMON: Final[OptionSpec] = {
    'no-index': bool_option,
    'no-index-entry': bool_option,
}
_OPTION_SPEC_HAS_MEMBERS: Final[OptionSpec] = _OPTION_SPEC_COMMON | {
    'members': members_option,
    'exclude-members': exclude_members_option,
    'undoc-members': bool_option,
    'private-members': members_option,
    'special-members': members_option,
    'member-order': member_order_option,
}
_OPTION_SPEC_MODULE_SPECIFIC: Final[OptionSpec] = {
    'ignore-module-all': bool_option,
    'imported-members': bool_option,
    'deprecated': bool_option,
    'platform': identity,
    'synopsis': identity,
}
_OPTION_SPEC_CLASS_SPECIFIC: Final[OptionSpec] = {
    'class-doc-from': class_doc_from_option,
    'show-inheritance': bool_option,
    'inherited-members': inherited_members_option,
}
_OPTION_SPEC_ASSIGNMENT: Final[OptionSpec] = _OPTION_SPEC_COMMON | {
    'annotation': annotation_option,
    'no-value': bool_option,
}
_OPTION_SPEC_DEPRECATED: Final[OptionSpec] = {
    'noindex': bool_option,
}
_OPTION_SPEC_FUNCTION_DEF: Final = _OPTION_SPEC_COMMON | _OPTION_SPEC_DEPRECATED
_OPTION_SPECS: Final[Mapping[_AutodocObjType, OptionSpec]] = {
    'module': _OPTION_SPEC_HAS_MEMBERS
    | _OPTION_SPEC_MODULE_SPECIFIC
    | {'show-inheritance': bool_option}  # special case
    | {'inherited-members': inherited_members_option}  # special case
    | {'no-value': bool_option}  # special case
    | _OPTION_SPEC_DEPRECATED,
    'class': _OPTION_SPEC_HAS_MEMBERS
    | _OPTION_SPEC_CLASS_SPECIFIC
    | _OPTION_SPEC_DEPRECATED,
    'exception': _OPTION_SPEC_HAS_MEMBERS
    | _OPTION_SPEC_CLASS_SPECIFIC
    | _OPTION_SPEC_DEPRECATED,
    'function': _OPTION_SPEC_FUNCTION_DEF,
    'decorator': _OPTION_SPEC_FUNCTION_DEF,
    'method': _OPTION_SPEC_FUNCTION_DEF,
    'property': _OPTION_SPEC_FUNCTION_DEF,
    'attribute': _OPTION_SPEC_ASSIGNMENT | _OPTION_SPEC_DEPRECATED,
    'data': _OPTION_SPEC_ASSIGNMENT | _OPTION_SPEC_DEPRECATED,
    'type': _OPTION_SPEC_ASSIGNMENT,
}


def _process_documenter_options(
    *,
    obj_type: _AutodocObjType,
    default_options: Mapping[str, str | bool],
    options: dict[str, str | None],
) -> _AutoDocumenterOptions:
    """Recognize options of object type from user input."""
    option_spec = _OPTION_SPECS[obj_type]
    for name in AUTODOC_DEFAULT_OPTIONS:
        if name not in option_spec:
            continue

        negated = options.pop(f'no-{name}', True) is None
        if name in default_options and not negated:
            if name in options and isinstance(default_options[name], str):
                # take value from options if present or extend it
                # with autodoc_default_options if necessary
                if name in AUTODOC_EXTENDABLE_OPTIONS:
                    opt_value = options[name]
                    if opt_value is not None and opt_value.startswith('+'):
                        options[name] = f'{default_options[name]},{opt_value[1:]}'
            else:
                options[name] = default_options[name]  # type: ignore[assignment]
        elif (opt_value := options.get(name)) is not None:
            # remove '+' from option argument if there's nothing to merge it with
            options[name] = opt_value.removeprefix('+')

    opts = assemble_option_dict(options.items(), option_spec)  # type: ignore[arg-type]
    return _AutoDocumenterOptions.from_directive_options(opts)
