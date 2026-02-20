"""Signature utilities for autodoc"""

from __future__ import annotations

import sys
from inspect import Parameter, Signature
from typing import TYPE_CHECKING, NewType, TypeVar

from sphinx.errors import PycodeError
from sphinx.ext.autodoc._dynamic._preserve_defaults import update_default_value
from sphinx.ext.autodoc._dynamic._type_annotations import _record_typehints
from sphinx.ext.autodoc._dynamic._type_comments import (
    _update_annotations_using_type_comments,
)
from sphinx.ext.autodoc._names import py_ext_sig_re
from sphinx.ext.autodoc._property_types import _AssignStatementProperties
from sphinx.ext.autodoc._shared import LOGGER
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import inspect
from sphinx.util.docstrings import prepare_docstring
from sphinx.util.inspect import (
    _stringify_signature_to_parts,
    evaluate_signature,
    safe_getattr,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping
    from typing import Any

    from sphinx.events import EventManager
    from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
    from sphinx.ext.autodoc._property_types import _ItemProperties
    from sphinx.ext.autodoc._shared import _AttrGetter, _AutodocConfig

    type _FormattedSignature = tuple[str, str]


def _format_signatures(
    *,
    autodoc_annotations: dict[str, dict[str, str]],
    config: _AutodocConfig,
    docstrings: list[list[str]] | None,
    events: EventManager,
    get_attr: _AttrGetter,
    parent: Any,
    options: _AutoDocumenterOptions,
    props: _ItemProperties,
    args: str | None = None,
    retann: str | None = '',
    **kwargs: Any,
) -> list[_FormattedSignature]:
    """Format the signature (arguments and return annotation) of the object.

    Let the user process it via the ``autodoc-process-signature`` event.
    """
    if props.obj_type in {'class', 'exception'}:
        from sphinx.ext.autodoc._property_types import _ClassDefProperties

        assert isinstance(props, _ClassDefProperties)
        if props.doc_as_attr:
            return []
        if config.autodoc_class_signature == 'separated':
            # do not show signatures
            return []

    if config.autodoc_typehints_format == 'short':
        kwargs.setdefault('unqualified_typehints', True)
    if config.python_display_short_literal_types:
        kwargs.setdefault('short_literals', True)

    if args is None:
        signatures: list[_FormattedSignature] = []
    else:
        signatures = [(args, retann or '')]

    if (
        not signatures
        and config.autodoc_docstring_signature
        and props.obj_type not in {'module', 'data', 'type'}
        and docstrings is not None
    ):
        # only act if a signature is not explicitly given already,
        # and if the feature is enabled
        signatures[:] = _extract_signatures_from_docstrings(
            docstrings, props=props, tab_width=options._tab_width
        )

    if not signatures:
        # try to introspect the signature
        try:
            signatures[:] = _extract_signature_from_object(
                config=config,
                events=events,
                get_attr=get_attr,
                parent=parent,
                props=props,
                **kwargs,
            )
        except Exception as exc:
            msg = __('error while formatting arguments for %s: %s')
            LOGGER.warning(msg, props.full_name, exc, type='autodoc')

    if props.obj_type in {'attribute', 'property'}:
        # Only keep the return annotation
        signatures = [('', retann) for _args, retann in signatures]

    _record_typehints(
        autodoc_annotations=autodoc_annotations,
        name=props.full_name,
        obj=props._obj,
        short_literals=kwargs.get('short_literals', False),
        type_aliases=config.autodoc_type_aliases,
        unqualified_typehints=kwargs.get('unqualified_typehints', False),
    )
    if result := events.emit_firstresult(
        'autodoc-process-signature',
        props.obj_type,
        props.full_name,
        props._obj,
        options,
        signatures[0][0] if signatures else None,  # args
        signatures[0][1] if signatures else '',  # retann
    ):
        if len(result) == 2 and isinstance(result[0], str):
            args, retann = result
            signatures[0] = (args, retann if isinstance(retann, str) else '')

    if props.obj_type in {'module', 'data', 'type'}:
        signatures[1:] = ()  # discard all signatures save the first

    analyzer_overloads: dict[str, list[Signature]] = {}
    try:
        analyzer = ModuleAnalyzer.for_module(props.canonical_module_name)
        # parse right now, to get PycodeErrors on parsing (results will
        # be cached anyway)
        analyzer.analyze()
    except PycodeError as exc:
        LOGGER.debug('[autodoc] module analyzer failed: %s', exc)
        # no source file -- e.g. for builtin and C modules
    else:
        analyzer_overloads = analyzer.overloads

    if props.obj_type in {'function', 'decorator'}:
        overloaded = (
            props.dotted_parts in analyzer_overloads
            and config.autodoc_typehints != 'none'
        )
        is_singledispatch = inspect.is_singledispatch_function(props._obj)

        if overloaded:
            # Use signatures for overloaded functions and methods instead of
            # their implementations.
            signatures.clear()
        elif not is_singledispatch:
            return signatures

        if is_singledispatch:
            from sphinx.ext.autodoc._property_types import _FunctionDefProperties

            # append signature of singledispatch'ed functions
            for typ, func in props._obj.registry.items():
                if typ is object:
                    continue  # default implementation. skipped.
                dispatch_func = _annotate_to_first_argument(
                    func, typ, config=config, props=props
                )
                if not dispatch_func:
                    continue
                dispatch_props = _FunctionDefProperties(
                    obj_type='function',
                    module_name='',
                    parts=('',),
                    docstring_lines=(),
                    signatures=(),
                    _obj=dispatch_func,
                    _obj___module__=None,
                    _obj___qualname__=None,
                    _obj___name__=None,
                    properties=frozenset(),
                )
                signatures += _format_signatures(
                    autodoc_annotations=autodoc_annotations,
                    config=config,
                    docstrings=None,
                    events=events,
                    get_attr=get_attr,
                    parent=None,
                    options=options,
                    props=dispatch_props,
                )
        if overloaded:
            actual = inspect.signature(
                props._obj, type_aliases=config.autodoc_type_aliases
            )
            obj_globals = safe_getattr(props._obj, '__globals__', {})
            overloads = analyzer_overloads[props.dotted_parts]
            for overload in overloads:
                overload = _merge_default_value(actual, overload)
                overload = evaluate_signature(
                    overload, obj_globals, config.autodoc_type_aliases
                )
                signatures.append(_stringify_signature_to_parts(overload, **kwargs))

        return signatures

    if props.obj_type in {'class', 'exception'}:
        from sphinx.ext.autodoc._property_types import _ClassDefProperties

        assert isinstance(props, _ClassDefProperties)
        method_name = props._signature_method_name
        if method_name == '__call__':
            signature_cls = type(props._obj)
        else:
            signature_cls = props._obj
        overloads = []
        overloaded = False
        if method_name:
            for cls in signature_cls.__mro__:
                try:
                    analyzer = ModuleAnalyzer.for_module(cls.__module__)
                    analyzer.analyze()
                except PycodeError:
                    pass
                else:
                    qualname = f'{cls.__qualname__}.{method_name}'
                    if qualname in analyzer.overloads:
                        overloads = analyzer.overloads[qualname]
                        overloaded = True
                        break
                    if qualname in analyzer.tagorder:
                        # the constructor is defined in the class, but not overridden.
                        break
        if overloaded and config.autodoc_typehints != 'none':
            # Use signatures for overloaded methods instead of the implementation method.
            signatures.clear()
            method = safe_getattr(signature_cls, method_name, None)
            method_globals = safe_getattr(method, '__globals__', {})
            for overload in overloads:
                overload = evaluate_signature(
                    overload, method_globals, config.autodoc_type_aliases
                )

                parameters = list(overload.parameters.values())
                overload = overload.replace(
                    parameters=parameters[1:], return_annotation=Parameter.empty
                )
                signatures.append(_stringify_signature_to_parts(overload, **kwargs))
            return signatures

        return signatures

    if props.obj_type == 'method':
        overloaded = (
            props.dotted_parts in analyzer_overloads
            and config.autodoc_typehints != 'none'
        )
        meth = parent.__dict__.get(props.name)
        is_singledispatch = inspect.is_singledispatch_method(meth)

        if overloaded:
            # Use signatures for overloaded functions and methods instead of
            # their implementations.
            signatures.clear()
        elif not is_singledispatch:
            return signatures

        if is_singledispatch:
            from sphinx.ext.autodoc._property_types import _FunctionDefProperties

            # append signature of singledispatch'ed methods
            for typ, func in meth.dispatcher.registry.items():
                if typ is object:
                    continue  # default implementation. skipped.
                if inspect.isclassmethod(func):
                    func = func.__func__
                dispatch_meth = _annotate_to_first_argument(
                    func, typ, config=config, props=props
                )
                if not dispatch_meth:
                    continue
                dispatch_props = _FunctionDefProperties(
                    obj_type='method',
                    module_name='',
                    parts=('',),
                    docstring_lines=(),
                    signatures=(),
                    _obj=dispatch_meth,
                    _obj___module__=None,
                    _obj___qualname__=None,
                    _obj___name__=None,
                    properties=frozenset(),
                )
                signatures += _format_signatures(
                    autodoc_annotations=autodoc_annotations,
                    config=config,
                    docstrings=None,
                    events=events,
                    get_attr=get_attr,
                    parent=parent,
                    options=options,
                    props=dispatch_props,
                )
        if overloaded:
            from sphinx.ext.autodoc._property_types import _FunctionDefProperties

            assert isinstance(props, _FunctionDefProperties)
            actual = inspect.signature(
                props._obj,
                bound_method=not props.is_staticmethod,
                type_aliases=config.autodoc_type_aliases,
            )

            obj_globals = safe_getattr(props._obj, '__globals__', {})
            overloads = analyzer_overloads[props.dotted_parts]
            for overload in overloads:
                overload = _merge_default_value(actual, overload)
                overload = evaluate_signature(
                    overload, obj_globals, config.autodoc_type_aliases
                )

                if not props.is_staticmethod:
                    # hide the first argument (e.g. 'self')
                    parameters = list(overload.parameters.values())
                    overload = overload.replace(parameters=parameters[1:])
                signatures.append(_stringify_signature_to_parts(overload, **kwargs))

        return signatures

    return signatures


def _extract_signatures_from_docstrings(
    docstrings: list[list[str]],
    /,
    props: _ItemProperties,
    tab_width: int,
) -> list[_FormattedSignature]:
    signatures: list[_FormattedSignature] = []

    # candidates of the object name
    valid_names = {props.name}
    if props.obj_type in {'class', 'exception'}:
        valid_names.add('__init__')
        if hasattr(props._obj, '__mro__'):
            valid_names |= {cls.__name__ for cls in props._obj.__mro__}

    stripped_docstrings = [list(l) for l in (docstrings or ())]
    for i, doclines in enumerate(docstrings):
        j = 0
        for j, line in enumerate(doclines):  # NoQA: B007
            if not line:
                # no lines in docstring, no match
                break
            line = line.rstrip('\\').rstrip()

            # match first line of docstring against signature RE
            match = py_ext_sig_re.match(line)
            if not match:
                break
            _exmod, _path, base, _tp_list, args, retann = match.groups()
            if args is not None:
                args = f'({args})'
            else:
                args = ''  # i.e. property or attribute

            # the base name must match ours
            if base not in valid_names:
                break

            if props.obj_type in {'class', 'exception'} and retann == 'None':
                # Strip a return value from signatures of constructor in docstring
                signatures.append((args, ''))
            else:
                signatures.append((args, retann or ''))

        if signatures:
            # re-prepare docstring to ignore more leading indentation
            stripped_docstrings[i] = prepare_docstring(
                '\n'.join(doclines[j:]), tab_width
            )

            # finish the loop after finding at least one signature
            break

    if not signatures:
        return []

    # Update docstrings from stripped_docstrings if needed
    if props.obj_type in {
        'class',
        'exception',
        'function',
        'method',
        'property',
        'decorator',
    } or (
        props.obj_type == 'attribute'
        and isinstance(props, _AssignStatementProperties)
        and props._obj_is_attribute_descriptor
    ):
        docstrings[:] = stripped_docstrings

    return signatures


def _extract_signature_from_object(
    config: _AutodocConfig,
    events: EventManager,
    get_attr: _AttrGetter,
    parent: Any,
    props: _ItemProperties,
    **kwargs: Any,
) -> list[_FormattedSignature]:
    """Format the signature using runtime introspection."""
    sig = _get_signature_object(
        events=events,
        get_attr=get_attr,
        parent=parent,
        preserve_defaults=config.autodoc_preserve_defaults,
        props=props,
        type_aliases=config.autodoc_type_aliases,
        use_type_comments=config.autodoc_use_type_comments,
    )
    if sig is None:
        return []

    if props.obj_type == 'decorator' and len(sig.parameters) == 1:
        # Special case for single-argument decorators
        return [('', '')]

    if config.autodoc_typehints in {'none', 'description'}:
        kwargs.setdefault('show_annotation', False)
    if config.autodoc_typehints_format == 'short':
        kwargs.setdefault('unqualified_typehints', True)
    if config.python_display_short_literal_types:
        kwargs.setdefault('short_literals', True)
    if props.obj_type in {'class', 'exception'}:
        kwargs['show_return_annotation'] = False

    args, retann = _stringify_signature_to_parts(sig, **kwargs)
    if config.strip_signature_backslash:
        # escape backslashes for reST
        args = args.replace('\\', '\\\\')
        retann = retann.replace('\\', '\\\\')

    return [(args, retann)]


# Types which have confusing metaclass signatures it would be best not to show.
# These are listed by name, rather than storing the objects themselves, to avoid
# needing to import the modules.
_METACLASS_CALL_BLACKLIST = frozenset({
    'enum.EnumType.__call__',
})


# Types whose __new__ signature is a pass-through.
_CLASS_NEW_BLACKLIST = frozenset({
    'typing.Generic.__new__',
})


def _get_signature_object(
    events: EventManager,
    get_attr: _AttrGetter,
    parent: Any,
    preserve_defaults: bool,
    props: _ItemProperties,
    type_aliases: Mapping[str, str] | None,
    use_type_comments: bool,
) -> Signature | None:
    """Return a Signature for *obj*, or None on failure."""
    obj, is_bound_method = _get_object_for_signature(
        props=props, get_attr=get_attr, parent=parent, type_aliases=type_aliases
    )
    if obj is None or isinstance(obj, Signature):
        return obj

    if preserve_defaults:
        update_default_value(obj, bound_method=is_bound_method)
    if use_type_comments:
        _update_annotations_using_type_comments(obj, bound_method=is_bound_method)
    events.emit('autodoc-before-process-signature', obj, is_bound_method)

    if props.obj_type in {'class', 'exception', 'function', 'method', 'decorator'}:
        try:
            return inspect.signature(
                obj, bound_method=is_bound_method, type_aliases=type_aliases
            )
        except TypeError as exc:
            if props.obj_type in {'class', 'exception'}:
                msg = __('Failed to get a constructor signature for %s: %s')
            elif props.obj_type in {'function', 'decorator'}:
                msg = __('Failed to get a function signature for %s: %s')
            elif props.obj_type == 'method':
                msg = __('Failed to get a method signature for %s: %s')
            else:
                msg = __('Failed to get a signature for %s: %s')
            LOGGER.warning(msg, props.full_name, exc)
            return None
        except ValueError:
            # Still no signature: happens e.g. for old-style classes
            # with __init__ in C and no `__text_signature__`.
            return None

    return None


def _get_object_for_signature(
    props: _ItemProperties,
    get_attr: _AttrGetter,
    parent: Any,
    type_aliases: Mapping[str, str] | None,
) -> tuple[Any, bool]:
    """Return the object from which we will obtain the signature."""
    obj = props._obj
    if props.obj_type in {'function', 'decorator'}:
        return obj, False

    if props.obj_type in {'class', 'exception'}:
        if isinstance(obj, (NewType, TypeVar)):
            # Suppress signature
            return None, False

        try:
            object_sig = obj.__signature__
        except AttributeError:
            pass
        else:
            if isinstance(object_sig, Signature):
                return object_sig, False
            if sys.version_info[:2] <= (3, 14) and callable(object_sig):
                # Support for enum.Enum.__signature__ in Python 3.12 & 3.13
                if isinstance(object_sig_str := object_sig(), str):
                    return inspect.signature_from_str(object_sig_str), False

        def get_user_defined_function_or_method(obj: Any, attr: str) -> Any:
            """Get the `attr` function or method from `obj`, if it is user-defined."""
            if inspect.is_builtin_class_method(obj, attr):
                return None
            attr = get_attr(obj, attr, None)
            if not (inspect.ismethod(attr) or inspect.isfunction(attr)):
                return None
            return attr

        # This sequence is copied from inspect._signature_from_callable.
        # ValueError means that no signature could be found, so we keep going.

        # Let's see if it has an overloaded __call__ defined in its metaclass,
        # or if the 'obj' class has a '__new__' or '__init__' method
        for obj_, meth_name, blacklist in (
            (type(obj), '__call__', _METACLASS_CALL_BLACKLIST),
            (obj, '__new__', _CLASS_NEW_BLACKLIST),
            (obj, '__init__', frozenset()),
        ):
            meth = get_user_defined_function_or_method(obj_, meth_name)
            if meth is None:
                continue
            if blacklist:
                if f'{meth.__module__}.{meth.__qualname__}' in blacklist:
                    continue

            try:
                inspect.signature(meth, bound_method=True, type_aliases=type_aliases)
            except TypeError:
                return meth, True  # _get_signature_object() needs to log the failure
            except ValueError:
                continue
            else:
                from sphinx.ext.autodoc._property_types import _ClassDefProperties

                assert isinstance(props, _ClassDefProperties)
                props._signature_method_name = meth_name
                return meth, True

        # None of the attributes are user-defined, so fall back to let inspect
        # handle it.
        # We don't know the exact method that inspect.signature will read
        # the signature from, so just return the object itself to be passed
        # to the ``autodoc-before-process-signature`` hook.
        return obj, False

    if props.obj_type == 'method':
        if obj == object.__init__ and parent != object:  # NoQA: E721
            # Classes not having own __init__() method are shown as no arguments.
            #
            # Note: The signature of object.__init__() is (self, /, *args, **kwargs).
            #       But it makes users confused.
            return Signature(), False

        is_bound_method = not inspect.isstaticmethod(
            obj, cls=parent, name=props.object_name
        )
        return obj, is_bound_method

    return None, False


def _annotate_to_first_argument(
    func: Callable[..., Any],
    typ: type,
    *,
    config: _AutodocConfig,
    props: _ItemProperties,
) -> Callable[..., Any] | None:
    """Annotate type hint to the first argument of function if needed."""
    try:
        sig = inspect.signature(func, type_aliases=config.autodoc_type_aliases)
    except TypeError as exc:
        msg = __('Failed to get a function signature for %s: %s')
        LOGGER.warning(msg, props.full_name, exc)
        return None
    except ValueError:
        return None

    first_arg_idx = 1 * (props.obj_type == 'method')
    if len(sig.parameters) == first_arg_idx:
        return None

    def dummy():  # type: ignore[no-untyped-def]  # NoQA: ANN202
        pass

    params = list(sig.parameters.values())
    if params[first_arg_idx].annotation is Parameter.empty:
        params[first_arg_idx] = params[first_arg_idx].replace(annotation=typ)
        try:
            dummy.__signature__ = sig.replace(parameters=params)  # type: ignore[attr-defined]
            return dummy
        except (AttributeError, TypeError):
            # failed to update signature (ex. built-in or extension types)
            return None

    return func


def _merge_default_value(actual: Signature, overload: Signature) -> Signature:
    """Merge default values of actual implementation to the overload variants."""
    parameters = list(overload.parameters.values())
    for i, param in enumerate(parameters):
        actual_param = actual.parameters.get(param.name)
        if actual_param and param.default == '...':
            parameters[i] = param.replace(default=actual_param.default)

    return overload.replace(parameters=parameters)
