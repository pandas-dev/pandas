"""Update annotations info of living objects using type_comments."""

from __future__ import annotations

import ast
import contextlib
import sys
from inspect import Parameter, Signature, getsource
from types import ModuleType
from typing import TYPE_CHECKING, cast
from weakref import WeakSet

from sphinx.errors import PycodeError
from sphinx.ext.autodoc._shared import LOGGER
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer
from sphinx.pycode.ast import unparse as ast_unparse
from sphinx.util import inspect
from sphinx.util.inspect import safe_getattr

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any


_objects_with_type_comment_annotations: WeakSet[Any] = WeakSet()
"""Cache of objects with annotations updated from type comments."""


def _ensure_annotations_from_type_comments(obj: Any) -> None:
    """Ensures `obj.__annotations__` includes type comment information.

    Failures to assign to `__annotations__` are silently ignored.

    If `obj` is a class type, this also ensures that type comment
    information is incorporated into the `__annotations__` member of
    all parent classes, if possible.

    This mutates the `__annotations__` of existing imported objects,
    in order to allow the existing `typing.get_type_hints` method to
    take the modified annotations into account.

    Modifying existing imported objects is unfortunate but avoids the
    need to reimplement `typing.get_type_hints` in order to take into
    account type comment information.

    Note that this does not directly include type comment information
    from parent classes, but `typing.get_type_hints` takes that into
    account.
    """
    if obj in _objects_with_type_comment_annotations:
        return
    with contextlib.suppress(TypeError):
        _objects_with_type_comment_annotations.add(obj)

    if isinstance(obj, type):
        for cls in inspect.getmro(obj):
            modname = safe_getattr(cls, '__module__')
            mod = sys.modules.get(modname)
            if mod is not None:
                _ensure_annotations_from_type_comments(mod)

    elif isinstance(obj, ModuleType):
        _update_module_annotations_from_type_comments(obj)


def _update_module_annotations_from_type_comments(mod: ModuleType) -> None:
    """Adds type comment annotations for a single module.

    Both module-level and class-level annotations are added.
    """
    mod_annotations = dict(inspect.getannotations(mod))
    mod.__annotations__ = mod_annotations

    class_annotations: dict[str, dict[str, Any]] = {}

    try:
        analyzer = ModuleAnalyzer.for_module(mod.__name__)
        analyzer.analyze()
        anns = analyzer.annotations
        for (classname, attrname), annotation in anns.items():
            if not classname:
                annotations = mod_annotations
            else:
                cls_annotations = class_annotations.get(classname)
                if cls_annotations is None:
                    try:
                        cls = mod
                        for part in classname.split('.'):
                            cls = safe_getattr(cls, part)
                        annotations = dict(inspect.getannotations(cls))
                        # Ignore errors setting __annotations__
                        with contextlib.suppress(TypeError, AttributeError):
                            cls.__annotations__ = annotations
                    except AttributeError:
                        annotations = {}
                    class_annotations[classname] = annotations
                else:
                    annotations = cls_annotations
            annotations.setdefault(attrname, annotation)
    except PycodeError:
        pass


def _update_annotations_using_type_comments(obj: Any, bound_method: bool) -> None:
    """Update annotations info of *obj* using type_comments."""
    try:
        type_sig = get_type_comment(obj, bound_method)
        if type_sig:
            sig = inspect.signature(obj, bound_method)
            for param in sig.parameters.values():
                if param.name not in obj.__annotations__:
                    annotation = type_sig.parameters[param.name].annotation
                    if annotation is not Parameter.empty:
                        obj.__annotations__[param.name] = ast_unparse(annotation)

            if 'return' not in obj.__annotations__:
                obj.__annotations__['return'] = type_sig.return_annotation
    except KeyError as exc:
        LOGGER.warning(
            __('Failed to update signature for %r: parameter not found: %s'), obj, exc
        )
    except NotImplementedError as exc:  # failed to ast.unparse()
        LOGGER.warning(__('Failed to parse type_comment for %r: %s'), obj, exc)


def get_type_comment(obj: Any, bound_method: bool = False) -> Signature | None:
    """Get type_comment'ed FunctionDef object from living object.

    This tries to parse original code for living object and returns
    Signature for given *obj*.
    """
    try:
        source = getsource(obj)
        if source.startswith((' ', r'\t')):
            # subject is placed inside class or block.  To read its docstring,
            # this adds if-block before the declaration.
            module = ast.parse('if True:\n' + source, type_comments=True)
            subject = cast('ast.FunctionDef', module.body[0].body[0])  # type: ignore[attr-defined]
        else:
            module = ast.parse(source, type_comments=True)
            subject = cast('ast.FunctionDef', module.body[0])

        type_comment = getattr(subject, 'type_comment', None)
        if type_comment:
            function = ast.parse(type_comment, mode='func_type', type_comments=True)
            return signature_from_ast(subject, bound_method, function)
        else:
            return None
    except (OSError, TypeError):  # failed to load source code
        return None
    except SyntaxError:  # failed to parse type_comments
        return None


def signature_from_ast(
    node: ast.FunctionDef, bound_method: bool, type_comment: ast.FunctionDef
) -> Signature:
    """Return a Signature object for the given *node*.

    :param bound_method: Specify *node* is a bound method or not
    """
    params = []
    for arg in node.args.posonlyargs:
        param = Parameter(
            arg.arg,
            Parameter.POSITIONAL_ONLY,
            annotation=arg.type_comment,
        )
        params.append(param)

    for arg in node.args.args:
        param = Parameter(
            arg.arg,
            Parameter.POSITIONAL_OR_KEYWORD,
            annotation=arg.type_comment or Parameter.empty,
        )
        params.append(param)

    if node.args.vararg:
        param = Parameter(
            node.args.vararg.arg,
            Parameter.VAR_POSITIONAL,
            annotation=node.args.vararg.type_comment or Parameter.empty,
        )
        params.append(param)

    for arg in node.args.kwonlyargs:
        param = Parameter(
            arg.arg,
            Parameter.KEYWORD_ONLY,
            annotation=arg.type_comment or Parameter.empty,
        )
        params.append(param)

    if node.args.kwarg:
        param = Parameter(
            node.args.kwarg.arg,
            Parameter.VAR_KEYWORD,
            annotation=node.args.kwarg.type_comment or Parameter.empty,
        )
        params.append(param)

    # Remove first parameter when *obj* is bound_method
    if bound_method and params:
        params.pop(0)

    # merge type_comment into signature
    if not_suppressed(type_comment.argtypes):  # type: ignore[attr-defined]
        for i, param in enumerate(params):
            params[i] = param.replace(annotation=type_comment.argtypes[i])  # type: ignore[attr-defined]

    if node.returns:
        return Signature(params, return_annotation=node.returns)
    elif type_comment.returns:
        return Signature(params, return_annotation=ast_unparse(type_comment.returns))
    else:
        return Signature(params)


def not_suppressed(argtypes: Sequence[ast.expr] = ()) -> bool:
    """Check given *argtypes* is suppressed type_comment or not."""
    if len(argtypes) == 0:  # no argtypees
        return False
    if len(argtypes) == 1:
        arg = argtypes[0]
        if isinstance(arg, ast.Constant) and arg.value is ...:  # suppressed
            return False
    # not suppressed
    return True
