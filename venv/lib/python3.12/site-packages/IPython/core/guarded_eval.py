from copy import copy
from inspect import isclass, signature, Signature, getmodule
from typing import (
    Annotated,
    AnyStr,
    Literal,
    NamedTuple,
    NewType,
    Optional,
    Protocol,
    TypeGuard,
    Union,
    get_args,
    get_origin,
    is_typeddict,
)
from collections.abc import Callable, Sequence
import ast
import builtins
import collections
import dataclasses
import operator
import sys
import typing
import warnings
from functools import cached_property
from dataclasses import dataclass, field
from types import MethodDescriptorType, ModuleType, MethodType

from IPython.utils.decorators import undoc

import types
from typing import Self, LiteralString, get_type_hints

if sys.version_info < (3, 12):
    from typing_extensions import TypeAliasType
else:
    from typing import TypeAliasType


@undoc
class HasGetItem(Protocol):
    def __getitem__(self, key) -> None:
        ...


@undoc
class InstancesHaveGetItem(Protocol):
    def __call__(self, *args, **kwargs) -> HasGetItem:
        ...


@undoc
class HasGetAttr(Protocol):
    def __getattr__(self, key) -> None:
        ...


@undoc
class DoesNotHaveGetAttr(Protocol):
    pass


# By default `__getattr__` is not explicitly implemented on most objects
MayHaveGetattr = Union[HasGetAttr, DoesNotHaveGetAttr]


def _unbind_method(func: Callable) -> Union[Callable, None]:
    """Get unbound method for given bound method.

    Returns None if cannot get unbound method, or method is already unbound.
    """
    owner = getattr(func, "__self__", None)
    owner_class = type(owner)
    name = getattr(func, "__name__", None)
    instance_dict_overrides = getattr(owner, "__dict__", None)
    if (
        owner is not None
        and name
        and (
            not instance_dict_overrides
            or (instance_dict_overrides and name not in instance_dict_overrides)
        )
    ):
        return getattr(owner_class, name)
    return None


@undoc
@dataclass
class EvaluationPolicy:
    """Definition of evaluation policy."""

    allow_locals_access: bool = False
    allow_globals_access: bool = False
    allow_item_access: bool = False
    allow_attr_access: bool = False
    allow_builtins_access: bool = False
    allow_all_operations: bool = False
    allow_any_calls: bool = False
    allow_auto_import: bool = False
    allowed_calls: set[Callable] = field(default_factory=set)

    def can_get_item(self, value, item):
        return self.allow_item_access

    def can_get_attr(self, value, attr):
        return self.allow_attr_access

    def can_operate(self, dunders: tuple[str, ...], a, b=None):
        if self.allow_all_operations:
            return True

    def can_call(self, func):
        if self.allow_any_calls:
            return True

        if func in self.allowed_calls:
            return True

        owner_method = _unbind_method(func)

        if owner_method and owner_method in self.allowed_calls:
            return True


def _get_external(module_name: str, access_path: Sequence[str]):
    """Get value from external module given a dotted access path.

    Only gets value if the module is already imported.

    Raises:
    * `KeyError` if module is removed not found, and
    * `AttributeError` if access path does not match an exported object
    """
    try:
        member_type = sys.modules[module_name]
        # standard module
        for attr in access_path:
            member_type = getattr(member_type, attr)
        return member_type
    except (KeyError, AttributeError):
        # handle modules in namespace packages
        module_path = ".".join([module_name, *access_path])
        if module_path in sys.modules:
            return sys.modules[module_path]
        raise


def _has_original_dunder_external(
    value,
    module_name: str,
    access_path: Sequence[str],
    method_name: str,
):
    if module_name not in sys.modules:
        full_module_path = ".".join([module_name, *access_path])
        if full_module_path not in sys.modules:
            # LBYLB as it is faster
            return False
    try:
        member_type = _get_external(module_name, access_path)
        value_type = type(value)
        if type(value) == member_type:
            return True
        if isinstance(member_type, ModuleType):
            value_module = getmodule(value_type)
            if not value_module or not value_module.__name__:
                return False
            if (
                value_module.__name__ == member_type.__name__
                or value_module.__name__.startswith(member_type.__name__ + ".")
            ):
                return True
        if method_name == "__getattribute__":
            # we have to short-circuit here due to an unresolved issue in
            # `isinstance` implementation: https://bugs.python.org/issue32683
            return False
        if not isinstance(member_type, ModuleType) and isinstance(value, member_type):
            method = getattr(value_type, method_name, None)
            member_method = getattr(member_type, method_name, None)
            if member_method == method:
                return True
        if isinstance(member_type, ModuleType):
            method = getattr(value_type, method_name, None)
            for base_class in value_type.__mro__[1:]:
                base_module = getmodule(base_class)
                if base_module and (
                    base_module.__name__ == member_type.__name__
                    or base_module.__name__.startswith(member_type.__name__ + ".")
                ):
                    # Check if the method comes from this trusted base class
                    base_method = getattr(base_class, method_name, None)
                    if base_method is not None and base_method == method:
                        return True
    except (AttributeError, KeyError):
        return False


def _has_original_dunder(
    value, allowed_types, allowed_methods, allowed_external, method_name
):
    # note: Python ignores `__getattr__`/`__getitem__` on instances,
    # we only need to check at class level
    value_type = type(value)

    # strict type check passes → no need to check method
    if value_type in allowed_types:
        return True

    method = getattr(value_type, method_name, None)

    if method is None:
        return None

    if method in allowed_methods:
        return True

    for module_name, *access_path in allowed_external:
        if _has_original_dunder_external(value, module_name, access_path, method_name):
            return True

    return False


def _coerce_path_to_tuples(
    allow_list: set[tuple[str, ...] | str],
) -> set[tuple[str, ...]]:
    """Replace dotted paths on the provided allow-list with tuples."""
    return {
        path if isinstance(path, tuple) else tuple(path.split("."))
        for path in allow_list
    }


@undoc
@dataclass
class SelectivePolicy(EvaluationPolicy):
    allowed_getitem: set[InstancesHaveGetItem] = field(default_factory=set)
    allowed_getitem_external: set[tuple[str, ...] | str] = field(default_factory=set)

    allowed_getattr: set[MayHaveGetattr] = field(default_factory=set)
    allowed_getattr_external: set[tuple[str, ...] | str] = field(default_factory=set)

    allowed_operations: set = field(default_factory=set)
    allowed_operations_external: set[tuple[str, ...] | str] = field(default_factory=set)

    allow_getitem_on_types: bool = field(default_factory=bool)

    _operation_methods_cache: dict[str, set[Callable]] = field(
        default_factory=dict, init=False
    )

    def can_get_attr(self, value, attr):
        allowed_getattr_external = _coerce_path_to_tuples(self.allowed_getattr_external)

        has_original_attribute = _has_original_dunder(
            value,
            allowed_types=self.allowed_getattr,
            allowed_methods=self._getattribute_methods,
            allowed_external=allowed_getattr_external,
            method_name="__getattribute__",
        )
        has_original_attr = _has_original_dunder(
            value,
            allowed_types=self.allowed_getattr,
            allowed_methods=self._getattr_methods,
            allowed_external=allowed_getattr_external,
            method_name="__getattr__",
        )

        accept = False

        # Many objects do not have `__getattr__`, this is fine.
        if has_original_attr is None and has_original_attribute:
            accept = True
        else:
            # Accept objects without modifications to `__getattr__` and `__getattribute__`
            accept = has_original_attr and has_original_attribute

        if accept:
            # We still need to check for overridden properties.

            value_class = type(value)
            if not hasattr(value_class, attr):
                return True

            class_attr_val = getattr(value_class, attr)
            is_property = isinstance(class_attr_val, property)

            if not is_property:
                return True

            # Properties in allowed types are ok (although we do not include any
            # properties in our default allow list currently).
            if type(value) in self.allowed_getattr:
                return True  # pragma: no cover

            # Properties in subclasses of allowed types may be ok if not changed
            for module_name, *access_path in allowed_getattr_external:
                try:
                    external_class = _get_external(module_name, access_path)
                    external_class_attr_val = getattr(external_class, attr)
                except (KeyError, AttributeError):
                    return False  # pragma: no cover
                return class_attr_val == external_class_attr_val

        return False

    def can_get_item(self, value, item):
        """Allow accessing `__getiitem__` of allow-listed instances unless it was not modified."""
        allowed_getitem_external = _coerce_path_to_tuples(self.allowed_getitem_external)
        if self.allow_getitem_on_types:
            # e.g. Union[str, int] or Literal[True, 1]
            if isinstance(value, (typing._SpecialForm, typing._BaseGenericAlias)):
                return True
            # PEP 560 e.g. list[str]
            if isinstance(value, type) and hasattr(value, "__class_getitem__"):
                return True
        return _has_original_dunder(
            value,
            allowed_types=self.allowed_getitem,
            allowed_methods=self._getitem_methods,
            allowed_external=allowed_getitem_external,
            method_name="__getitem__",
        )

    def can_operate(self, dunders: tuple[str, ...], a, b=None):
        allowed_operations_external = _coerce_path_to_tuples(
            self.allowed_operations_external
        )
        objects = [a]
        if b is not None:
            objects.append(b)
        return all(
            [
                _has_original_dunder(
                    obj,
                    allowed_types=self.allowed_operations,
                    allowed_methods=self._operator_dunder_methods(dunder),
                    allowed_external=allowed_operations_external,
                    method_name=dunder,
                )
                for dunder in dunders
                for obj in objects
            ]
        )

    def _operator_dunder_methods(self, dunder: str) -> set[Callable]:
        if dunder not in self._operation_methods_cache:
            self._operation_methods_cache[dunder] = self._safe_get_methods(
                self.allowed_operations, dunder
            )
        return self._operation_methods_cache[dunder]

    @cached_property
    def _getitem_methods(self) -> set[Callable]:
        return self._safe_get_methods(self.allowed_getitem, "__getitem__")

    @cached_property
    def _getattr_methods(self) -> set[Callable]:
        return self._safe_get_methods(self.allowed_getattr, "__getattr__")

    @cached_property
    def _getattribute_methods(self) -> set[Callable]:
        return self._safe_get_methods(self.allowed_getattr, "__getattribute__")

    def _safe_get_methods(self, classes, name) -> set[Callable]:
        return {
            method
            for class_ in classes
            for method in [getattr(class_, name, None)]
            if method
        }


class _DummyNamedTuple(NamedTuple):
    """Used internally to retrieve methods of named tuple instance."""


EvaluationPolicyName = Literal["forbidden", "minimal", "limited", "unsafe", "dangerous"]


@dataclass
class EvaluationContext:
    #: Local namespace
    locals: dict
    #: Global namespace
    globals: dict
    #: Evaluation policy identifier
    evaluation: EvaluationPolicyName = "forbidden"
    #: Whether the evaluation of code takes place inside of a subscript.
    #: Useful for evaluating ``:-1, 'col'`` in ``df[:-1, 'col']``.
    in_subscript: bool = False
    #: Auto import method
    auto_import: Callable[[Sequence[str]], ModuleType] | None = None
    #: Overrides for evaluation policy
    policy_overrides: dict = field(default_factory=dict)
    #: Transient local namespace used to store mocks
    transient_locals: dict = field(default_factory=dict)
    #: Transients of class level
    class_transients: dict | None = None
    #: Instance variable name used in the method definition
    instance_arg_name: str | None = None
    #: Currently associated value
    #: Useful for adding items to _Duck on annotated assignment
    current_value: ast.AST | None = None

    def replace(self, /, **changes):
        """Return a new copy of the context, with specified changes"""
        return dataclasses.replace(self, **changes)


class _IdentitySubscript:
    """Returns the key itself when item is requested via subscript."""

    def __getitem__(self, key):
        return key


IDENTITY_SUBSCRIPT = _IdentitySubscript()
SUBSCRIPT_MARKER = "__SUBSCRIPT_SENTINEL__"
UNKNOWN_SIGNATURE = Signature()
NOT_EVALUATED = object()


class GuardRejection(Exception):
    """Exception raised when guard rejects evaluation attempt."""

    pass


def guarded_eval(code: str, context: EvaluationContext):
    """Evaluate provided code in the evaluation context.

    If evaluation policy given by context is set to ``forbidden``
    no evaluation will be performed; if it is set to ``dangerous``
    standard :func:`eval` will be used; finally, for any other,
    policy :func:`eval_node` will be called on parsed AST.
    """
    locals_ = context.locals

    if context.evaluation == "forbidden":
        raise GuardRejection("Forbidden mode")

    # note: not using `ast.literal_eval` as it does not implement
    # getitem at all, for example it fails on simple `[0][1]`

    if context.in_subscript:
        # syntactic sugar for ellipsis (:) is only available in subscripts
        # so we need to trick the ast parser into thinking that we have
        # a subscript, but we need to be able to later recognise that we did
        # it so we can ignore the actual __getitem__ operation
        if not code:
            return tuple()
        locals_ = locals_.copy()
        locals_[SUBSCRIPT_MARKER] = IDENTITY_SUBSCRIPT
        code = SUBSCRIPT_MARKER + "[" + code + "]"
        context = context.replace(locals=locals_)

    if context.evaluation == "dangerous":
        return eval(code, context.globals, context.locals)

    node = ast.parse(code, mode="exec")

    return eval_node(node, context)


BINARY_OP_DUNDERS: dict[type[ast.operator], tuple[str]] = {
    ast.Add: ("__add__",),
    ast.Sub: ("__sub__",),
    ast.Mult: ("__mul__",),
    ast.Div: ("__truediv__",),
    ast.FloorDiv: ("__floordiv__",),
    ast.Mod: ("__mod__",),
    ast.Pow: ("__pow__",),
    ast.LShift: ("__lshift__",),
    ast.RShift: ("__rshift__",),
    ast.BitOr: ("__or__",),
    ast.BitXor: ("__xor__",),
    ast.BitAnd: ("__and__",),
    ast.MatMult: ("__matmul__",),
}

COMP_OP_DUNDERS: dict[type[ast.cmpop], tuple[str, ...]] = {
    ast.Eq: ("__eq__",),
    ast.NotEq: ("__ne__", "__eq__"),
    ast.Lt: ("__lt__", "__gt__"),
    ast.LtE: ("__le__", "__ge__"),
    ast.Gt: ("__gt__", "__lt__"),
    ast.GtE: ("__ge__", "__le__"),
    ast.In: ("__contains__",),
    # Note: ast.Is, ast.IsNot, ast.NotIn are handled specially
}

UNARY_OP_DUNDERS: dict[type[ast.unaryop], tuple[str, ...]] = {
    ast.USub: ("__neg__",),
    ast.UAdd: ("__pos__",),
    # we have to check both __inv__ and __invert__!
    ast.Invert: ("__invert__", "__inv__"),
    ast.Not: ("__not__",),
}

GENERIC_CONTAINER_TYPES = (dict, list, set, tuple, frozenset)


class ImpersonatingDuck:
    """A dummy class used to create objects of other classes without calling their ``__init__``"""

    # no-op: override __class__ to impersonate


class _Duck:
    """A dummy class used to create objects pretending to have given attributes"""

    def __init__(self, attributes: Optional[dict] = None, items: Optional[dict] = None):
        self.attributes = attributes if attributes is not None else {}
        self.items = items if items is not None else {}

    def __getattr__(self, attr: str):
        return self.attributes[attr]

    def __hasattr__(self, attr: str):
        return attr in self.attributes

    def __dir__(self):
        return [*dir(super), *self.attributes]

    def __getitem__(self, key: str):
        return self.items[key]

    def __hasitem__(self, key: str):
        return self.items[key]

    def _ipython_key_completions_(self):
        return self.items.keys()


def _find_dunder(node_op, dunders) -> Union[tuple[str, ...], None]:
    dunder = None
    for op, candidate_dunder in dunders.items():
        if isinstance(node_op, op):
            dunder = candidate_dunder
    return dunder


def get_policy(context: EvaluationContext) -> EvaluationPolicy:
    policy = copy(EVALUATION_POLICIES[context.evaluation])

    for key, value in context.policy_overrides.items():
        if hasattr(policy, key):
            setattr(policy, key, value)
    return policy


def _validate_policy_overrides(
    policy_name: EvaluationPolicyName, policy_overrides: dict
) -> bool:
    policy = EVALUATION_POLICIES[policy_name]

    all_good = True
    for key, value in policy_overrides.items():
        if not hasattr(policy, key):
            warnings.warn(
                f"Override {key!r} is not valid with {policy_name!r} evaluation policy"
            )
            all_good = False
    return all_good


def _is_type_annotation(obj) -> bool:
    """
    Returns True if obj is a type annotation, False otherwise.
    """
    if isinstance(obj, type):
        return True
    if isinstance(obj, types.GenericAlias):
        return True
    if hasattr(types, "UnionType") and isinstance(obj, types.UnionType):
        return True
    if isinstance(obj, (typing._SpecialForm, typing._BaseGenericAlias)):
        return True
    if isinstance(obj, typing.TypeVar):
        return True
    # Types that support __class_getitem__
    if isinstance(obj, type) and hasattr(obj, "__class_getitem__"):
        return True
    # Fallback: check if get_origin returns something
    if hasattr(typing, "get_origin") and get_origin(obj) is not None:
        return True

    return False


def _handle_assign(node: ast.Assign, context: EvaluationContext):
    value = eval_node(node.value, context)
    transient_locals = context.transient_locals
    policy = get_policy(context)
    class_transients = context.class_transients
    for target in node.targets:
        if isinstance(target, (ast.Tuple, ast.List)):
            # Handle unpacking assignment
            values = list(value)
            targets = target.elts
            starred = [i for i, t in enumerate(targets) if isinstance(t, ast.Starred)]

            # Unified handling: treat no starred as starred at end
            star_or_last_idx = starred[0] if starred else len(targets)

            # Before starred
            for i in range(star_or_last_idx):
                # Check for self.x assignment
                if _is_instance_attribute_assignment(targets[i], context):
                    class_transients[targets[i].attr] = values[i]
                else:
                    transient_locals[targets[i].id] = values[i]

            # Starred if exists
            if starred:
                end = len(values) - (len(targets) - star_or_last_idx - 1)
                if _is_instance_attribute_assignment(
                    targets[star_or_last_idx], context
                ):
                    class_transients[targets[star_or_last_idx].attr] = values[
                        star_or_last_idx:end
                    ]
                else:
                    transient_locals[targets[star_or_last_idx].value.id] = values[
                        star_or_last_idx:end
                    ]

                # After starred
                for i in range(star_or_last_idx + 1, len(targets)):
                    if _is_instance_attribute_assignment(targets[i], context):
                        class_transients[targets[i].attr] = values[
                            len(values) - (len(targets) - i)
                        ]
                    else:
                        transient_locals[targets[i].id] = values[
                            len(values) - (len(targets) - i)
                        ]
        elif isinstance(target, ast.Subscript):
            if isinstance(target.value, ast.Name):
                name = target.value.id
                container = transient_locals.get(name)
                if container is None:
                    container = context.locals.get(name)
                if container is None:
                    container = context.globals.get(name)
                if container is None:
                    raise NameError(
                        f"{name} not found in locals, globals, nor builtins"
                    )
                storage_dict = transient_locals
                storage_key = name
            elif isinstance(
                target.value, ast.Attribute
            ) and _is_instance_attribute_assignment(target.value, context):
                attr = target.value.attr
                container = class_transients.get(attr, None)
                if container is None:
                    raise NameError(f"{attr} not found in class transients")
                storage_dict = class_transients
                storage_key = attr
            else:
                return

            key = eval_node(target.slice, context)
            attributes = (
                dict.fromkeys(dir(container))
                if policy.can_call(container.__dir__)
                else {}
            )
            items = {}

            if policy.can_get_item(container, None):
                try:
                    items = dict(container.items())
                except Exception:
                    pass

            items[key] = value
            duck_container = _Duck(attributes=attributes, items=items)
            storage_dict[storage_key] = duck_container
        elif _is_instance_attribute_assignment(target, context):
            class_transients[target.attr] = value
        else:
            transient_locals[target.id] = value
    return None


def _handle_annassign(node, context):
    context_with_value = context.replace(current_value=getattr(node, "value", None))
    annotation_result = eval_node(node.annotation, context_with_value)
    if _is_type_annotation(annotation_result):
        annotation_value = _resolve_annotation(annotation_result, context)
        # Use Value for generic types
        use_value = (
            isinstance(annotation_value, GENERIC_CONTAINER_TYPES) and node.value is not None
        )
    else:
        annotation_value = annotation_result
        use_value = False

    # LOCAL VARIABLE
    if getattr(node, "simple", False) and isinstance(node.target, ast.Name):
        name = node.target.id
        if use_value:
            return _handle_assign(
                ast.Assign(targets=[node.target], value=node.value), context
            )
        context.transient_locals[name] = annotation_value
        return None

    # INSTANCE ATTRIBUTE
    if _is_instance_attribute_assignment(node.target, context):
        attr = node.target.attr
        if use_value:
            return _handle_assign(
                ast.Assign(targets=[node.target], value=node.value), context
            )
        context.class_transients[attr] = annotation_value
        return None

    return None

def _extract_args_and_kwargs(node: ast.Call, context: EvaluationContext):
    args = [eval_node(arg, context) for arg in node.args]
    kwargs = {
        k: v
        for kw in node.keywords
        for k, v in (
            {kw.arg: eval_node(kw.value, context)}
            if kw.arg
            else eval_node(kw.value, context)
        ).items()
    }
    return args, kwargs


def _is_instance_attribute_assignment(
    target: ast.AST, context: EvaluationContext
) -> bool:
    """Return True if target is an attribute access on the instance argument."""
    return (
        context.class_transients is not None
        and context.instance_arg_name is not None
        and isinstance(target, ast.Attribute)
        and isinstance(getattr(target, "value", None), ast.Name)
        and getattr(target.value, "id", None) == context.instance_arg_name
    )


def _get_coroutine_attributes() -> dict[str, Optional[object]]:
    async def _dummy():
        return None

    coro = _dummy()
    try:
        return {attr: getattr(coro, attr, None) for attr in dir(coro)}
    finally:
        coro.close()


def eval_node(node: Union[ast.AST, None], context: EvaluationContext):
    """Evaluate AST node in provided context.

    Applies evaluation restrictions defined in the context. Currently does not support evaluation of functions with keyword arguments.

    Does not evaluate actions that always have side effects:

    - class definitions (``class sth: ...``)
    - function definitions (``def sth: ...``)
    - variable assignments (``x = 1``)
    - augmented assignments (``x += 1``)
    - deletions (``del x``)

    Does not evaluate operations which do not return values:

    - assertions (``assert x``)
    - pass (``pass``)
    - imports (``import x``)
    - control flow:

        - conditionals (``if x:``) except for ternary IfExp (``a if x else b``)
        - loops (``for`` and ``while``)
        - exception handling

    The purpose of this function is to guard against unwanted side-effects;
    it does not give guarantees on protection from malicious code execution.
    """
    policy = get_policy(context)

    if node is None:
        return None
    if isinstance(node, (ast.Interactive, ast.Module)):
        result = None
        for child_node in node.body:
            result = eval_node(child_node, context)
        return result
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        is_async = isinstance(node, ast.AsyncFunctionDef)
        func_locals = context.transient_locals.copy()
        func_context = context.replace(transient_locals=func_locals)
        is_property = False
        is_static = False
        is_classmethod = False
        for decorator_node in node.decorator_list:
            try:
                decorator = eval_node(decorator_node, context)
            except NameError:
                # if the decorator is not yet defined this is fine
                # especialy because we don't handle imports yet
                continue
            if decorator is property:
                is_property = True
            elif decorator is staticmethod:
                is_static = True
            elif decorator is classmethod:
                is_classmethod = True

        if func_context.class_transients is not None:
            if not is_static and not is_classmethod:
                func_context.instance_arg_name = (
                    node.args.args[0].arg if node.args.args else None
                )

        return_type = eval_node(node.returns, context=context)

        for child_node in node.body:
            eval_node(child_node, func_context)

        if is_property:
            if return_type is not None:
                if _is_type_annotation(return_type):
                    context.transient_locals[node.name] = _resolve_annotation(
                        return_type, context
                    )
                else:
                    context.transient_locals[node.name] = return_type
            else:
                return_value = _infer_return_value(node, func_context)
                context.transient_locals[node.name] = return_value

            return None

        def dummy_function(*args, **kwargs):
            pass

        if return_type is not None:
            if _is_type_annotation(return_type):
                dummy_function.__annotations__["return"] = return_type
            else:
                dummy_function.__inferred_return__ = return_type
        else:
            inferred_return = _infer_return_value(node, func_context)
            if inferred_return is not None:
                dummy_function.__inferred_return__ = inferred_return

        dummy_function.__name__ = node.name
        dummy_function.__node__ = node
        dummy_function.__is_async__ = is_async
        context.transient_locals[node.name] = dummy_function
        return None
    if isinstance(node, ast.Lambda):

        def dummy_function(*args, **kwargs):
            pass

        dummy_function.__inferred_return__ = eval_node(node.body, context)
        return dummy_function
    if isinstance(node, ast.ClassDef):
        # TODO support class decorators?
        class_locals = {}
        outer_locals = context.locals.copy()
        outer_locals.update(context.transient_locals)
        class_context = context.replace(
            transient_locals=class_locals, locals=outer_locals
        )
        class_context.class_transients = class_locals
        for child_node in node.body:
            eval_node(child_node, class_context)
        bases = tuple([eval_node(base, context) for base in node.bases])
        dummy_class = type(node.name, bases, class_locals)
        context.transient_locals[node.name] = dummy_class
        return None
    if isinstance(node, ast.Await):
        value = eval_node(node.value, context)
        if hasattr(value, "__awaited_type__"):
            return value.__awaited_type__
        return value
    if isinstance(node, ast.While):
        loop_locals = context.transient_locals.copy()
        loop_context = context.replace(transient_locals=loop_locals)

        result = None
        for stmt in node.body:
            result = eval_node(stmt, loop_context)

        policy = get_policy(context)
        merged_locals = _merge_dicts_by_key(
            [loop_locals, context.transient_locals.copy()], policy
        )
        context.transient_locals.update(merged_locals)

        return result
    if isinstance(node, ast.For):
        try:
            iterable = eval_node(node.iter, context)
        except Exception:
            iterable = None

        sample = None
        if iterable is not None:
            try:
                if policy.can_call(getattr(iterable, "__iter__", None)):
                    sample = next(iter(iterable))
            except Exception:
                sample = None

        loop_locals = context.transient_locals.copy()
        loop_context = context.replace(transient_locals=loop_locals)

        if sample is not None:
            try:
                fake_assign = ast.Assign(
                    targets=[node.target], value=ast.Constant(value=sample)
                )
                _handle_assign(fake_assign, loop_context)
            except Exception:
                pass

        result = None
        for stmt in node.body:
            result = eval_node(stmt, loop_context)

        policy = get_policy(context)
        merged_locals = _merge_dicts_by_key(
            [loop_locals, context.transient_locals.copy()], policy
        )
        context.transient_locals.update(merged_locals)

        return result
    if isinstance(node, ast.If):
        branches = []
        current = node
        result = None
        while True:
            branch_locals = context.transient_locals.copy()
            branch_context = context.replace(transient_locals=branch_locals)
            for stmt in current.body:
                result = eval_node(stmt, branch_context)
            branches.append(branch_locals)
            if not current.orelse:
                break
            elif len(current.orelse) == 1 and isinstance(current.orelse[0], ast.If):
                # It's an elif - continue loop
                current = current.orelse[0]
            else:
                # It's an else block - process and break
                else_locals = context.transient_locals.copy()
                else_context = context.replace(transient_locals=else_locals)
                for stmt in current.orelse:
                    result = eval_node(stmt, else_context)
                branches.append(else_locals)
                break
        branches.append(context.transient_locals.copy())
        policy = get_policy(context)
        merged_locals = _merge_dicts_by_key(branches, policy)
        context.transient_locals.update(merged_locals)
        return result
    if isinstance(node, ast.Assign):
        return _handle_assign(node, context)
    if isinstance(node, ast.AnnAssign):
        return _handle_annassign(node, context)
    if isinstance(node, ast.Expression):
        return eval_node(node.body, context)
    if isinstance(node, ast.Expr):
        return eval_node(node.value, context)
    if isinstance(node, ast.Pass):
        return None
    if isinstance(node, ast.Import):
        # TODO: populate transient_locals
        return None
    if isinstance(node, (ast.AugAssign, ast.Delete)):
        return None
    if isinstance(node, (ast.Global, ast.Nonlocal)):
        return None
    if isinstance(node, ast.BinOp):
        left = eval_node(node.left, context)
        right = eval_node(node.right, context)
        if (
            isinstance(node.op, ast.BitOr)
            and _is_type_annotation(left)
            and _is_type_annotation(right)
        ):
            left_duck = (
                _Duck(dict.fromkeys(dir(left)))
                if policy.can_call(left.__dir__)
                else _Duck()
            )
            right_duck = (
                _Duck(dict.fromkeys(dir(right)))
                if policy.can_call(right.__dir__)
                else _Duck()
            )
            value_node = context.current_value
            if value_node is not None and isinstance(value_node, ast.Dict):
                if dict in [left, right]:
                    return _merge_values(
                        [left_duck, right_duck, ast.literal_eval(value_node)],
                        policy=get_policy(context),
                    )
            return _merge_values([left_duck, right_duck], policy=get_policy(context))
        dunders = _find_dunder(node.op, BINARY_OP_DUNDERS)
        if dunders:
            if policy.can_operate(dunders, left, right):
                return getattr(left, dunders[0])(right)
            else:
                raise GuardRejection(
                    f"Operation (`{dunders}`) for",
                    type(left),
                    f"not allowed in {context.evaluation} mode",
                )
    if isinstance(node, ast.Compare):
        left = eval_node(node.left, context)
        all_true = True
        negate = False
        for op, right in zip(node.ops, node.comparators):
            right = eval_node(right, context)
            dunder = None
            dunders = _find_dunder(op, COMP_OP_DUNDERS)
            if not dunders:
                if isinstance(op, ast.NotIn):
                    dunders = COMP_OP_DUNDERS[ast.In]
                    negate = True
                if isinstance(op, ast.Is):
                    dunder = "is_"
                if isinstance(op, ast.IsNot):
                    dunder = "is_"
                    negate = True
            if not dunder and dunders:
                dunder = dunders[0]
            if dunder:
                a, b = (right, left) if dunder == "__contains__" else (left, right)
                if dunder == "is_" or dunders and policy.can_operate(dunders, a, b):
                    result = getattr(operator, dunder)(a, b)
                    if negate:
                        result = not result
                    if not result:
                        all_true = False
                    left = right
                else:
                    raise GuardRejection(
                        f"Comparison (`{dunder}`) for",
                        type(left),
                        f"not allowed in {context.evaluation} mode",
                    )
            else:
                raise ValueError(
                    f"Comparison `{dunder}` not supported"
                )  # pragma: no cover
        return all_true
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Tuple):
        return tuple(eval_node(e, context) for e in node.elts)
    if isinstance(node, ast.List):
        return [eval_node(e, context) for e in node.elts]
    if isinstance(node, ast.Set):
        return {eval_node(e, context) for e in node.elts}
    if isinstance(node, ast.Dict):
        return dict(
            zip(
                [eval_node(k, context) for k in node.keys],
                [eval_node(v, context) for v in node.values],
            )
        )
    if isinstance(node, ast.Slice):
        return slice(
            eval_node(node.lower, context),
            eval_node(node.upper, context),
            eval_node(node.step, context),
        )
    if isinstance(node, ast.UnaryOp):
        value = eval_node(node.operand, context)
        dunders = _find_dunder(node.op, UNARY_OP_DUNDERS)
        if dunders:
            if policy.can_operate(dunders, value):
                try:
                    return getattr(value, dunders[0])()
                except AttributeError:
                    raise TypeError(
                        f"bad operand type for unary {node.op}: {type(value)}"
                    )
            else:
                raise GuardRejection(
                    f"Operation (`{dunders}`) for",
                    type(value),
                    f"not allowed in {context.evaluation} mode",
                )
    if isinstance(node, ast.Subscript):
        value = eval_node(node.value, context)
        slice_ = eval_node(node.slice, context)
        if policy.can_get_item(value, slice_):
            return value[slice_]
        raise GuardRejection(
            "Subscript access (`__getitem__`) for",
            type(value),  # not joined to avoid calling `repr`
            f" not allowed in {context.evaluation} mode",
        )
    if isinstance(node, ast.Name):
        return _eval_node_name(node.id, context)
    if isinstance(node, ast.Attribute):
        if (
            context.class_transients is not None
            and isinstance(node.value, ast.Name)
            and node.value.id == context.instance_arg_name
        ):
            return context.class_transients.get(node.attr)
        value = eval_node(node.value, context)
        if policy.can_get_attr(value, node.attr):
            return getattr(value, node.attr)
        try:
            cls = (
                value if isinstance(value, type) else getattr(value, "__class__", None)
            )
            if cls is not None:
                resolved_hints = get_type_hints(
                    cls,
                    globalns=(context.globals or {}),
                    localns=(context.locals or {}),
                )
                if node.attr in resolved_hints:
                    annotated = resolved_hints[node.attr]
                    return _resolve_annotation(annotated, context)
        except Exception:
            # Fall through to the guard rejection
            pass
        raise GuardRejection(
            "Attribute access (`__getattr__`) for",
            type(value),  # not joined to avoid calling `repr`
            f"not allowed in {context.evaluation} mode",
        )
    if isinstance(node, ast.IfExp):
        test = eval_node(node.test, context)
        if test:
            return eval_node(node.body, context)
        else:
            return eval_node(node.orelse, context)
    if isinstance(node, ast.Call):
        func = eval_node(node.func, context)
        if policy.can_call(func):
            args, kwargs = _extract_args_and_kwargs(node, context)
            return func(*args, **kwargs)
        if isclass(func):
            # this code path gets entered when calling class e.g. `MyClass()`
            # or `my_instance.__class__()` - in both cases `func` is `MyClass`.
            # Should return `MyClass` if `__new__` is not overridden,
            # otherwise whatever `__new__` return type is.
            overridden_return_type = _eval_return_type(func.__new__, node, context)
            if overridden_return_type is not NOT_EVALUATED:
                return overridden_return_type
            return _create_duck_for_heap_type(func)
        else:
            inferred_return = getattr(func, "__inferred_return__", NOT_EVALUATED)
            return_type = _eval_return_type(func, node, context)
            if getattr(func, "__is_async__", False):
                awaited_type = (
                    inferred_return if inferred_return is not None else return_type
                )
                coroutine_duck = _Duck(attributes=_get_coroutine_attributes())
                coroutine_duck.__awaited_type__ = awaited_type
                return coroutine_duck
            if inferred_return is not NOT_EVALUATED:
                return inferred_return
            if return_type is not NOT_EVALUATED:
                return return_type
        raise GuardRejection(
            "Call for",
            func,  # not joined to avoid calling `repr`
            f"not allowed in {context.evaluation} mode",
        )
    if isinstance(node, ast.Assert):
        # message is always the second item, so if it is defined user would be completing
        # on the message, not on the assertion test
        if node.msg:
            return eval_node(node.msg, context)
        return eval_node(node.test, context)
    return None


def _merge_dicts_by_key(dicts: list, policy: EvaluationPolicy):
    """Merge multiple dictionaries, combining values for each key."""
    if len(dicts) == 1:
        return dicts[0]

    all_keys = set()
    for d in dicts:
        all_keys.update(d.keys())

    merged = {}
    for key in all_keys:
        values = [d[key] for d in dicts if key in d]
        if values:
            merged[key] = _merge_values(values, policy)

    return merged


def _merge_values(values, policy: EvaluationPolicy):
    """Recursively merge multiple values, combining attributes and dict items."""
    if len(values) == 1:
        return values[0]

    types = {type(v) for v in values}
    merged_items = None
    key_values = {}
    attributes = set()
    for v in values:
        if policy.can_call(v.__dir__):
            attributes.update(dir(v))
        try:
            if policy.can_call(v.items):
                try:
                    for k, val in v.items():
                        key_values.setdefault(k, []).append(val)
                except Exception as e:
                    pass
            elif policy.can_call(v.keys):
                try:
                    for k in v.keys():
                        key_values.setdefault(k, []).append(None)
                except Exception as e:
                    pass
        except Exception as e:
            pass

    if key_values:
        merged_items = {
            k: _merge_values(vals, policy) if vals[0] is not None else None
            for k, vals in key_values.items()
        }

    if len(types) == 1:
        t = next(iter(types))
        if t not in (dict,) and not (
            hasattr(next(iter(values)), "__getitem__")
            and (
                hasattr(next(iter(values)), "items")
                or hasattr(next(iter(values)), "keys")
            )
        ):
            if t in (list, set, tuple):
                return t
            return values[0]

    return _Duck(attributes=dict.fromkeys(attributes), items=merged_items)


def _infer_return_value(node: ast.FunctionDef, context: EvaluationContext):
    """Infer the return value(s) of a function by evaluating all return statements."""
    return_values = _collect_return_values(node.body, context)

    if not return_values:
        return None
    if len(return_values) == 1:
        return return_values[0]

    policy = get_policy(context)
    return _merge_values(return_values, policy)


def _collect_return_values(body, context):
    """Recursively collect return values from a list of AST statements."""
    return_values = []
    for stmt in body:
        if isinstance(stmt, ast.Return):
            if stmt.value is None:
                continue
            try:
                value = eval_node(stmt.value, context)
                if value is not None and value is not NOT_EVALUATED:
                    return_values.append(value)
            except Exception:
                pass
        if isinstance(
            stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)
        ):
            continue
        elif hasattr(stmt, "body") and isinstance(stmt.body, list):
            return_values.extend(_collect_return_values(stmt.body, context))
        if isinstance(stmt, ast.Try):
            for h in stmt.handlers:
                if hasattr(h, "body"):
                    return_values.extend(_collect_return_values(h.body, context))
            if hasattr(stmt, "orelse"):
                return_values.extend(_collect_return_values(stmt.orelse, context))
            if hasattr(stmt, "finalbody"):
                return_values.extend(_collect_return_values(stmt.finalbody, context))
        if hasattr(stmt, "orelse") and isinstance(stmt.orelse, list):
            return_values.extend(_collect_return_values(stmt.orelse, context))
    return return_values


def _eval_return_type(func: Callable, node: ast.Call, context: EvaluationContext):
    """Evaluate return type of a given callable function.

    Returns the built-in type, a duck or NOT_EVALUATED sentinel.
    """
    try:
        sig = signature(func)
    except ValueError:
        sig = UNKNOWN_SIGNATURE
    # if annotation was not stringized, or it was stringized
    # but resolved by signature call we know the return type
    not_empty = sig.return_annotation is not Signature.empty
    if not_empty:
        return _resolve_annotation(sig.return_annotation, context, sig, func, node)
    return NOT_EVALUATED


def _eval_annotation(
    annotation: str,
    context: EvaluationContext,
):
    return (
        _eval_node_name(annotation, context)
        if isinstance(annotation, str)
        else annotation
    )


class _GetItemDuck(dict):
    """A dict subclass that always returns the factory instance and claims to have any item."""

    def __init__(self, factory, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._factory = factory

    def __getitem__(self, key):
        return self._factory()

    def __contains__(self, key):
        return True


def _resolve_annotation(
    annotation: object | str,
    context: EvaluationContext,
    sig: Signature | None = None,
    func: Callable | None = None,
    node: ast.Call | None = None,
):
    """Resolve annotation created by user with `typing` module and custom objects."""
    if annotation is None:
        return None
    annotation = _eval_annotation(annotation, context)
    origin = get_origin(annotation)
    if annotation is Self and func and hasattr(func, "__self__"):
        return func.__self__
    elif origin is Literal:
        type_args = get_args(annotation)
        if len(type_args) == 1:
            return type_args[0]
    elif annotation is LiteralString:
        return ""
    elif annotation is AnyStr:
        index = None
        if func and hasattr(func, "__node__"):
            def_node = func.__node__
            for i, arg in enumerate(def_node.args.args):
                if not arg.annotation:
                    continue
                annotation = _eval_annotation(arg.annotation.id, context)
                if annotation is AnyStr:
                    index = i
                    break
            is_bound_method = (
                isinstance(func, MethodType) and getattr(func, "__self__") is not None
            )
            if index and is_bound_method:
                index -= 1
        elif sig:
            for i, (key, value) in enumerate(sig.parameters.items()):
                if value.annotation is AnyStr:
                    index = i
                    break
        if index is None:
            return None
        if index < 0 or index >= len(node.args):
            return None
        return eval_node(node.args[index], context)
    elif origin is TypeGuard:
        return False
    elif origin is set or origin is list:
        # only one type argument allowed
        attributes = [
            attr
            for attr in dir(
                _resolve_annotation(get_args(annotation)[0], context, sig, func, node)
            )
        ]
        duck = _Duck(attributes=dict.fromkeys(attributes))
        return _Duck(
            attributes=dict.fromkeys(dir(origin())),
            # items are not strrictly needed for set
            items=_GetItemDuck(lambda: duck),
        )
    elif origin is tuple:
        # multiple type arguments
        return tuple(
            _resolve_annotation(arg, context, sig, func, node)
            for arg in get_args(annotation)
        )
    elif origin is Union:
        # multiple type arguments
        attributes = [
            attr
            for type_arg in get_args(annotation)
            for attr in dir(_resolve_annotation(type_arg, context, sig, func, node))
        ]
        return _Duck(attributes=dict.fromkeys(attributes))
    elif is_typeddict(annotation):
        return _Duck(
            attributes=dict.fromkeys(dir(dict())),
            items={
                k: _resolve_annotation(v, context, sig, func, node)
                for k, v in annotation.__annotations__.items()
            },
        )
    elif hasattr(annotation, "_is_protocol"):
        return _Duck(attributes=dict.fromkeys(dir(annotation)))
    elif origin is Annotated:
        type_arg = get_args(annotation)[0]
        return _resolve_annotation(type_arg, context, sig, func, node)
    elif isinstance(annotation, NewType):
        return _eval_or_create_duck(annotation.__supertype__, context)
    elif isinstance(annotation, TypeAliasType):
        return _eval_or_create_duck(annotation.__value__, context)
    else:
        return _eval_or_create_duck(annotation, context)


def _eval_node_name(node_id: str, context: EvaluationContext):
    policy = get_policy(context)
    if node_id in context.transient_locals:
        return context.transient_locals[node_id]
    if policy.allow_locals_access and node_id in context.locals:
        return context.locals[node_id]
    if policy.allow_globals_access and node_id in context.globals:
        return context.globals[node_id]
    if policy.allow_builtins_access and hasattr(builtins, node_id):
        # note: do not use __builtins__, it is implementation detail of cPython
        return getattr(builtins, node_id)
    if policy.allow_auto_import and context.auto_import:
        return context.auto_import(node_id)
    if not policy.allow_globals_access and not policy.allow_locals_access:
        raise GuardRejection(
            f"Namespace access not allowed in {context.evaluation} mode"
        )
    else:
        raise NameError(f"{node_id} not found in locals, globals, nor builtins")


def _eval_or_create_duck(duck_type, context: EvaluationContext):
    policy = get_policy(context)
    # if allow-listed builtin is on type annotation, instantiate it
    if policy.can_call(duck_type):
        return duck_type()
    # if custom class is in type annotation, mock it
    return _create_duck_for_heap_type(duck_type)


def _create_duck_for_heap_type(duck_type):
    """Create an imitation of an object of a given type (a duck).

    Returns the duck or NOT_EVALUATED sentinel if duck could not be created.
    """
    duck = ImpersonatingDuck()
    try:
        # this only works for heap types, not builtins
        duck.__class__ = duck_type
        return duck
    except TypeError:
        pass
    return NOT_EVALUATED


SUPPORTED_EXTERNAL_GETITEM = {
    ("pandas", "core", "indexing", "_iLocIndexer"),
    ("pandas", "core", "indexing", "_LocIndexer"),
    ("pandas", "DataFrame"),
    ("pandas", "Series"),
    ("numpy", "ndarray"),
    ("numpy", "void"),
}


BUILTIN_GETITEM: set[InstancesHaveGetItem] = {
    dict,
    str,  # type: ignore[arg-type]
    bytes,  # type: ignore[arg-type]
    list,
    tuple,
    type,  # for type annotations like list[str]
    _Duck,
    collections.defaultdict,
    collections.deque,
    collections.OrderedDict,
    collections.ChainMap,
    collections.UserDict,
    collections.UserList,
    collections.UserString,  # type: ignore[arg-type]
    _DummyNamedTuple,
    _IdentitySubscript,
}


def _list_methods(cls, source=None):
    """For use on immutable objects or with methods returning a copy"""
    return [getattr(cls, k) for k in (source if source else dir(cls))]


dict_non_mutating_methods = ("copy", "keys", "values", "items")
list_non_mutating_methods = ("copy", "index", "count")
set_non_mutating_methods = set(dir(set)) & set(dir(frozenset))


dict_keys: type[collections.abc.KeysView] = type({}.keys())
dict_values: type = type({}.values())
dict_items: type = type({}.items())

NUMERICS = {int, float, complex}

ALLOWED_CALLS = {
    bytes,
    *_list_methods(bytes),
    bytes.__iter__,
    dict,
    *_list_methods(dict, dict_non_mutating_methods),
    dict.__iter__,
    dict_keys.__iter__,
    dict_values.__iter__,
    dict_items.__iter__,
    dict_keys.isdisjoint,
    list,
    *_list_methods(list, list_non_mutating_methods),
    list.__iter__,
    set,
    *_list_methods(set, set_non_mutating_methods),
    set.__iter__,
    frozenset,
    *_list_methods(frozenset),
    frozenset.__iter__,
    range,
    range.__iter__,
    str,
    *_list_methods(str),
    str.__iter__,
    tuple,
    *_list_methods(tuple),
    tuple.__iter__,
    bool,
    *_list_methods(bool),
    *NUMERICS,
    *[method for numeric_cls in NUMERICS for method in _list_methods(numeric_cls)],
    collections.deque,
    *_list_methods(collections.deque, list_non_mutating_methods),
    collections.deque.__iter__,
    collections.defaultdict,
    *_list_methods(collections.defaultdict, dict_non_mutating_methods),
    collections.defaultdict.__iter__,
    collections.OrderedDict,
    *_list_methods(collections.OrderedDict, dict_non_mutating_methods),
    collections.OrderedDict.__iter__,
    collections.UserDict,
    *_list_methods(collections.UserDict, dict_non_mutating_methods),
    collections.UserDict.__iter__,
    collections.UserList,
    *_list_methods(collections.UserList, list_non_mutating_methods),
    collections.UserList.__iter__,
    collections.UserString,
    *_list_methods(collections.UserString, dir(str)),
    collections.UserString.__iter__,
    collections.Counter,
    *_list_methods(collections.Counter, dict_non_mutating_methods),
    collections.Counter.__iter__,
    collections.Counter.elements,
    collections.Counter.most_common,
    object.__dir__,
    type.__dir__,
    _Duck.__dir__,
}

BUILTIN_GETATTR: set[MayHaveGetattr] = {
    *BUILTIN_GETITEM,
    set,
    frozenset,
    object,
    type,  # `type` handles a lot of generic cases, e.g. numbers as in `int.real`.
    *NUMERICS,
    dict_keys,
    MethodDescriptorType,
    ModuleType,
}


BUILTIN_OPERATIONS = {*BUILTIN_GETATTR}

EVALUATION_POLICIES = {
    "minimal": EvaluationPolicy(
        allow_builtins_access=True,
        allow_locals_access=False,
        allow_globals_access=False,
        allow_item_access=False,
        allow_attr_access=False,
        allowed_calls=set(),
        allow_any_calls=False,
        allow_all_operations=False,
    ),
    "limited": SelectivePolicy(
        allowed_getitem=BUILTIN_GETITEM,
        allowed_getitem_external=SUPPORTED_EXTERNAL_GETITEM,
        allowed_getattr=BUILTIN_GETATTR,
        allowed_getattr_external={
            # pandas Series/Frame implements custom `__getattr__`
            ("pandas", "DataFrame"),
            ("pandas", "Series"),
        },
        allowed_operations=BUILTIN_OPERATIONS,
        allow_builtins_access=True,
        allow_locals_access=True,
        allow_globals_access=True,
        allow_getitem_on_types=True,
        allowed_calls=ALLOWED_CALLS,
    ),
    "unsafe": EvaluationPolicy(
        allow_builtins_access=True,
        allow_locals_access=True,
        allow_globals_access=True,
        allow_attr_access=True,
        allow_item_access=True,
        allow_any_calls=True,
        allow_all_operations=True,
    ),
}


__all__ = [
    "guarded_eval",
    "eval_node",
    "GuardRejection",
    "EvaluationContext",
    "_unbind_method",
]
