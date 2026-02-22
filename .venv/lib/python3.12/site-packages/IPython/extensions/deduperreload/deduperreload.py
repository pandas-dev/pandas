from __future__ import annotations
import ast
import builtins
import contextlib
import itertools
import logging
import os
import pickle
import platform
import sys
import textwrap
from types import ModuleType
from typing import TYPE_CHECKING, Any, NamedTuple, cast
from collections.abc import Generator, Iterable

from IPython.extensions.deduperreload.deduperreload_patching import (
    DeduperReloaderPatchingMixin,
)

if TYPE_CHECKING:
    TDefinitionAst = (
        ast.FunctionDef
        | ast.AsyncFunctionDef
        | ast.Import
        | ast.ImportFrom
        | ast.Assign
        | ast.AnnAssign
    )


def get_module_file_name(module: ModuleType | str) -> str | None:
    """Returns the module's file path, or the empty string if it's inaccessible"""
    if (mod := sys.modules.get(module) if isinstance(module, str) else module) is None:
        return ""
    return getattr(mod, "__file__", "") or ""


def compare_ast(node1: ast.AST | list[ast.AST], node2: ast.AST | list[ast.AST]) -> bool:
    """Checks if node1 and node2 have identical AST structure/values, apart from some attributes"""
    if type(node1) is not type(node2):
        return False

    if isinstance(node1, ast.AST):
        for k, v in node1.__dict__.items():
            if k in (
                "lineno",
                "end_lineno",
                "col_offset",
                "end_col_offset",
                "ctx",
                "parent",
            ):
                continue
            if not hasattr(node2, k) or not compare_ast(v, getattr(node2, k)):
                return False
        return True

    elif isinstance(node1, list) and isinstance(  # type:ignore [redundant-expr]
        node2, list
    ):
        return len(node1) == len(node2) and all(
            compare_ast(n1, n2) for n1, n2 in zip(node1, node2)
        )
    else:
        return node1 == node2


class DependencyNode(NamedTuple):
    """
    Each node represents a function.
    qualified_name: string which represents the namespace/name of the function
    abstract_syntax_tree: subtree of the overall module which corresponds to this function

    qualified_name is of the structure: (namespace1, namespace2, ..., name)

    For example, foo() in the following would be represented as (A, B, foo):

    class A:
        class B:
            def foo():
                pass
    """

    qualified_name: tuple[str, ...]
    abstract_syntax_tree: ast.AST


class GatherResult(NamedTuple):
    import_defs: list[tuple[tuple[str, ...], ast.Import | ast.ImportFrom]] = []
    assign_defs: list[tuple[tuple[str, ...], ast.Assign | ast.AnnAssign]] = []
    function_defs: list[
        tuple[tuple[str, ...], ast.FunctionDef | ast.AsyncFunctionDef]
    ] = []
    classes: dict[str, ast.ClassDef] = {}
    unfixable: list[ast.AST] = []

    @classmethod
    def create(cls) -> GatherResult:
        return cls([], [], [], {}, [])

    def all_defs(self) -> Iterable[tuple[tuple[str, ...], TDefinitionAst]]:
        return itertools.chain(self.import_defs, self.assign_defs, self.function_defs)

    def inplace_merge(self, other: GatherResult) -> None:
        self.import_defs.extend(other.import_defs)
        self.assign_defs.extend(other.assign_defs)
        self.function_defs.extend(other.function_defs)
        self.classes.update(other.classes)
        self.unfixable.extend(other.unfixable)


class ConstexprDetector(ast.NodeVisitor):
    def __init__(self) -> None:
        self.is_constexpr = True
        self._allow_builtins_exceptions = True

    @contextlib.contextmanager
    def disallow_builtins_exceptions(self) -> Generator[None, None, None]:
        prev_allow = self._allow_builtins_exceptions
        self._allow_builtins_exceptions = False
        try:
            yield
        finally:
            self._allow_builtins_exceptions = prev_allow

    def visit_Attribute(self, node: ast.Attribute) -> None:
        with self.disallow_builtins_exceptions():
            self.visit(node.value)

    def visit_Name(self, node: ast.Name) -> None:
        if self._allow_builtins_exceptions and hasattr(builtins, node.id):
            return
        self.is_constexpr = False

    def visit(self, node: ast.AST) -> None:
        if not self.is_constexpr:
            # can short-circuit if we've already detected that it's not a constexpr
            return
        super().visit(node)

    def __call__(self, node: ast.AST) -> bool:
        self.is_constexpr = True
        self.visit(node)
        return self.is_constexpr


class AutoreloadTree:
    """
    Recursive data structure to keep track of reloadable functions/methods. Each object corresponds to a specific scope level.
    children: classes inside given scope, maps class name to autoreload tree for that class's scope
    funcs_to_autoreload: list of function names that can be autoreloaded in given scope.
    new_nested_classes: Classes getting added in new autoreload cycle
    """

    def __init__(self) -> None:
        self.children: dict[str, AutoreloadTree] = {}
        self.defs_to_reload: list[tuple[tuple[str, ...], ast.AST]] = []
        self.defs_to_delete: set[str] = set()
        self.new_nested_classes: dict[str, ast.AST] = {}

    def traverse_prefixes(self, prefixes: list[str]) -> AutoreloadTree:
        """
        Return ref to the AutoreloadTree at the namespace specified by prefixes
        """
        cur = self
        for prefix in prefixes:
            if prefix not in cur.children:
                cur.children[prefix] = AutoreloadTree()
            cur = cur.children[prefix]
        return cur


class DeduperReloader(DeduperReloaderPatchingMixin):
    """
    This version of autoreload detects when we can leverage targeted recompilation of a subset of a module and patching
    existing function/method objects to reflect these changes.

    Detects what functions/methods can be reloaded by recursively comparing the old/new AST of module-level classes,
    module-level classes' methods, recursing through nested classes' methods. If other changes are made, original
    autoreload algorithm is called directly.
    """

    def __init__(self) -> None:
        self._to_autoreload: AutoreloadTree = AutoreloadTree()
        self.source_by_modname: dict[str, str] = {}
        self.dependency_graph: dict[tuple[str, ...], list[DependencyNode]] = {}
        self._enabled = True

    @property
    def enabled(self) -> bool:
        return self._enabled and platform.python_implementation() == "CPython"

    @enabled.setter
    def enabled(self, value: bool) -> None:
        self._enabled = value

    def update_sources(self) -> None:
        """
        Update dictionary source_by_modname with current modules' source codes.
        """
        if not self.enabled:
            return
        for new_modname in sys.modules.keys() - self.source_by_modname.keys():
            new_module = sys.modules[new_modname]
            if (
                (fname := get_module_file_name(new_module)) is None
                or "site-packages" in fname
                or "dist-packages" in fname
                or not os.access(fname, os.R_OK)
            ):
                self.source_by_modname[new_modname] = ""
                continue
            # Skip binary files (e.g., .so, .pyd, .pyo)
            if not fname.endswith(".py"):
                self.source_by_modname[new_modname] = ""
                continue
            try:
                with open(fname, "r") as f:
                    self.source_by_modname[new_modname] = f.read()
            except Exception as e:
                logger = logging.getLogger("autoreload")
                logger.exception(
                    f"Failed to read module file '{fname}' for module '{new_modname}': {type(e).__name__}"
                )
                self.source_by_modname[new_modname] = ""

    constexpr_detector = ConstexprDetector()

    @staticmethod
    def is_enum_subclass(node: ast.Module | ast.ClassDef) -> bool:
        if isinstance(node, ast.Module):
            return False
        for base in node.bases:
            if isinstance(base, ast.Name) and base.id == "Enum":
                return True
            elif (
                isinstance(base, ast.Attribute)
                and base.attr == "Enum"
                and isinstance(base.value, ast.Name)
                and base.value.id == "enum"
            ):
                return True
        return False

    @classmethod
    def is_constexpr_assign(
        cls, node: ast.AST, parent_node: ast.Module | ast.ClassDef
    ) -> bool:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)) or node.value is None:
            return False
        if cls.is_enum_subclass(parent_node):
            return False
        for target in node.targets if isinstance(node, ast.Assign) else [node.target]:
            if not isinstance(target, ast.Name):
                return False
        return cls.constexpr_detector(node.value)

    @classmethod
    def _gather_children(
        cls, body: list[ast.stmt], parent_node: ast.Module | ast.ClassDef
    ) -> GatherResult:
        """
        Given list of ast elements, return:
        1. dict mapping function names to their ASTs.
        2. dict mapping class names to their ASTs.
        3. list of any other ASTs.
        """
        result = GatherResult.create()
        for ast_node in body:
            ast_elt: ast.expr | ast.stmt = ast_node
            while isinstance(ast_elt, ast.Expr):
                ast_elt = ast_elt.value
            if isinstance(ast_elt, (ast.FunctionDef, ast.AsyncFunctionDef)):
                result.function_defs.append(((ast_elt.name,), ast_elt))
            elif isinstance(ast_elt, (ast.Import, ast.ImportFrom)):
                result.import_defs.append(
                    (tuple(name.asname or name.name for name in ast_elt.names), ast_elt)
                )
            elif isinstance(ast_elt, ast.ClassDef):
                result.classes[ast_elt.name] = ast_elt
            elif isinstance(ast_elt, ast.If):
                result.unfixable.append(ast_elt.test)
                result.inplace_merge(cls._gather_children(ast_elt.body, parent_node))
                result.inplace_merge(cls._gather_children(ast_elt.orelse, parent_node))
            elif isinstance(ast_elt, (ast.AsyncWith, ast.With)):
                result.unfixable.extend(ast_elt.items)
                result.inplace_merge(cls._gather_children(ast_elt.body, parent_node))
            elif isinstance(ast_elt, ast.Try):
                result.inplace_merge(cls._gather_children(ast_elt.body, parent_node))
                result.inplace_merge(cls._gather_children(ast_elt.orelse, parent_node))
                result.inplace_merge(
                    cls._gather_children(ast_elt.finalbody, parent_node)
                )
                for handler in ast_elt.handlers:
                    if handler.type is not None:
                        result.unfixable.append(handler.type)
                    result.inplace_merge(
                        cls._gather_children(handler.body, parent_node)
                    )
            elif not isinstance(ast_elt, (ast.Constant, ast.Pass)):
                if cls.is_constexpr_assign(ast_elt, parent_node):
                    assert isinstance(ast_elt, (ast.Assign, ast.AnnAssign))
                    targets = (
                        ast_elt.targets
                        if isinstance(ast_elt, ast.Assign)
                        else [ast_elt.target]
                    )
                    result.assign_defs.append(
                        (
                            tuple(cast(ast.Name, target).id for target in targets),
                            ast_elt,
                        )
                    )
                else:
                    result.unfixable.append(ast_elt)
        return result

    def detect_autoreload(
        self,
        old_node: ast.Module | ast.ClassDef,
        new_node: ast.Module | ast.ClassDef,
        prefixes: list[str] | None = None,
    ) -> bool:
        """
        Returns
        -------
        `True` if we can run our targeted autoreload algorithm safely.
        `False` if we should instead use IPython's original autoreload implementation.
        """
        if not self.enabled:
            return False
        prefixes = prefixes or []

        old_result = self._gather_children(old_node.body, old_node)
        new_result = self._gather_children(new_node.body, new_node)
        old_defs_by_name: dict[str, ast.AST] = {
            name: ast_def for names, ast_def in old_result.all_defs() for name in names
        }
        new_defs_by_name: dict[str, ast.AST] = {
            name: ast_def for names, ast_def in new_result.all_defs() for name in names
        }

        if not compare_ast(old_result.unfixable, new_result.unfixable):
            return False

        cur = self._to_autoreload.traverse_prefixes(prefixes)
        for names, new_ast_def in new_result.all_defs():
            names_to_reload = []
            for name in names:
                if new_defs_by_name[name] is not new_ast_def:
                    continue
                if name not in old_defs_by_name or not compare_ast(
                    new_ast_def, old_defs_by_name[name]
                ):
                    names_to_reload.append(name)
            if names_to_reload:
                cur.defs_to_reload.append((tuple(names), new_ast_def))
        cur.defs_to_delete |= set(old_defs_by_name.keys()) - set(
            new_defs_by_name.keys()
        )
        for name, new_ast_def_class in new_result.classes.items():
            if name not in old_result.classes:
                cur.new_nested_classes[name] = new_ast_def_class
            elif not compare_ast(
                new_ast_def_class, old_result.classes[name]
            ) and not self.detect_autoreload(
                old_result.classes[name], new_ast_def_class, prefixes + [name]
            ):
                return False
        return True

    def _check_dependents(self) -> bool:
        """
        If a decorator function is modified, we should similarly reload the functions which are decorated by this
        decorator. Iterate through the Dependency Graph to find such cases in the given AutoreloadTree.
        """
        for node in self._check_dependents_inner():
            self._add_node_to_autoreload_tree(node)
        return True

    def _add_node_to_autoreload_tree(self, node: DependencyNode) -> None:
        """
        Given a node of the dependency graph, add decorator dependencies to the autoreload tree.
        """
        if len(node.qualified_name) == 0:
            return
        cur = self._to_autoreload.traverse_prefixes(list(node.qualified_name[:-1]))
        if node.abstract_syntax_tree is not None:
            cur.defs_to_reload.append(
                ((node.qualified_name[-1],), node.abstract_syntax_tree)
            )

    def _check_dependents_inner(
        self, prefixes: list[str] | None = None
    ) -> list[DependencyNode]:
        prefixes = prefixes or []
        cur = self._to_autoreload.traverse_prefixes(prefixes)
        ans = []
        for (func_name, *_), _ in cur.defs_to_reload:
            node = tuple(prefixes + [func_name])
            ans.extend(self._gen_dependents(node))
        for class_name in cur.new_nested_classes:
            ans.extend(self._check_dependents_inner(prefixes + [class_name]))
        return ans

    def _gen_dependents(self, qualname: tuple[str, ...]) -> list[DependencyNode]:
        ans = []
        if qualname not in self.dependency_graph:
            return []
        for elt in self.dependency_graph[qualname]:
            ans.extend(self._gen_dependents(elt.qualified_name))
            ans.append(elt)
        return ans

    def _patch_namespace_inner(
        self, ns: ModuleType | type, prefixes: list[str] | None = None
    ) -> bool:
        """
        This function patches module functions and methods. Specifically, only objects with their name in
        self.to_autoreload will be considered for patching. If an object has been marked to be autoreloaded,
        new_source_code gets executed in the old version's global environment. Then, replace the old function's
        attributes with the new function's attributes.
        """
        prefixes = prefixes or []
        cur = self._to_autoreload.traverse_prefixes(prefixes)
        namespace_to_check = ns
        for prefix in prefixes:
            namespace_to_check = namespace_to_check.__dict__[prefix]
        seen_names: set[str] = set()
        for names, new_ast_def in cur.defs_to_reload:
            if len(names) == 1 and names[0] in seen_names:
                continue
            seen_names.update(names)
            local_env: dict[str, Any] = {}
            if (
                isinstance(new_ast_def, (ast.FunctionDef, ast.AsyncFunctionDef))
                and (name := names[0]) in namespace_to_check.__dict__
            ):
                assert len(names) == 1
                to_patch_to = namespace_to_check.__dict__[name]
                if isinstance(to_patch_to, (staticmethod, classmethod)):
                    to_patch_to = to_patch_to.__func__
                # exec new source code using old function's (obj) globals environment.
                func_code = textwrap.dedent(ast.unparse(new_ast_def))
                if is_method := (len(prefixes) > 0):
                    func_code = "class __autoreload_class__:\n" + textwrap.indent(
                        func_code, "    "
                    )
                global_env = ns.__dict__
                if not isinstance(global_env, dict):
                    global_env = dict(global_env)
                # Compile with correct filename to preserve in traceback
                filename = (
                    getattr(to_patch_to, "__code__", None)
                    and to_patch_to.__code__.co_filename
                    or "<string>"
                )
                func_asts = [ast.parse(func_code)]
                if len(cast(ast.FunctionDef, func_asts[0].body[0]).decorator_list) > 0:
                    without_decorator_list = pickle.loads(pickle.dumps(func_asts[0]))
                    cast(
                        ast.FunctionDef, without_decorator_list.body[0]
                    ).decorator_list = []
                    func_asts.insert(0, without_decorator_list)
                for func_ast in func_asts:
                    compiled_code = compile(
                        func_ast, filename, mode="exec", dont_inherit=True
                    )
                    exec(compiled_code, global_env, local_env)  # type: ignore[arg-type]
                    # local_env contains the function exec'd from  new version of function
                    if is_method:
                        to_patch_from = getattr(local_env["__autoreload_class__"], name)
                    else:
                        to_patch_from = local_env[name]
                    if isinstance(to_patch_from, (staticmethod, classmethod)):
                        to_patch_from = to_patch_from.__func__
                    if isinstance(to_patch_to, property) and isinstance(
                        to_patch_from, property
                    ):
                        for attr in ("fget", "fset", "fdel"):
                            if (
                                getattr(to_patch_to, attr) is None
                                or getattr(to_patch_from, attr) is None
                            ):
                                self.try_patch_attr(to_patch_to, to_patch_from, attr)
                            else:
                                self.patch_function(
                                    getattr(to_patch_to, attr),
                                    getattr(to_patch_from, attr),
                                    is_method,
                                )
                    elif not isinstance(to_patch_to, property) and not isinstance(
                        to_patch_from, property
                    ):
                        self.patch_function(to_patch_to, to_patch_from, is_method)
                    else:
                        raise ValueError(
                            "adding or removing property decorations not supported"
                        )
            else:
                exec(
                    ast.unparse(new_ast_def),
                    ns.__dict__ | namespace_to_check.__dict__,
                    local_env,
                )
                for name in names:
                    setattr(namespace_to_check, name, local_env[name])
        cur.defs_to_reload.clear()
        for name in cur.defs_to_delete:
            try:
                delattr(namespace_to_check, name)
            except (AttributeError, TypeError, ValueError):
                # give up on deleting the attribute, let the stale one dangle
                pass
        cur.defs_to_delete.clear()
        for class_name, class_ast_node in cur.new_nested_classes.items():
            local_env_class: dict[str, Any] = {}
            exec(
                ast.unparse(class_ast_node),
                ns.__dict__ | namespace_to_check.__dict__,
                local_env_class,
            )
            setattr(namespace_to_check, class_name, local_env_class[class_name])
        cur.new_nested_classes.clear()
        for class_name in cur.children.keys():
            if not self._patch_namespace(ns, prefixes + [class_name]):
                return False
        cur.children.clear()
        return True

    def _patch_namespace(
        self, ns: ModuleType | type, prefixes: list[str] | None = None
    ) -> bool:
        """
        Wrapper for patching all elements in a namespace as specified by the to_autoreload member variable.
        Returns `true` if patching was successful, and `false` if unsuccessful.
        """
        try:
            return self._patch_namespace_inner(ns, prefixes=prefixes)
        except Exception:
            return False

    def maybe_reload_module(self, module: ModuleType) -> bool:
        """
        Uses Deduperreload to try to update a module.
        Returns `true` on success and `false` on failure.
        """
        if not self.enabled:
            return False
        if not (modname := getattr(module, "__name__", None)):
            return False
        if (fname := get_module_file_name(module)) is None:
            return False
        try:
            with open(fname, "r") as f:
                new_source_code = f.read()
        except Exception as e:
            logger = logging.getLogger("autoreload")
            logger.exception(
                f"Failed to read module file '{fname}' for module '{modname}': {type(e).__name__}"
            )
            return False
        patched_flag = False
        if old_source_code := self.source_by_modname.get(modname):
            # get old/new module ast
            try:
                old_module_ast = ast.parse(old_source_code)
                new_module_ast = ast.parse(new_source_code)
            except Exception:
                return False
            # detect if we are able to use our autoreload algorithm
            ctx = contextlib.suppress()
            with ctx:
                self._build_dependency_graph(new_module_ast)
                if (
                    self.detect_autoreload(old_module_ast, new_module_ast)
                    and self._check_dependents()
                    and self._patch_namespace(module)
                ):
                    patched_flag = True

        self.source_by_modname[modname] = new_source_code
        self._to_autoreload = AutoreloadTree()
        return patched_flag

    def _separate_name(
        self,
        decorator: ast.Attribute | ast.Name | ast.Call | ast.expr,
        accept_calls: bool,
    ) -> list[str] | None:
        """
        Generates a qualified name for a given decorator by finding its relative namespace.
        """
        if isinstance(decorator, ast.Name):
            return [decorator.id]
        elif isinstance(decorator, ast.Call):
            if accept_calls:
                return self._separate_name(decorator.func, False)
            else:
                return None
        if not isinstance(decorator, ast.Attribute):
            return None
        if pref := self._separate_name(decorator.value, False):
            return pref + [decorator.attr]
        else:
            return None

    def _gather_dependents(
        self, body: list[ast.stmt], body_prefixes: list[str] | None = None
    ) -> bool:
        body_prefixes = body_prefixes or []
        for ast_node in body:
            ast_elt: ast.expr | ast.stmt = ast_node
            if isinstance(ast_elt, ast.ClassDef):
                self._gather_dependents(ast_elt.body, body_prefixes + [ast_elt.name])
                continue
            if not isinstance(ast_elt, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            qualified_name = tuple(body_prefixes + [ast_elt.name])
            cur_dependency_node = DependencyNode(qualified_name, ast_elt)
            for decorator in ast_elt.decorator_list:
                decorator_path = self._separate_name(decorator, True)
                if not decorator_path:
                    continue
                decorator_path_tuple = tuple(decorator_path)
                self.dependency_graph.setdefault(decorator_path_tuple, []).append(
                    cur_dependency_node
                )
        return True

    def _build_dependency_graph(self, new_ast: ast.Module | ast.ClassDef) -> bool:
        """
        Wrapper function for generating dependency graph given some AST.
        Returns `true` on success. Returns `false` on failure.
        Currently, only returns `true` as we do not block on failure to build this graph.
        """
        return self._gather_dependents(new_ast.body)
