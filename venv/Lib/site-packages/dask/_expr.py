from __future__ import annotations

import functools
import os
import weakref
from collections import defaultdict
from collections.abc import Generator
from typing import TYPE_CHECKING, Literal

import toolz

import dask
from dask._task_spec import Task
from dask.tokenize import _tokenize_deterministic
from dask.typing import Key
from dask.utils import funcname, import_required

if TYPE_CHECKING:
    # TODO import from typing (requires Python >=3.10)
    from typing import Any, TypeAlias

OptimizerStage: TypeAlias = Literal[
    "logical",
    "simplified-logical",
    "tuned-logical",
    "physical",
    "simplified-physical",
    "fused",
]


def _unpack_collections(o):
    if isinstance(o, Expr):
        return o

    if hasattr(o, "expr"):
        return o.expr
    else:
        return o


class Expr:
    _parameters: list[str] = []
    _defaults: dict[str, Any] = {}
    _instances: weakref.WeakValueDictionary[str, Expr] = weakref.WeakValueDictionary()

    operands: list

    _determ_token: str | None

    def __new__(cls, *args, _determ_token=None, **kwargs):
        operands = list(args)
        for parameter in cls._parameters[len(operands) :]:
            try:
                operands.append(kwargs.pop(parameter))
            except KeyError:
                operands.append(cls._defaults[parameter])
        assert not kwargs, kwargs
        inst = object.__new__(cls)

        inst._determ_token = _determ_token
        inst.operands = [_unpack_collections(o) for o in operands]
        _name = inst._name
        if _name in Expr._instances:
            return Expr._instances[_name]

        Expr._instances[_name] = inst
        return inst

    def _tune_down(self):
        return None

    def _tune_up(self, parent):
        return None

    def _operands_for_repr(self):
        raise NotImplementedError("Subclasses should implement this method")

    def __str__(self):
        s = ", ".join(self._operands_for_repr())
        return f"{type(self).__name__}({s})"

    def __repr__(self):
        return str(self)

    def _tree_repr_argument_construction(self, i, op, header):
        try:
            param = self._parameters[i]
            default = self._defaults[param]
        except (IndexError, KeyError):
            param = self._parameters[i] if i < len(self._parameters) else ""
            default = "--no-default--"

        if repr(op) != repr(default):
            if param:
                header += f" {param}={repr(op)}"
            else:
                header += repr(op)
        return header

    def _tree_repr_lines(self, indent=0, recursive=True):
        raise NotImplementedError("Subclasses should implement this method")

    def tree_repr(self):
        return os.linesep.join(self._tree_repr_lines())

    def analyze(self, filename: str | None = None, format: str | None = None) -> None:
        from dask.dataframe.dask_expr._expr import Expr as DFExpr
        from dask.dataframe.dask_expr.diagnostics import analyze

        if not isinstance(self, DFExpr):
            raise TypeError(
                "analyze is only supported for dask.dataframe.Expr objects."
            )
        return analyze(self, filename=filename, format=format)

    def explain(
        self, stage: OptimizerStage = "fused", format: str | None = None
    ) -> None:
        from dask.dataframe.dask_expr.diagnostics import explain

        return explain(self, stage, format)

    def pprint(self):
        for line in self._tree_repr_lines():
            print(line)

    def __hash__(self):
        return hash(self._name)

    def __dask_tokenize__(self):
        return self._name

    @staticmethod
    def _reconstruct(*args):
        typ, *operands, token = args
        return typ(*operands, _determ_token=token)

    def __reduce__(self):
        if dask.config.get("dask-expr-no-serialize", False):
            raise RuntimeError(f"Serializing a {type(self)} object")
        return Expr._reconstruct, tuple(
            [type(self)] + self.operands + [self.deterministic_token]
        )

    def _depth(self, cache=None):
        """Depth of the expression tree

        Returns
        -------
        depth: int
        """
        if cache is None:
            cache = {}
        if not self.dependencies():
            return 1
        else:
            result = []
            for expr in self.dependencies():
                if expr._name in cache:
                    result.append(cache[expr._name])
                else:
                    result.append(expr._depth(cache) + 1)
                    cache[expr._name] = result[-1]
            return max(result)

    def __setattr__(self, name: str, value: Any) -> None:
        if name in ["operands", "_determ_token"]:
            object.__setattr__(self, name, value)
            return
        try:
            params = object.__getattribute__(type(self), "_parameters")
            operands = object.__getattribute__(self, "operands")
            operands[params.index(name)] = value
        except ValueError:
            raise AttributeError(
                f"{type(self).__name__} object has no attribute {name}"
            )

    def operand(self, key):
        # Access an operand unambiguously
        # (e.g. if the key is reserved by a method/property)
        return self.operands[type(self)._parameters.index(key)]

    def dependencies(self):
        # Dependencies are `Expr` operands only
        return [operand for operand in self.operands if isinstance(operand, Expr)]

    def _task(self, key: Key, index: int) -> Task:
        """The task for the i'th partition

        Parameters
        ----------
        index:
            The index of the partition of this dataframe

        Examples
        --------
        >>> class Add(Expr):
        ...     def _task(self, i):
        ...         return Task(
        ...            name,
        ...            operator.add,
        ...            TaskRef((self.left._name, i)),
        ...            TaskRef((self.right._name, i))
        ...        )

        Returns
        -------
        task:
            The Dask task to compute this partition

        See Also
        --------
        Expr._layer
        """
        raise NotImplementedError(
            "Expressions should define either _layer (full dictionary) or _task"
            f" (single task).  This expression {type(self)} defines neither"
        )

    def _layer(self) -> dict:
        """The graph layer added by this expression

        Examples
        --------
        >>> class Add(Expr):
        ...     def _layer(self):
        ...         return {
        ...             (self._name, i): (operator.add, (self.left._name, i), (self.right._name, i))
        ...             for i in range(self.npartitions)
        ...         }

        Returns
        -------
        layer: dict
            The Dask task graph added by this expression

        See Also
        --------
        Expr._task
        Expr.__dask_graph__
        """

        return {
            (self._name, i): self._task((self._name, i), i)
            for i in range(self.npartitions)
        }

    def rewrite(self, kind: str, rewritten):
        """Rewrite an expression

        This leverages the ``._{kind}_down`` and ``._{kind}_up``
        methods defined on each class

        Returns
        -------
        expr:
            output expression
        changed:
            whether or not any change occured
        """
        if self._name in rewritten:
            return rewritten[self._name]

        expr = self
        down_name = f"_{kind}_down"
        up_name = f"_{kind}_up"
        while True:
            _continue = False

            # Rewrite this node
            out = getattr(expr, down_name)()
            if out is None:
                out = expr
            if not isinstance(out, Expr):
                return out
            if out._name != expr._name:
                expr = out
                continue

            # Allow children to rewrite their parents
            for child in expr.dependencies():
                out = getattr(child, up_name)(expr)
                if out is None:
                    out = expr
                if not isinstance(out, Expr):
                    return out
                if out is not expr and out._name != expr._name:
                    expr = out
                    _continue = True
                    break

            if _continue:
                continue

            # Rewrite all of the children
            new_operands = []
            changed = False
            for operand in expr.operands:
                if isinstance(operand, Expr):
                    new = operand.rewrite(kind=kind, rewritten=rewritten)
                    rewritten[operand._name] = new
                    if new._name != operand._name:
                        changed = True
                else:
                    new = operand
                new_operands.append(new)

            if changed:
                expr = type(expr)(*new_operands)
                continue
            else:
                break

        return expr

    def simplify_once(self, dependents: defaultdict, simplified: dict):
        """Simplify an expression

        This leverages the ``._simplify_down`` and ``._simplify_up``
        methods defined on each class

        Parameters
        ----------

        dependents: defaultdict[list]
            The dependents for every node.
        simplified: dict
            Cache of simplified expressions for these dependents.

        Returns
        -------
        expr:
            output expression
        """
        # Check if we've already simplified for these dependents
        if self._name in simplified:
            return simplified[self._name]

        expr = self

        while True:
            out = expr._simplify_down()
            if out is None:
                out = expr
            if not isinstance(out, Expr):
                return out
            if out._name != expr._name:
                expr = out

            # Allow children to simplify their parents
            for child in expr.dependencies():
                out = child._simplify_up(expr, dependents)
                if out is None:
                    out = expr

                if not isinstance(out, Expr):
                    return out
                if out is not expr and out._name != expr._name:
                    expr = out
                    break

            # Rewrite all of the children
            new_operands = []
            changed = False
            for operand in expr.operands:
                if isinstance(operand, Expr):
                    # Bandaid for now, waiting for Singleton
                    dependents[operand._name].append(weakref.ref(expr))
                    new = operand.simplify_once(
                        dependents=dependents, simplified=simplified
                    )
                    simplified[operand._name] = new
                    if new._name != operand._name:
                        changed = True
                else:
                    new = operand
                new_operands.append(new)

            if changed:
                expr = type(expr)(*new_operands)

            break

        return expr

    def optimize(self, fuse: bool = False) -> Expr:
        stage: OptimizerStage = "fused" if fuse else "simplified-physical"

        return optimize_until(self, stage)

    def fuse(self) -> Expr:
        return self

    def simplify(self) -> Expr:
        expr = self
        seen = set()
        while True:
            dependents = collect_dependents(expr)
            new = expr.simplify_once(dependents=dependents, simplified={})
            if new._name == expr._name:
                break
            if new._name in seen:
                raise RuntimeError(
                    f"Optimizer does not converge. {expr!r} simplified to {new!r} which was already seen. "
                    "Please report this issue on the dask issue tracker with a minimal reproducer."
                )
            seen.add(new._name)
            expr = new
        return expr

    def _simplify_down(self):
        return

    def _simplify_up(self, parent, dependents):
        return

    def lower_once(self, lowered: dict):
        # Check for a chached result
        try:
            return lowered[self._name]
        except KeyError:
            pass

        expr = self

        # Lower this node
        out = expr._lower()
        if out is None:
            out = expr
        if not isinstance(out, Expr):
            return out

        # Lower all children
        new_operands = []
        changed = False
        for operand in out.operands:
            if isinstance(operand, Expr):
                new = operand.lower_once(lowered)
                if new._name != operand._name:
                    changed = True
            else:
                new = operand
            new_operands.append(new)

        if changed:
            out = type(out)(*new_operands)

        # Cache the result and return
        return lowered.setdefault(self._name, out)

    def lower_completely(self) -> Expr:
        """Lower an expression completely

        This calls the ``lower_once`` method in a loop
        until nothing changes. This function does not
        apply any other optimizations (like ``simplify``).

        Returns
        -------
        expr:
            output expression

        See Also
        --------
        Expr.lower_once
        Expr._lower
        """
        # Lower until nothing changes
        expr = self
        lowered: dict = {}
        while True:
            new = expr.lower_once(lowered)
            if new._name == expr._name:
                break
            expr = new
        return expr

    def _lower(self):
        return

    @functools.cached_property
    def _funcname(self) -> str:
        return funcname(type(self)).lower()

    @property
    def deterministic_token(self):
        if not self._determ_token:
            self._determ_token = _tokenize_deterministic(*self.operands)
        return self._determ_token

    @functools.cached_property
    def _name(self) -> str:
        return self._funcname + "-" + self.deterministic_token

    @property
    def _meta(self):
        raise NotImplementedError()

    def __dask_graph__(self):
        """Traverse expression tree, collect layers"""
        stack = [self]
        seen = set()
        layers = []
        while stack:
            expr = stack.pop()

            if expr._name in seen:
                continue
            seen.add(expr._name)

            layers.append(expr._layer())
            for operand in expr.dependencies():
                stack.append(operand)

        return toolz.merge(layers)

    @property
    def dask(self):
        return self.__dask_graph__()

    def substitute(self, old, new) -> Expr:
        """Substitute a specific term within the expression

        Note that replacing non-`Expr` terms may produce
        unexpected results, and is not recommended.
        Substituting boolean values is not allowed.

        Parameters
        ----------
        old:
            Old term to find and replace.
        new:
            New term to replace instances of `old` with.

        Examples
        --------
        >>> (df + 10).substitute(10, 20)  # doctest: +SKIP
        df + 20
        """
        return self._substitute(old, new, _seen=set())

    def _substitute(self, old, new, _seen):
        if self._name in _seen:
            return self
        # Check if we are replacing a literal
        if isinstance(old, Expr):
            substitute_literal = False
            if self._name == old._name:
                return new
        else:
            substitute_literal = True
            if isinstance(old, bool):
                raise TypeError("Arguments to `substitute` cannot be bool.")

        new_exprs = []
        update = False
        for operand in self.operands:
            if isinstance(operand, Expr):
                val = operand._substitute(old, new, _seen)
                if operand._name != val._name:
                    update = True
                new_exprs.append(val)
            elif (
                "Fused" in type(self).__name__
                and isinstance(operand, list)
                and all(isinstance(op, Expr) for op in operand)
            ):
                # Special handling for `Fused`.
                # We make no promise to dive through a
                # list operand in general, but NEED to
                # do so for the `Fused.exprs` operand.
                val = []
                for op in operand:
                    val.append(op._substitute(old, new, _seen))
                    if val[-1]._name != op._name:
                        update = True
                new_exprs.append(val)
            elif (
                substitute_literal
                and not isinstance(operand, bool)
                and isinstance(operand, type(old))
                and operand == old
            ):
                new_exprs.append(new)
                update = True
            else:
                new_exprs.append(operand)

        if update:  # Only recreate if something changed
            return type(self)(*new_exprs)
        else:
            _seen.add(self._name)
        return self

    def substitute_parameters(self, substitutions: dict) -> Expr:
        """Substitute specific `Expr` parameters

        Parameters
        ----------
        substitutions:
            Mapping of parameter keys to new values. Keys that
            are not found in ``self._parameters`` will be ignored.
        """
        if not substitutions:
            return self

        changed = False
        new_operands = []
        for i, operand in enumerate(self.operands):
            if i < len(self._parameters) and self._parameters[i] in substitutions:
                new_operands.append(substitutions[self._parameters[i]])
                changed = True
            else:
                new_operands.append(operand)
        if changed:
            return type(self)(*new_operands)
        return self

    def _node_label_args(self):
        """Operands to include in the node label by `visualize`"""
        return self.dependencies()

    def _to_graphviz(
        self,
        rankdir="BT",
        graph_attr=None,
        node_attr=None,
        edge_attr=None,
        **kwargs,
    ):
        from dask.dot import label, name

        graphviz = import_required(
            "graphviz",
            "Drawing dask graphs with the graphviz visualization engine requires the `graphviz` "
            "python library and the `graphviz` system library.\n\n"
            "Please either conda or pip install as follows:\n\n"
            "  conda install python-graphviz     # either conda install\n"
            "  python -m pip install graphviz    # or pip install and follow installation instructions",
        )

        graph_attr = graph_attr or {}
        node_attr = node_attr or {}
        edge_attr = edge_attr or {}

        graph_attr["rankdir"] = rankdir
        node_attr["shape"] = "box"
        node_attr["fontname"] = "helvetica"

        graph_attr.update(kwargs)
        g = graphviz.Digraph(
            graph_attr=graph_attr,
            node_attr=node_attr,
            edge_attr=edge_attr,
        )

        stack = [self]
        seen = set()
        dependencies = {}
        while stack:
            expr = stack.pop()

            if expr._name in seen:
                continue
            seen.add(expr._name)

            dependencies[expr] = set(expr.dependencies())
            for dep in expr.dependencies():
                stack.append(dep)

        cache = {}
        for expr in dependencies:
            expr_name = name(expr)
            attrs = {}

            # Make node label
            deps = [
                funcname(type(dep)) if isinstance(dep, Expr) else str(dep)
                for dep in expr._node_label_args()
            ]
            _label = funcname(type(expr))
            if deps:
                _label = f"{_label}({', '.join(deps)})" if deps else _label
            node_label = label(_label, cache=cache)

            attrs.setdefault("label", str(node_label))
            attrs.setdefault("fontsize", "20")
            g.node(expr_name, **attrs)

        for expr, deps in dependencies.items():
            expr_name = name(expr)
            for dep in deps:
                dep_name = name(dep)
                g.edge(dep_name, expr_name)

        return g

    def visualize(self, filename="dask-expr.svg", format=None, **kwargs):
        """
        Visualize the expression graph.
        Requires ``graphviz`` to be installed.

        Parameters
        ----------
        filename : str or None, optional
            The name of the file to write to disk. If the provided `filename`
            doesn't include an extension, '.png' will be used by default.
            If `filename` is None, no file will be written, and the graph is
            rendered in the Jupyter notebook only.
        format : {'png', 'pdf', 'dot', 'svg', 'jpeg', 'jpg'}, optional
            Format in which to write output file. Default is 'svg'.
        **kwargs
           Additional keyword arguments to forward to ``to_graphviz``.
        """
        from dask.dot import graphviz_to_file

        g = self._to_graphviz(**kwargs)
        graphviz_to_file(g, filename, format)
        return g

    def walk(self) -> Generator[Expr]:
        """Iterate through all expressions in the tree

        Returns
        -------
        nodes
            Generator of Expr instances in the graph.
            Ordering is a depth-first search of the expression tree
        """
        stack = [self]
        seen = set()
        while stack:
            node = stack.pop()
            if node._name in seen:
                continue
            seen.add(node._name)

            for dep in node.dependencies():
                stack.append(dep)

            yield node

    def find_operations(self, operation: type | tuple[type]) -> Generator[Expr]:
        """Search the expression graph for a specific operation type

        Parameters
        ----------
        operation
            The operation type to search for.

        Returns
        -------
        nodes
            Generator of `operation` instances. Ordering corresponds
            to a depth-first search of the expression graph.
        """
        assert (
            isinstance(operation, tuple)
            and all(issubclass(e, Expr) for e in operation)
            or issubclass(operation, Expr)  # type: ignore
        ), "`operation` must be`Expr` subclass)"
        return (expr for expr in self.walk() if isinstance(expr, operation))

    def __getattr__(self, key):
        try:
            return object.__getattribute__(self, key)
        except AttributeError as err:
            if key.startswith("_meta"):
                # Avoid a recursive loop if/when `self._meta*`
                # produces an `AttributeError`
                raise RuntimeError(
                    f"Failed to generate metadata for {self}. "
                    "This operation may not be supported by the current backend."
                )

            # Allow operands to be accessed as attributes
            # as long as the keys are not already reserved
            # by existing methods/properties
            _parameters = type(self)._parameters
            if key in _parameters:
                idx = _parameters.index(key)
                return self.operands[idx]

            raise AttributeError(
                f"{err}\n\n"
                "This often means that you are attempting to use an unsupported "
                f"API function.."
            )


def collect_dependents(expr) -> defaultdict:
    dependents = defaultdict(list)
    stack = [expr]
    seen = set()
    while stack:
        node = stack.pop()
        if node._name in seen:
            continue
        seen.add(node._name)

        for dep in node.dependencies():
            stack.append(dep)
            dependents[dep._name].append(weakref.ref(node))
    return dependents


def optimize(expr: Expr, fuse: bool = True) -> Expr:
    """High level query optimization

    This leverages three optimization passes:

    1.  Class based simplification using the ``_simplify`` function and methods
    2.  Blockwise fusion

    Parameters
    ----------
    expr:
        Input expression to optimize
    fuse:
        whether or not to turn on blockwise fusion

    See Also
    --------
    simplify
    optimize_blockwise_fusion
    """
    stage: OptimizerStage = "fused" if fuse else "simplified-physical"

    return optimize_until(expr, stage)


def optimize_until(expr: Expr, stage: OptimizerStage) -> Expr:
    result = expr
    if stage == "logical":
        return result

    # Simplify
    expr = result.simplify()
    if stage == "simplified-logical":
        return expr

    # Manipulate Expression to make it more efficient
    expr = expr.rewrite(kind="tune", rewritten={})
    if stage == "tuned-logical":
        return expr

    # Lower
    expr = expr.lower_completely()
    if stage == "physical":
        return expr

    # Simplify again
    expr = expr.simplify()
    if stage == "simplified-physical":
        return expr

    # Final graph-specific optimizations
    expr = expr.fuse()
    if stage == "fused":
        return expr

    raise ValueError(f"Stage {stage!r} not supported.")
