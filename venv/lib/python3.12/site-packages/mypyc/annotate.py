"""Generate source code formatted as HTML, with bottlenecks annotated and highlighted.

Various heuristics are used to detect common issues that cause slower than
expected performance.
"""

from __future__ import annotations

import os.path
import sys
from html import escape
from typing import Final

from mypy.build import BuildResult
from mypy.nodes import (
    AssignmentStmt,
    CallExpr,
    ClassDef,
    Decorator,
    DictionaryComprehension,
    Expression,
    ForStmt,
    FuncDef,
    GeneratorExpr,
    IndexExpr,
    LambdaExpr,
    MemberExpr,
    MypyFile,
    NamedTupleExpr,
    NameExpr,
    NewTypeExpr,
    Node,
    OpExpr,
    RefExpr,
    TupleExpr,
    TypedDictExpr,
    TypeInfo,
    TypeVarExpr,
    Var,
    WithStmt,
)
from mypy.traverser import TraverserVisitor
from mypy.types import AnyType, Instance, ProperType, Type, TypeOfAny, get_proper_type
from mypy.util import FancyFormatter
from mypyc.ir.func_ir import FuncIR
from mypyc.ir.module_ir import ModuleIR
from mypyc.ir.ops import CallC, LoadLiteral, LoadStatic, Value
from mypyc.irbuild.mapper import Mapper


class Annotation:
    """HTML annotation for compiled source code"""

    def __init__(self, message: str, priority: int = 1) -> None:
        # Message as HTML that describes an issue and/or how to fix it.
        # Multiple messages on a line may be concatenated.
        self.message = message
        # If multiple annotations are generated for a single line, only report
        # the highest-priority ones. Some use cases generate multiple annotations,
        # and this can be used to reduce verbosity by hiding the lower-priority
        # ones.
        self.priority = priority


op_hints: Final = {
    "PyNumber_Add": Annotation('Generic "+" operation.'),
    "PyNumber_Subtract": Annotation('Generic "-" operation.'),
    "PyNumber_Multiply": Annotation('Generic "*" operation.'),
    "PyNumber_TrueDivide": Annotation('Generic "/" operation.'),
    "PyNumber_FloorDivide": Annotation('Generic "//" operation.'),
    "PyNumber_Positive": Annotation('Generic unary "+" operation.'),
    "PyNumber_Negative": Annotation('Generic unary "-" operation.'),
    "PyNumber_And": Annotation('Generic "&" operation.'),
    "PyNumber_Or": Annotation('Generic "|" operation.'),
    "PyNumber_Xor": Annotation('Generic "^" operation.'),
    "PyNumber_Lshift": Annotation('Generic "<<" operation.'),
    "PyNumber_Rshift": Annotation('Generic ">>" operation.'),
    "PyNumber_Invert": Annotation('Generic "~" operation.'),
    "PyObject_Call": Annotation("Generic call operation."),
    "PyObject_RichCompare": Annotation("Generic comparison operation."),
    "PyObject_GetItem": Annotation("Generic indexing operation."),
    "PyObject_SetItem": Annotation("Generic indexed assignment."),
}

stdlib_hints: Final = {
    "functools.partial": Annotation(
        '"functools.partial" is inefficient in compiled code.', priority=3
    ),
    "itertools.chain": Annotation(
        '"itertools.chain" is inefficient in compiled code (hint: replace with for loops).',
        priority=3,
    ),
    "itertools.groupby": Annotation(
        '"itertools.groupby" is inefficient in compiled code.', priority=3
    ),
    "itertools.islice": Annotation(
        '"itertools.islice" is inefficient in compiled code (hint: replace with for loop over index range).',
        priority=3,
    ),
    "copy.deepcopy": Annotation(
        '"copy.deepcopy" tends to be slow. Make a shallow copy if possible.', priority=2
    ),
}

CSS = """\
.collapsible {
    cursor: pointer;
}

.content {
    display: block;
    margin-top: 10px;
    margin-bottom: 10px;
}

.hint {
    display: inline;
    border: 1px solid #ccc;
    padding: 5px;
}
"""

JS = """\
document.querySelectorAll('.collapsible').forEach(function(collapsible) {
    collapsible.addEventListener('click', function() {
        const content = this.nextElementSibling;
        if (content.style.display === 'none') {
            content.style.display = 'block';
        } else {
            content.style.display = 'none';
        }
    });
});
"""


class AnnotatedSource:
    """Annotations for a single compiled source file."""

    def __init__(self, path: str, annotations: dict[int, list[Annotation]]) -> None:
        self.path = path
        self.annotations = annotations


def generate_annotated_html(
    html_fnam: str, result: BuildResult, modules: dict[str, ModuleIR], mapper: Mapper
) -> None:
    annotations = []
    for mod, mod_ir in modules.items():
        path = result.graph[mod].path
        tree = result.graph[mod].tree
        assert tree is not None
        annotations.append(
            generate_annotations(path or "<source>", tree, mod_ir, result.types, mapper)
        )
    html = generate_html_report(annotations)
    with open(html_fnam, "w") as f:
        f.write(html)

    formatter = FancyFormatter(sys.stdout, sys.stderr, False)
    formatted = formatter.style(os.path.abspath(html_fnam), "none", underline=True, bold=True)
    print(f"\nWrote {formatted} -- open in browser to view\n")


def generate_annotations(
    path: str, tree: MypyFile, ir: ModuleIR, type_map: dict[Expression, Type], mapper: Mapper
) -> AnnotatedSource:
    anns = {}
    for func_ir in ir.functions:
        anns.update(function_annotations(func_ir, tree))
    visitor = ASTAnnotateVisitor(type_map, mapper)
    for defn in tree.defs:
        defn.accept(visitor)
    anns.update(visitor.anns)
    for line in visitor.ignored_lines:
        if line in anns:
            del anns[line]
    return AnnotatedSource(path, anns)


def function_annotations(func_ir: FuncIR, tree: MypyFile) -> dict[int, list[Annotation]]:
    """Generate annotations based on mypyc IR."""
    # TODO: check if func_ir.line is -1
    anns: dict[int, list[Annotation]] = {}
    for block in func_ir.blocks:
        for op in block.ops:
            if isinstance(op, CallC):
                name = op.function_name
                ann: str | Annotation | None = None
                if name == "CPyObject_GetAttr":
                    attr_name = get_str_literal(op.args[1])
                    if attr_name in ("__prepare__", "GeneratorExit", "StopIteration"):
                        # These attributes are internal to mypyc/CPython, and/or accessed
                        # implicitly in generated code. The user has little control over
                        # them.
                        ann = None
                    elif attr_name:
                        ann = f'Get non-native attribute "{attr_name}".'
                    else:
                        ann = "Dynamic attribute lookup."
                elif name == "PyObject_SetAttr":
                    attr_name = get_str_literal(op.args[1])
                    if attr_name == "__mypyc_attrs__":
                        # This is set implicitly and can't be avoided.
                        ann = None
                    elif attr_name:
                        ann = f'Set non-native attribute "{attr_name}".'
                    else:
                        ann = "Dynamic attribute set."
                elif name == "PyObject_VectorcallMethod":
                    method_name = get_str_literal(op.args[0])
                    if method_name:
                        ann = f'Call non-native method "{method_name}" (it may be defined in a non-native class, or decorated).'
                    else:
                        ann = "Dynamic method call."
                elif name in op_hints:
                    ann = op_hints[name]
                elif name in ("CPyDict_GetItem", "CPyDict_SetItem"):
                    if (
                        isinstance(op.args[0], LoadStatic)
                        and isinstance(op.args[1], LoadLiteral)
                        and func_ir.name != "__top_level__"
                    ):
                        load = op.args[0]
                        name = str(op.args[1].value)
                        sym = tree.names.get(name)
                        if (
                            sym
                            and sym.node
                            and load.namespace == "static"
                            and load.identifier == "globals"
                        ):
                            if sym.node.fullname in stdlib_hints:
                                ann = stdlib_hints[sym.node.fullname]
                            elif isinstance(sym.node, Var):
                                ann = (
                                    f'Access global "{name}" through namespace '
                                    + "dictionary (hint: access is faster if you can make it Final)."
                                )
                            else:
                                ann = f'Access "{name}" through global namespace dictionary.'
                if ann:
                    if isinstance(ann, str):
                        ann = Annotation(ann)
                    anns.setdefault(op.line, []).append(ann)
    return anns


class ASTAnnotateVisitor(TraverserVisitor):
    """Generate annotations from mypy AST and inferred types."""

    def __init__(self, type_map: dict[Expression, Type], mapper: Mapper) -> None:
        self.anns: dict[int, list[Annotation]] = {}
        self.ignored_lines: set[int] = set()
        self.func_depth = 0
        self.type_map = type_map
        self.mapper = mapper

    def visit_func_def(self, o: FuncDef, /) -> None:
        if self.func_depth > 0:
            self.annotate(
                o,
                "A nested function object is allocated each time statement is executed. "
                + "A module-level function would be faster.",
            )
        self.func_depth += 1
        super().visit_func_def(o)
        self.func_depth -= 1

    def visit_for_stmt(self, o: ForStmt, /) -> None:
        self.check_iteration([o.expr], "For loop")
        super().visit_for_stmt(o)

    def visit_dictionary_comprehension(self, o: DictionaryComprehension, /) -> None:
        self.check_iteration(o.sequences, "Comprehension")
        super().visit_dictionary_comprehension(o)

    def visit_generator_expr(self, o: GeneratorExpr, /) -> None:
        self.check_iteration(o.sequences, "Comprehension or generator")
        super().visit_generator_expr(o)

    def check_iteration(self, expressions: list[Expression], kind: str) -> None:
        for expr in expressions:
            typ = self.get_type(expr)
            if isinstance(typ, AnyType):
                self.annotate(expr, f'{kind} uses generic operations (iterable has type "Any").')
            elif isinstance(typ, Instance) and typ.type.fullname in (
                "typing.Iterable",
                "typing.Iterator",
                "typing.Sequence",
                "typing.MutableSequence",
            ):
                self.annotate(
                    expr,
                    f'{kind} uses generic operations (iterable has the abstract type "{typ.type.fullname}").',
                )

    def visit_class_def(self, o: ClassDef, /) -> None:
        super().visit_class_def(o)
        if self.func_depth == 0:
            # Don't complain about base classes at top level
            for base in o.base_type_exprs:
                self.ignored_lines.add(base.line)

            for s in o.defs.body:
                if isinstance(s, AssignmentStmt):
                    # Don't complain about attribute initializers
                    self.ignored_lines.add(s.line)
                elif isinstance(s, Decorator):
                    # Don't complain about decorator definitions that generate some
                    # dynamic operations. This is a bit heavy-handed.
                    self.ignored_lines.add(s.func.line)

    def visit_with_stmt(self, o: WithStmt, /) -> None:
        for expr in o.expr:
            if isinstance(expr, CallExpr) and isinstance(expr.callee, RefExpr):
                node = expr.callee.node
                if isinstance(node, Decorator):
                    if any(
                        isinstance(d, RefExpr)
                        and d.node
                        and d.node.fullname == "contextlib.contextmanager"
                        for d in node.decorators
                    ):
                        self.annotate(
                            expr,
                            f'"{node.name}" uses @contextmanager, which is slow '
                            + "in compiled code. Use a native class with "
                            + '"__enter__" and "__exit__" methods instead.',
                            priority=3,
                        )
        super().visit_with_stmt(o)

    def visit_assignment_stmt(self, o: AssignmentStmt, /) -> None:
        special_form = False
        if self.func_depth == 0:
            analyzed: Expression | None = o.rvalue
            if isinstance(o.rvalue, (CallExpr, IndexExpr, OpExpr)):
                analyzed = o.rvalue.analyzed
            if o.is_alias_def or isinstance(
                analyzed, (TypeVarExpr, NamedTupleExpr, TypedDictExpr, NewTypeExpr)
            ):
                special_form = True
            if special_form:
                # TODO: Ignore all lines if multi-line
                self.ignored_lines.add(o.line)
        super().visit_assignment_stmt(o)

    def visit_name_expr(self, o: NameExpr, /) -> None:
        if ann := stdlib_hints.get(o.fullname):
            self.annotate(o, ann)

    def visit_member_expr(self, o: MemberExpr, /) -> None:
        super().visit_member_expr(o)
        if ann := stdlib_hints.get(o.fullname):
            self.annotate(o, ann)

    def visit_call_expr(self, o: CallExpr, /) -> None:
        super().visit_call_expr(o)
        if (
            isinstance(o.callee, RefExpr)
            and o.callee.fullname == "builtins.isinstance"
            and len(o.args) == 2
        ):
            arg = o.args[1]
            self.check_isinstance_arg(arg)
        elif isinstance(o.callee, RefExpr) and isinstance(o.callee.node, TypeInfo):
            info = o.callee.node
            class_ir = self.mapper.type_to_ir.get(info)
            if (class_ir and not class_ir.is_ext_class) or (
                class_ir is None and not info.fullname.startswith("builtins.")
            ):
                self.annotate(
                    o, f'Creating an instance of non-native class "{info.name}" ' + "is slow.", 2
                )
            elif class_ir and class_ir.is_augmented:
                self.annotate(
                    o,
                    f'Class "{info.name}" is only partially native, and '
                    + "constructing an instance is slow.",
                    2,
                )
        elif isinstance(o.callee, RefExpr) and isinstance(o.callee.node, Decorator):
            decorator = o.callee.node
            if self.mapper.is_native_ref_expr(o.callee):
                self.annotate(
                    o,
                    f'Calling a decorated function ("{decorator.name}") is inefficient, even if it\'s native.',
                    2,
                )

    def check_isinstance_arg(self, arg: Expression) -> None:
        if isinstance(arg, RefExpr):
            if isinstance(arg.node, TypeInfo) and arg.node.is_protocol:
                self.annotate(
                    arg, f'Expensive isinstance() check against protocol "{arg.node.name}".'
                )
        elif isinstance(arg, TupleExpr):
            for item in arg.items:
                self.check_isinstance_arg(item)

    def visit_lambda_expr(self, o: LambdaExpr, /) -> None:
        self.annotate(
            o,
            "A new object is allocated for lambda each time it is evaluated. "
            + "A module-level function would be faster.",
        )
        super().visit_lambda_expr(o)

    def annotate(self, o: Node, ann: str | Annotation, priority: int = 1) -> None:
        if isinstance(ann, str):
            ann = Annotation(ann, priority=priority)
        self.anns.setdefault(o.line, []).append(ann)

    def get_type(self, e: Expression) -> ProperType:
        t = self.type_map.get(e)
        if t:
            return get_proper_type(t)
        return AnyType(TypeOfAny.unannotated)


def get_str_literal(v: Value) -> str | None:
    if isinstance(v, LoadLiteral) and isinstance(v.value, str):
        return v.value
    return None


def get_max_prio(anns: list[Annotation]) -> list[Annotation]:
    max_prio = max(a.priority for a in anns)
    return [a for a in anns if a.priority == max_prio]


def generate_html_report(sources: list[AnnotatedSource]) -> str:
    html = []
    html.append("<html>\n<head>\n")
    html.append(f"<style>\n{CSS}\n</style>")
    html.append("</head>\n")
    html.append("<body>\n")
    for src in sources:
        html.append(f"<h2><tt>{src.path}</tt></h2>\n")
        html.append("<pre>")
        src_anns = src.annotations
        with open(src.path) as f:
            lines = f.readlines()
        for i, s in enumerate(lines):
            s = escape(s)
            line = i + 1
            linenum = "%5d" % line
            if line in src_anns:
                anns = get_max_prio(src_anns[line])
                ann_strs = [a.message for a in anns]
                hint = " ".join(ann_strs)
                s = colorize_line(linenum, s, hint_html=hint)
            else:
                s = linenum + "  " + s
            html.append(s)
        html.append("</pre>")

    html.append("<script>")
    html.append(JS)
    html.append("</script>")

    html.append("</body></html>\n")
    return "".join(html)


def colorize_line(linenum: str, s: str, hint_html: str) -> str:
    hint_prefix = " " * len(linenum) + "  "
    line_span = f'<div class="collapsible" style="background-color: #fcc">{linenum}  {s}</div>'
    hint_div = f'<div class="content">{hint_prefix}<div class="hint">{hint_html}</div></div>'
    return f"<span>{line_span}{hint_div}</span>"
