from __future__ import annotations

from dask._expr import Expr as BaseExpr
from dask._expr import OptimizerStage, optimize_until
from dask.dataframe.dask_expr._expr import Expr, Projection
from dask.dataframe.dask_expr._merge import Merge
from dask.dataframe.dask_expr.io.parquet import ReadParquet
from dask.utils import funcname, import_required

STAGE_LABELS: dict[OptimizerStage, str] = {
    "logical": "Logical Plan",
    "simplified-logical": "Simplified Logical Plan",
    "tuned-logical": "Tuned Logical Plan",
    "physical": "Physical Plan",
    "simplified-physical": "Simplified Physical Plan",
    "fused": "Fused Physical Plan",
}


def explain(
    expr: BaseExpr, stage: OptimizerStage = "fused", format: str | None = None
) -> None:
    graphviz = import_required(
        "graphviz", "graphviz is a required dependency for using the explain method."
    )

    if format is None:
        format = "png"

    g = graphviz.Digraph(
        STAGE_LABELS[stage], filename=f"explain-{stage}-{expr._name}", format=format
    )
    g.node_attr.update(shape="record")

    expr = optimize_until(expr, stage)

    seen = set(expr._name)
    stack: list[BaseExpr] = [expr]

    while stack:
        node = stack.pop()
        explain_info = _explain_info(node)
        _add_graphviz_node(explain_info, g)
        _add_graphviz_edges(explain_info, g)

        for dep in node.operands:
            if not isinstance(dep, Expr) or dep._name in seen:
                continue
            seen.add(dep._name)
            stack.append(dep)

    g.view()


def _add_graphviz_node(explain_info, graph):
    label = "".join(
        [
            "<{<b>",
            explain_info["label"],
            "</b> | ",
            "<br />".join(
                [f"{key}: {value}" for key, value in explain_info["details"].items()]
            ),
            "}>",
        ]
    )

    graph.node(explain_info["name"], label)


def _add_graphviz_edges(explain_info, graph):
    name = explain_info["name"]
    for _, dep in explain_info["dependencies"]:
        graph.edge(dep, name)


def _explain_info(expr: BaseExpr):
    return {
        "name": expr._name,
        "label": funcname(type(expr)),
        "details": _explain_details(expr),
        "dependencies": _explain_dependencies(expr),
    }


def _explain_details(expr: BaseExpr):

    if not isinstance(expr, Expr):
        return {}
    details = {"npartitions": expr.npartitions}

    if isinstance(expr, Merge):
        details["how"] = expr.how
    elif isinstance(expr, ReadParquet):
        details["path"] = expr.path
    elif isinstance(expr, Projection):
        columns = expr.operand("columns")
        details["ncolumns"] = len(columns) if isinstance(columns, list) else "Series"

    return details


def _explain_dependencies(expr: BaseExpr) -> list[tuple[str, str]]:
    dependencies = []
    for i, operand in enumerate(expr.operands):
        if not isinstance(operand, BaseExpr):
            continue
        param = expr._parameters[i] if i < len(expr._parameters) else ""
        dependencies.append((str(param), operand._name))
    return dependencies
