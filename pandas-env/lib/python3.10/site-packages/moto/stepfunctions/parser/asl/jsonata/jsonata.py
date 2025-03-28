from __future__ import annotations

import json
import re
from typing import Any, Callable, Dict, List, Optional, Set

from moto.stepfunctions.parser.asl.utils.encoding import to_json_str

JSONataExpression = str
VariableReference = str
VariableDeclarations = str

_PATTERN_VARIABLE_REFERENCE: re.Pattern = re.compile(
    r"\$\$|\$[a-zA-Z0-9_$]+(?:\.[a-zA-Z0-9_][a-zA-Z0-9_$]*)*|\$"
)
_ILLEGAL_VARIABLE_REFERENCES: Set[str] = {"$", "$$"}
_VARIABLE_REFERENCE_ASSIGNMENT_OPERATOR: str = ":="
_VARIABLE_REFERENCE_ASSIGNMENT_STOP_SYMBOL: str = ";"
_EXPRESSION_OPEN_SYMBOL: str = "("
_EXPRESSION_CLOSE_SYMBOL: str = ")"


class JSONataException(Exception):
    error: str
    details: Optional[str]

    def __init__(self, error: str, details: Optional[str]):
        self.error = error
        self.details = details


class _JSONataJVMBridge:
    _java_OBJECT_MAPPER: "com.fasterxml.jackson.databind.ObjectMapper"  # noqa
    _java_JSONATA: "com.dashjoin.jsonata.Jsonata.jsonata"  # noqa

    def __init__(self):
        pass

    @staticmethod
    def get() -> _JSONataJVMBridge:
        return _JSONataJVMBridge()

    def eval_jsonata(self, jsonata_expression: JSONataExpression) -> Any:
        try:
            # Evaluate the JSONata expression with the JVM.
            # TODO: Investigate whether it is worth moving this chain of statements (java_*) to a
            #  Java program to reduce i/o between the JVM and this runtime.
            java_expression = self._java_JSONATA(jsonata_expression)
            java_output = java_expression.evaluate(None)
            java_output_string = self._java_OBJECT_MAPPER.writeValueAsString(
                java_output
            )

            # Compute a Python json object from the java string, this is to:
            #  1. Ensure we fully end interactions with the JVM about this value here;
            #  2. The output object may undergo under operations that are not compatible
            #     with jpype objects (such as json.dumps, equality, instanceof, etc.).
            result_str: str = str(java_output_string)
            result_json = json.loads(result_str)

            return result_json
        except Exception as ex:
            raise JSONataException("UNKNOWN", str(ex))


# Final reference to the java evaluation function.
_eval_jsonata: Callable[[JSONataExpression], Any] = None


def eval_jsonata_expression(jsonata_expression: JSONataExpression) -> Any:
    return _eval_jsonata(jsonata_expression)


class IllegalJSONataVariableReference(ValueError):
    variable_reference: VariableReference

    def __init__(self, variable_reference: VariableReference):
        self.variable_reference = variable_reference


def extract_jsonata_variable_references(
    jsonata_expression: JSONataExpression,
) -> set[VariableReference]:
    if not jsonata_expression:
        return set()
    variable_references: List[VariableReference] = _PATTERN_VARIABLE_REFERENCE.findall(
        jsonata_expression
    )
    for variable_reference in variable_references:
        if variable_reference in _ILLEGAL_VARIABLE_REFERENCES:
            raise IllegalJSONataVariableReference(variable_reference=variable_reference)
    return set(variable_references)


def encode_jsonata_variable_declarations(
    bindings: Dict[VariableReference, Any],
) -> VariableDeclarations:
    declarations_parts: List[str] = list()
    for variable_reference, value in bindings.items():
        if isinstance(value, str):
            value_str_lit = f'"{value}"'
        else:
            value_str_lit = to_json_str(value, separators=(",", ":"))
        declarations_parts.extend(
            [
                variable_reference,
                _VARIABLE_REFERENCE_ASSIGNMENT_OPERATOR,
                value_str_lit,
                _VARIABLE_REFERENCE_ASSIGNMENT_STOP_SYMBOL,
            ]
        )
    return "".join(declarations_parts)


def compose_jsonata_expression(
    final_jsonata_expression: JSONataExpression,
    variable_declarations_list: List[VariableDeclarations],
) -> JSONataExpression:
    # TODO: should be expanded to pack the intrinsic functions too.
    variable_declarations = "".join(variable_declarations_list)
    expression = "".join(
        [
            _EXPRESSION_OPEN_SYMBOL,
            variable_declarations,
            final_jsonata_expression,
            _EXPRESSION_CLOSE_SYMBOL,
        ]
    )
    return expression
