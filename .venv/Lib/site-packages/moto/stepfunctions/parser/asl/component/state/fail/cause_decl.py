import abc
import copy
from typing import Set

from moto.stepfunctions.parser.asl.component.common.jsonata.jsonata_template_value_terminal import (
    JSONataTemplateValueTerminalExpression,
)
from moto.stepfunctions.parser.asl.component.common.variable_sample import (
    VariableSample,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.component.intrinsic.function.function import Function
from moto.stepfunctions.parser.asl.component.intrinsic.functionname.state_function_name_types import (
    StatesFunctionNameType,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.parse.intrinsic.intrinsic_parser import (
    IntrinsicParser,
)
from moto.stepfunctions.parser.asl.utils.json_path import extract_json


class CauseDecl(EvalComponent, abc.ABC): ...


class CauseConst(CauseDecl):
    value: str

    def __init__(self, value: str):
        self.value = value

    def _eval_body(self, env: Environment) -> None:
        env.stack.append(self.value)


class CauseJSONata(CauseDecl):
    def __init__(
        self,
        jsonata_template_value_terminal_expression: JSONataTemplateValueTerminalExpression,
    ):
        super().__init__()
        self.jsonata_template_value_terminal_expression = (
            jsonata_template_value_terminal_expression
        )

    def _eval_body(self, env: Environment) -> None:
        self.jsonata_template_value_terminal_expression.eval(env=env)


class CauseVar(CauseDecl):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        self.variable_sample = variable_sample

    def _eval_body(self, env: Environment) -> None:
        self.variable_sample.eval(env=env)


_STRING_RETURN_FUNCTIONS: Set[str] = {
    typ.name()
    for typ in [
        StatesFunctionNameType.Format,
        StatesFunctionNameType.JsonToString,
        StatesFunctionNameType.ArrayGetItem,
        StatesFunctionNameType.Base64Decode,
        StatesFunctionNameType.Base64Encode,
        StatesFunctionNameType.Hash,
        StatesFunctionNameType.UUID,
    ]
}


class CausePathJsonPath(CauseConst):
    def _eval_body(self, env: Environment) -> None:
        current_output = env.stack[-1]
        cause = extract_json(self.value, current_output)
        env.stack.append(cause)


class CausePathContextObject(CauseConst):
    def _eval_body(self, env: Environment) -> None:
        value = extract_json(self.value, env.states.context_object.context_object_data)
        env.stack.append(copy.deepcopy(value))


class CausePathIntrinsicFunction(CauseConst):
    function: Function

    def __init__(self, value: str) -> None:
        super().__init__(value=value)
        self.function, _ = IntrinsicParser.parse(value)
        if self.function.name.name not in _STRING_RETURN_FUNCTIONS:
            raise ValueError(
                f"Unsupported Intrinsic Function for CausePath declaration: '{self.value}'."
            )

    def _eval_body(self, env: Environment) -> None:
        self.function.eval(env=env)
