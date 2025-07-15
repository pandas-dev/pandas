import abc
from typing import Final, Tuple

from moto.stepfunctions.parser.api import (
    EvaluationFailedEventDetails,
    HistoryEventType,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
    FailureEventException,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name import (
    StatesErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.states_error_name_type import (
    StatesErrorNameType,
)
from moto.stepfunctions.parser.asl.component.common.jsonata.jsonata_template_value_array import (
    JSONataTemplateValueArray,
)
from moto.stepfunctions.parser.asl.component.common.jsonata.jsonata_template_value_terminal import (
    JSONataTemplateValueTerminalExpression,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str


class Items(EvalComponent, abc.ABC): ...


class ItemsArray(Items):
    jsonata_template_value_array: JSONataTemplateValueArray

    def __init__(self, jsonata_template_value_array: JSONataTemplateValueArray):
        super().__init__()
        self.jsonata_template_value_array = jsonata_template_value_array

    def _eval_body(self, env: Environment) -> None:
        self.jsonata_template_value_array.eval(env=env)


class ItemsJSONata(Items):
    jsonata_template_value_terminal_expression: Final[
        JSONataTemplateValueTerminalExpression
    ]

    def __init__(
        self,
        jsonata_template_value_terminal_expression: JSONataTemplateValueTerminalExpression,
    ):
        self.jsonata_template_value_terminal_expression = (
            jsonata_template_value_terminal_expression
        )

    def _eval_body(self, env: Environment) -> None:
        self.jsonata_template_value_terminal_expression.eval(env=env)
        items = env.stack[-1]
        if not isinstance(items, list):
            # FIXME: If we pass in a 'function' type, the JSONata lib will return a dict and the
            # 'unsupported result type state' wont be reached.
            def _get_jsonata_value_type_pair(items) -> Tuple[str, str]:
                if items is None:
                    return "null", "null"
                elif isinstance(items, (int, float)):
                    if isinstance(items, bool):
                        return "true" if items else "false", "boolean"
                    return items, "number"
                elif isinstance(items, str):
                    return f'"{items}"', "string"
                elif isinstance(items, dict):
                    return to_json_str(items, separators=(",", ":")), "object"

            expr = self.jsonata_template_value_terminal_expression.expression
            if jsonata_pair := _get_jsonata_value_type_pair(items):
                jsonata_value, jsonata_type = jsonata_pair
                error_cause = (
                    f"The JSONata expression '{expr}' specified for the field 'Items' returned an unexpected result type. "
                    f"Expected 'array', but was '{jsonata_type}' for value: {jsonata_value}"
                )
            else:
                error_cause = f"The JSONata expression '{expr}' for the field 'Items' returned an unsupported result type."

            error_name = StatesErrorName(
                typ=StatesErrorNameType.StatesQueryEvaluationError
            )
            failure_event = FailureEvent(
                env=env,
                error_name=error_name,
                event_type=HistoryEventType.EvaluationFailed,
                event_details=EventDetails(
                    evaluationFailedEventDetails=EvaluationFailedEventDetails(
                        error=error_name.error_name, cause=error_cause, location="Items"
                    )
                ),
            )
            raise FailureEventException(failure_event=failure_event)
