import abc
from typing import Final

from moto.stepfunctions.parser.api import ExecutionFailedEventDetails, HistoryEventType
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
from moto.stepfunctions.parser.asl.component.common.jsonata.jsonata_template_value_terminal import (
    JSONataTemplateValueTerminalExpression,
)
from moto.stepfunctions.parser.asl.component.common.variable_sample import (
    VariableSample,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str
from moto.stepfunctions.parser.asl.utils.json_path import extract_json

DEFAULT_MAX_CONCURRENCY_VALUE: int = 0  # No limit.


class MaxConcurrencyDecl(EvalComponent, abc.ABC):
    @abc.abstractmethod
    def _eval_max_concurrency(self, env: Environment) -> int: ...

    def _eval_body(self, env: Environment) -> None:
        max_concurrency_value = self._eval_max_concurrency(env=env)
        env.stack.append(max_concurrency_value)


class MaxConcurrency(MaxConcurrencyDecl):
    max_concurrency_value: int

    def __init__(self, num: int = DEFAULT_MAX_CONCURRENCY_VALUE):
        super().__init__()
        self.max_concurrency_value = num

    def _eval_max_concurrency(self, env: Environment) -> int:
        return self.max_concurrency_value


class MaxConcurrencyJSONata(MaxConcurrencyDecl):
    jsonata_template_value_terminal_expression: Final[
        JSONataTemplateValueTerminalExpression
    ]

    def __init__(
        self,
        jsonata_template_value_terminal_expression: JSONataTemplateValueTerminalExpression,
    ):
        super().__init__()
        self.jsonata_template_value_terminal_expression = (
            jsonata_template_value_terminal_expression
        )

    def _eval_max_concurrency(self, env: Environment) -> int:
        self.jsonata_template_value_terminal_expression.eval(env=env)
        # TODO: add snapshot tests to verify AWS's behaviour about non integer values.
        seconds = int(env.stack.pop())
        return seconds


class MaxConcurrencyPathVar(MaxConcurrency):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        super().__init__()
        self.variable_sample = variable_sample

    def _eval_max_concurrency(self, env: Environment) -> int:
        self.variable_sample.eval(env=env)
        # TODO: add snapshot tests to verify AWS's behaviour about non integer values.
        max_concurrency: int = int(env.stack.pop())
        return max_concurrency


class MaxConcurrencyPath(MaxConcurrency):
    max_concurrency_path: str

    def __init__(self, max_concurrency_path: str):
        super().__init__()
        self.max_concurrency_path = max_concurrency_path

    def _eval_max_concurrency(self, env: Environment) -> int:
        inp = env.stack[-1]
        max_concurrency_value = extract_json(self.max_concurrency_path, inp)

        error_cause = None
        if not isinstance(max_concurrency_value, int):
            value_str = (
                to_json_str(max_concurrency_value)
                if not isinstance(max_concurrency_value, str)
                else max_concurrency_value
            )
            error_cause = f'The MaxConcurrencyPath field refers to value "{value_str}" which is not a valid integer: {self.max_concurrency_path}'
        elif max_concurrency_value < 0:
            error_cause = f"Expected non-negative integer for MaxConcurrency, got '{max_concurrency_value}' instead."

        if error_cause is not None:
            raise FailureEventException(
                failure_event=FailureEvent(
                    env=env,
                    error_name=StatesErrorName(typ=StatesErrorNameType.StatesRuntime),
                    event_type=HistoryEventType.ExecutionFailed,
                    event_details=EventDetails(
                        executionFailedEventDetails=ExecutionFailedEventDetails(
                            error=StatesErrorNameType.StatesRuntime.to_name(),
                            cause=error_cause,
                        )
                    ),
                )
            )

        return max_concurrency_value
