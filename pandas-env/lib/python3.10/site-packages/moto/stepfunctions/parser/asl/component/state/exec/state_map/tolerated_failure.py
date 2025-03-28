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

TOLERATED_FAILURE_COUNT_MIN: int = 0
TOLERATED_FAILURE_COUNT_DEFAULT: int = 0
TOLERATED_FAILURE_PERCENTAGE_MIN: float = 0.0
TOLERATED_FAILURE_PERCENTAGE_DEFAULT: float = 0.0
TOLERATED_FAILURE_PERCENTAGE_MAX: float = 100.0


class ToleratedFailureCountDecl(EvalComponent, abc.ABC):
    @abc.abstractmethod
    def _eval_tolerated_failure_count(self, env: Environment) -> int: ...

    def _eval_body(self, env: Environment) -> None:
        tolerated_failure_count = self._eval_tolerated_failure_count(env=env)
        env.stack.append(tolerated_failure_count)


class ToleratedFailureCount(ToleratedFailureCountDecl):
    tolerated_failure_count: int

    def __init__(self, tolerated_failure_count: int = TOLERATED_FAILURE_COUNT_DEFAULT):
        self.tolerated_failure_count = tolerated_failure_count

    def _eval_tolerated_failure_count(self, env: Environment) -> int:
        return self.tolerated_failure_count


class ToleratedFailureCountJSONata(ToleratedFailureCountDecl):
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

    def _eval_tolerated_failure_count(self, env: Environment) -> int:
        # TODO: add snapshot tests to verify AWS's behaviour about non integer values.
        self.jsonata_template_value_terminal_expression.eval(env=env)
        failure_count: int = int(env.stack.pop())
        return failure_count


class ToleratedFailureCountPathVar(ToleratedFailureCountDecl):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        super().__init__()
        self.variable_sample = variable_sample

    def _eval_tolerated_failure_count(self, env: Environment) -> int:
        self.variable_sample.eval(env=env)
        # TODO: add snapshot tests to verify AWS's behaviour about non integer values.
        tolerated_failure_count: int = int(env.stack.pop())
        return tolerated_failure_count


class ToleratedFailureCountPath(ToleratedFailureCountDecl):
    tolerated_failure_count_path: str

    def __init__(self, tolerated_failure_count_path: str):
        self.tolerated_failure_count_path = tolerated_failure_count_path

    def _eval_tolerated_failure_count(self, env: Environment) -> int:
        inp = env.stack[-1]
        tolerated_failure_count = extract_json(self.tolerated_failure_count_path, inp)

        error_cause = None
        if not isinstance(tolerated_failure_count, int):
            value_str = (
                to_json_str(tolerated_failure_count)
                if not isinstance(tolerated_failure_count, str)
                else tolerated_failure_count
            )
            error_cause = (
                f'The ToleratedFailureCountPath field refers to value "{value_str}" '
                f"which is not a valid integer: {self.tolerated_failure_count_path}"
            )

        elif tolerated_failure_count < TOLERATED_FAILURE_COUNT_MIN:
            error_cause = "ToleratedFailureCount cannot be negative."

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

        return tolerated_failure_count


class ToleratedFailurePercentageDecl(EvalComponent, abc.ABC):
    @abc.abstractmethod
    def _eval_tolerated_failure_percentage(self, env: Environment) -> float: ...

    def _eval_body(self, env: Environment) -> None:
        tolerated_failure_percentage = self._eval_tolerated_failure_percentage(env=env)
        env.stack.append(tolerated_failure_percentage)


class ToleratedFailurePercentage(ToleratedFailurePercentageDecl):
    tolerated_failure_percentage: float

    def __init__(
        self, tolerated_failure_percentage: float = TOLERATED_FAILURE_PERCENTAGE_DEFAULT
    ):
        self.tolerated_failure_percentage = tolerated_failure_percentage

    def _eval_tolerated_failure_percentage(self, env: Environment) -> float:
        return self.tolerated_failure_percentage


class ToleratedFailurePercentageJSONata(ToleratedFailurePercentageDecl):
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

    def _eval_tolerated_failure_percentage(self, env: Environment) -> float:
        # TODO: add snapshot tests to verify AWS's behaviour about non floating values.
        self.jsonata_template_value_terminal_expression.eval(env=env)
        failure_percentage: int = int(env.stack.pop())
        return failure_percentage


class ToleratedFailurePercentagePathVar(ToleratedFailurePercentageDecl):
    variable_sample: VariableSample

    def __init__(self, variable_sample: VariableSample):
        super().__init__()
        self.variable_sample = variable_sample

    def _eval_tolerated_failure_percentage(self, env: Environment) -> float:
        self.variable_sample.eval(env=env)
        # TODO: add snapshot tests to verify AWS's behaviour about non floating values.
        tolerated_failure_percentage: float = float(env.stack.pop())
        return tolerated_failure_percentage


class ToleratedFailurePercentagePath(ToleratedFailurePercentageDecl):
    tolerate_failure_percentage_path: str

    def __init__(self, tolerate_failure_percentage_path: str):
        self.tolerate_failure_percentage_path = tolerate_failure_percentage_path

    def _eval_tolerated_failure_percentage(self, env: Environment) -> float:
        inp = env.stack[-1]
        tolerated_failure_percentage = extract_json(
            self.tolerate_failure_percentage_path, inp
        )

        if isinstance(tolerated_failure_percentage, int):
            tolerated_failure_percentage = float(tolerated_failure_percentage)

        error_cause = None
        if not isinstance(tolerated_failure_percentage, float):
            value_str = (
                to_json_str(tolerated_failure_percentage)
                if not isinstance(tolerated_failure_percentage, str)
                else tolerated_failure_percentage
            )
            error_cause = (
                f'The ToleratedFailurePercentagePath field refers to value "{value_str}" '
                f"which is not a valid float: {self.tolerate_failure_percentage_path}"
            )
        elif (
            not TOLERATED_FAILURE_PERCENTAGE_MIN
            <= tolerated_failure_percentage
            <= TOLERATED_FAILURE_PERCENTAGE_MAX
        ):
            error_cause = "ToleratedFailurePercentage must be between 0 and 100."

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

        return tolerated_failure_percentage
