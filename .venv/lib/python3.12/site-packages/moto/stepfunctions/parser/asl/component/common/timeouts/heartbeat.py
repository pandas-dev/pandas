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
from moto.stepfunctions.parser.asl.component.common.string.string_expression import (
    StringJSONata,
    StringSampler,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.json_path import NoSuchJsonPathError


class Heartbeat(EvalComponent, abc.ABC):
    @abc.abstractmethod
    def _eval_seconds(self, env: Environment) -> int: ...

    def _eval_body(self, env: Environment) -> None:
        seconds = self._eval_seconds(env=env)
        env.stack.append(seconds)


class HeartbeatSeconds(Heartbeat):
    def __init__(self, heartbeat_seconds: int):
        if not isinstance(heartbeat_seconds, int) and heartbeat_seconds <= 0:
            raise ValueError(
                f"Expected non-negative integer for HeartbeatSeconds, got '{heartbeat_seconds}' instead."
            )
        self.heartbeat_seconds: Final[int] = heartbeat_seconds

    def _eval_seconds(self, env: Environment) -> int:
        return self.heartbeat_seconds


class HeartbeatSecondsJSONata(Heartbeat):
    string_jsonata: Final[StringJSONata]

    def __init__(self, string_jsonata: StringJSONata):
        super().__init__()
        self.string_jsonata = string_jsonata

    def _eval_seconds(self, env: Environment) -> int:
        self.string_jsonata.eval(env=env)
        # TODO: add snapshot tests to verify AWS's behaviour about non integer values.
        seconds = int(env.stack.pop())
        return seconds


class HeartbeatSecondsPath(Heartbeat):
    string_sampler: Final[StringSampler]

    def __init__(self, string_sampler: StringSampler):
        self.string_sampler = string_sampler

    def _eval_seconds(self, env: Environment) -> int:
        try:
            self.string_sampler.eval(env=env)
        except NoSuchJsonPathError as no_such_json_path_error:
            json_path = no_such_json_path_error.json_path
            cause = f"Invalid path '{json_path}' : No results for path: $['{json_path[2:]}']"
            raise FailureEventException(
                failure_event=FailureEvent(
                    env=env,
                    error_name=StatesErrorName(typ=StatesErrorNameType.StatesRuntime),
                    event_type=HistoryEventType.ExecutionFailed,
                    event_details=EventDetails(
                        executionFailedEventDetails=ExecutionFailedEventDetails(
                            error=StatesErrorNameType.StatesRuntime.to_name(),
                            cause=cause,
                        )
                    ),
                )
            )
        seconds = env.stack.pop()
        if not isinstance(seconds, int) and seconds <= 0:
            raise ValueError(
                f"Expected non-negative integer for HeartbeatSecondsPath, got '{seconds}' instead."
            )
        return seconds
