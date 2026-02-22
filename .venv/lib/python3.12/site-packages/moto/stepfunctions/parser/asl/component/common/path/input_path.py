from typing import Final, Optional

from moto.stepfunctions.parser.api import HistoryEventType, TaskFailedEventDetails
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
    StringJsonPath,
    StringSampler,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails
from moto.stepfunctions.parser.asl.utils.json_path import NoSuchJsonPathError


class InputPath(EvalComponent):
    string_sampler: Final[Optional[StringSampler]]

    def __init__(self, string_sampler: Optional[StringSampler]):
        self.string_sampler = string_sampler

    def _eval_body(self, env: Environment) -> None:
        if self.string_sampler is None:
            env.stack.append({})
            return
        if isinstance(self.string_sampler, StringJsonPath):
            # JsonPaths are sampled from a given state, hence pass the state's input.
            env.stack.append(env.states.get_input())
        try:
            self.string_sampler.eval(env=env)
        except NoSuchJsonPathError as no_such_json_path_error:
            json_path = no_such_json_path_error.json_path
            cause = f"Invalid path '{json_path}' : No results for path: $['{json_path[2:]}']"
            raise FailureEventException(
                failure_event=FailureEvent(
                    env=env,
                    error_name=StatesErrorName(typ=StatesErrorNameType.StatesRuntime),
                    event_type=HistoryEventType.TaskFailed,
                    event_details=EventDetails(
                        taskFailedEventDetails=TaskFailedEventDetails(
                            error=StatesErrorNameType.StatesRuntime.to_name(),
                            cause=cause,
                        )
                    ),
                )
            )
