from typing import Optional

from moto.stepfunctions.parser.api import HistoryEventType, TaskFailedEventDetails
from moto.stepfunctions.parser.asl.component.common.error_name.custom_error_name import (
    CustomErrorName,
)
from moto.stepfunctions.parser.asl.component.common.error_name.failure_event import (
    FailureEvent,
    FailureEventException,
)
from moto.stepfunctions.parser.asl.component.state.fail.cause_decl import CauseDecl
from moto.stepfunctions.parser.asl.component.state.fail.error_decl import ErrorDecl
from moto.stepfunctions.parser.asl.component.state.state import CommonStateField
from moto.stepfunctions.parser.asl.component.state.state_props import StateProps
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.event.event_detail import EventDetails


class StateFail(CommonStateField):
    def __init__(self):
        super().__init__(
            state_entered_event_type=HistoryEventType.FailStateEntered,
            state_exited_event_type=None,
        )
        self.cause: Optional[CauseDecl] = None
        self.error: Optional[ErrorDecl] = None

    def from_state_props(self, state_props: StateProps) -> None:
        super().from_state_props(state_props)
        self.cause = state_props.get(CauseDecl)
        self.error = state_props.get(ErrorDecl)

    def _eval_state(self, env: Environment) -> None:
        task_failed_event_details = TaskFailedEventDetails()

        error_value = None
        if self.error:
            self.error.eval(env=env)
            error_value = env.stack.pop()
            task_failed_event_details["error"] = error_value

        if self.cause:
            self.cause.eval(env=env)
            cause_value = env.stack.pop()
            task_failed_event_details["cause"] = cause_value

        error_name = CustomErrorName(error_value) if error_value else None
        failure_event = FailureEvent(
            env=env,
            error_name=error_name,
            event_type=HistoryEventType.TaskFailed,
            event_details=EventDetails(
                taskFailedEventDetails=task_failed_event_details
            ),
        )
        raise FailureEventException(failure_event=failure_event)
