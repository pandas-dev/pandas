from moto.stepfunctions.parser.api import HistoryEventType
from moto.stepfunctions.parser.asl.component.common.flow.end import End
from moto.stepfunctions.parser.asl.component.common.flow.next import Next
from moto.stepfunctions.parser.asl.component.state.state import CommonStateField
from moto.stepfunctions.parser.asl.component.state.state_continue_with import (
    ContinueWithSuccess,
)
from moto.stepfunctions.parser.asl.component.state.state_props import StateProps
from moto.stepfunctions.parser.asl.eval.environment import Environment


class StateSucceed(CommonStateField):
    def __init__(self):
        super().__init__(
            state_entered_event_type=HistoryEventType.SucceedStateEntered,
            state_exited_event_type=HistoryEventType.SucceedStateExited,
        )

    def from_state_props(self, state_props: StateProps) -> None:
        super(StateSucceed, self).from_state_props(state_props)
        # TODO: assert all other fields are undefined?

        # No Next or End field: Succeed states are terminal states.
        if state_props.get(Next) or state_props.get(End):
            raise ValueError(
                f"No Next or End field: Succeed states are terminal states: with state '{self}'."
            )
        self.continue_with = ContinueWithSuccess()

    def _eval_state(self, env: Environment) -> None:
        pass
