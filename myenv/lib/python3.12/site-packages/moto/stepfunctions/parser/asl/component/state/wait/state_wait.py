from __future__ import annotations

from moto.stepfunctions.parser.api import HistoryEventType
from moto.stepfunctions.parser.asl.component.state.state import CommonStateField
from moto.stepfunctions.parser.asl.component.state.state_props import StateProps
from moto.stepfunctions.parser.asl.component.state.wait.wait_function.wait_function import (
    WaitFunction,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class StateWait(CommonStateField):
    wait_function: WaitFunction

    def __init__(self):
        super().__init__(
            state_entered_event_type=HistoryEventType.WaitStateEntered,
            state_exited_event_type=HistoryEventType.WaitStateExited,
        )

    def from_state_props(self, state_props: StateProps) -> None:
        super(StateWait, self).from_state_props(state_props)
        self.wait_function = state_props.get(
            typ=WaitFunction,
            raise_on_missing=ValueError(
                f"Undefined WaitFunction for StateWait: '{self}'."
            ),
        )

    def _eval_state(self, env: Environment) -> None:
        self.wait_function.eval(env)
