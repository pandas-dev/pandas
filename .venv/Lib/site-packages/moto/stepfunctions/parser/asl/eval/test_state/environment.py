from __future__ import annotations

from typing import Dict, Optional

from moto.stepfunctions.parser.api import Arn, InspectionData, StateMachineType
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.eval.evaluation_details import AWSExecutionDetails
from moto.stepfunctions.parser.asl.eval.event.event_manager import (
    EventHistoryContext,
)
from moto.stepfunctions.parser.asl.eval.event.logging import (
    CloudWatchLoggingSession,
)
from moto.stepfunctions.parser.asl.eval.program_state import (
    ProgramRunning,
)
from moto.stepfunctions.parser.asl.eval.states import ContextObjectData
from moto.stepfunctions.parser.asl.eval.test_state.program_state import (
    ProgramChoiceSelected,
)
from moto.stepfunctions.parser.asl.eval.variable_store import VariableStore
from moto.stepfunctions.parser.backend.activity import Activity


class TestStateEnvironment(Environment):
    inspection_data: InspectionData

    def __init__(
        self,
        aws_execution_details: AWSExecutionDetails,
        execution_type: StateMachineType,
        context: ContextObjectData,
        event_history_context: EventHistoryContext,
        activity_store: Dict[Arn, Activity],
        cloud_watch_logging_session: Optional[CloudWatchLoggingSession] = None,
    ):
        super().__init__(
            aws_execution_details=aws_execution_details,
            execution_type=execution_type,
            context=context,
            event_history_context=event_history_context,
            cloud_watch_logging_session=cloud_watch_logging_session,
            activity_store=activity_store,
        )
        self.inspection_data = InspectionData()

    def as_frame_of(  # pylint: disable=arguments-renamed
        cls,
        env: TestStateEnvironment,
        event_history_frame_cache: Optional[EventHistoryContext] = None,
    ) -> Environment:
        frame = super().as_frame_of(
            env=env, event_history_frame_cache=event_history_frame_cache
        )
        frame.inspection_data = env.inspection_data
        return frame

    def as_inner_frame_of(  # pylint: disable=arguments-renamed
        cls,
        env: TestStateEnvironment,
        variable_store: VariableStore,
        event_history_frame_cache: Optional[EventHistoryContext] = None,
    ) -> Environment:
        frame = super().as_inner_frame_of(
            env=env,
            event_history_frame_cache=event_history_frame_cache,
            variable_store=variable_store,
        )
        frame.inspection_data = env.inspection_data
        return frame

    def set_choice_selected(self, next_state_name: str) -> None:
        with self._state_mutex:
            if isinstance(self._program_state, ProgramRunning):
                self._program_state = ProgramChoiceSelected(
                    next_state_name=next_state_name
                )
                self.program_state_event.set()
                self.program_state_event.clear()
            else:
                raise RuntimeError(
                    "Cannot set choice selected for non running ProgramState."
                )
