from __future__ import annotations

import datetime
import logging
from typing import Optional

from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.stepfunctions.models import Execution as SimpleExecution
from moto.stepfunctions.models import StateMachine
from moto.stepfunctions.parser.api import (
    CloudWatchEventsExecutionDataDetails,
    DescribeStateMachineForExecutionOutput,
    ExecutionListItem,
    ExecutionStatus,
    InvalidName,
    TraceHeader,
)
from moto.stepfunctions.parser.asl.eval.aws_execution_details import AWSExecutionDetails
from moto.stepfunctions.parser.asl.eval.contextobject.contex_object import (
    ContextObjectInitData,
)
from moto.stepfunctions.parser.asl.eval.contextobject.contex_object import (
    Execution as ContextObjectExecution,
)
from moto.stepfunctions.parser.asl.eval.contextobject.contex_object import (
    StateMachine as ContextObjectStateMachine,
)
from moto.stepfunctions.parser.asl.eval.program_state import (
    ProgramEnded,
    ProgramError,
    ProgramState,
    ProgramStopped,
    ProgramTimedOut,
)
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str
from moto.stepfunctions.parser.backend.execution_worker import ExecutionWorker
from moto.stepfunctions.parser.backend.execution_worker_comm import ExecutionWorkerComm
from moto.stepfunctions.parser.backend.state_machine import (
    StateMachineInstance,
    StateMachineVersion,
)

LOG = logging.getLogger(__name__)


class BaseExecutionWorkerComm(ExecutionWorkerComm):
    def __init__(self, execution: Execution):
        self.execution: Execution = execution

    def terminated(self) -> None:
        exit_program_state: ProgramState = (
            self.execution.exec_worker.env.program_state()
        )
        self.execution.stop_date = iso_8601_datetime_with_milliseconds()
        if isinstance(exit_program_state, ProgramEnded):
            self.execution.status = ExecutionStatus.SUCCEEDED
            self.execution.output = to_json_str(
                self.execution.exec_worker.env.inp, separators=(",", ":")
            )
        elif isinstance(exit_program_state, ProgramStopped):
            self.execution.status = ExecutionStatus.ABORTED
        elif isinstance(exit_program_state, ProgramError):
            self.execution.status = ExecutionStatus.FAILED
            self.execution.error = exit_program_state.error.get("error")
            self.execution.cause = exit_program_state.error.get("cause")
        elif isinstance(exit_program_state, ProgramTimedOut):
            self.execution.status = ExecutionStatus.TIMED_OUT
        else:
            raise RuntimeWarning(
                f"Execution ended with unsupported ProgramState type '{type(exit_program_state)}'."
            )


class Execution(SimpleExecution):
    def __init__(
        self,
        name: str,
        role_arn: str,
        account_id: str,
        region_name: str,
        state_machine: StateMachine,
        input_data: Optional[dict] = None,
        trace_header: Optional[TraceHeader] = None,
    ):
        super().__init__(
            region_name=region_name,
            account_id=account_id,
            state_machine_name=state_machine.name,
            state_machine_arn=state_machine.arn,
            execution_name=name,
            execution_input=input_data,
        )
        self.role_arn = role_arn
        self.state_machine = state_machine
        self.input_data = input_data
        self.input_details = CloudWatchEventsExecutionDataDetails(included=True)
        self.trace_header = trace_header
        self.status = None
        self.output_details = CloudWatchEventsExecutionDataDetails(included=True)
        self.exec_worker: Optional[ExecutionWorker] = None

    def to_describe_state_machine_for_execution_output(
        self,
    ) -> DescribeStateMachineForExecutionOutput:
        state_machine: StateMachineInstance = self.state_machine
        state_machine_arn = (
            state_machine.source_arn
            if isinstance(state_machine, StateMachineVersion)
            else state_machine.arn
        )
        out = DescribeStateMachineForExecutionOutput(
            stateMachineArn=state_machine_arn,
            name=state_machine.name,
            definition=state_machine.definition,
            roleArn=self.role_arn,
            # The date and time the state machine associated with an execution was updated.
            updateDate=state_machine.create_date,
        )
        revision_id = self.state_machine.revision_id
        if self.state_machine.revision_id:
            out["revisionId"] = revision_id
        return out

    def to_execution_list_item(self) -> ExecutionListItem:
        if isinstance(self.state_machine, StateMachineVersion):
            state_machine_arn = self.state_machine.source_arn
            state_machine_version_arn = self.state_machine.arn
        else:
            state_machine_arn = self.state_machine.arn
            state_machine_version_arn = None

        item = ExecutionListItem(
            executionArn=self.execution_arn,
            stateMachineArn=state_machine_arn,
            name=self.name,
            status=self.status,
            startDate=self.start_date,
            stopDate=self.stop_date,
        )
        if state_machine_version_arn is not None:
            item["stateMachineVersionArn"] = state_machine_version_arn
        return item

    def to_history_output(self):
        return self.exec_worker.env.event_history.get_event_history()

    @staticmethod
    def _to_serialized_date(timestamp: datetime.datetime) -> str:
        """See test in tests.aws.services.stepfunctions.v2.base.test_base.TestSnfBase.test_execution_dateformat"""
        return f'{timestamp.astimezone(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]}Z'

    def start(self) -> None:
        # TODO: checks exec_worker does not exists already?
        if self.exec_worker:
            raise InvalidName()  # TODO.
        self.exec_worker = ExecutionWorker(
            definition=self.state_machine.definition,
            input_data=self.input_data,
            exec_comm=BaseExecutionWorkerComm(self),
            context_object_init=ContextObjectInitData(
                Execution=ContextObjectExecution(
                    Id=self.execution_arn,
                    Input=self.input_data,
                    Name=self.name,
                    RoleArn=self.role_arn,
                    StartTime=self.start_date,
                ),
                StateMachine=ContextObjectStateMachine(
                    Id=self.state_machine.arn,
                    Name=self.state_machine.name,
                ),
            ),
            aws_execution_details=AWSExecutionDetails(
                account=self.account_id, region=self.region_name, role_arn=self.role_arn
            ),
        )
        self.status = ExecutionStatus.RUNNING
        self.exec_worker.start()

    def stop(
        self, stop_date: datetime.datetime, error: Optional[str], cause: Optional[str]
    ):
        exec_worker: Optional[ExecutionWorker] = self.exec_worker
        if not exec_worker:
            raise RuntimeError("No running executions.")
        exec_worker.stop(stop_date=stop_date, cause=cause, error=error)
