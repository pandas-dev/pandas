from __future__ import annotations

import datetime
import json
import logging
from typing import Optional

from moto.stepfunctions.parser.api import (
    Arn,
    CloudWatchEventsExecutionDataDetails,
    DescribeExecutionOutput,
    DescribeStateMachineForExecutionOutput,
    ExecutionListItem,
    ExecutionStatus,
    GetExecutionHistoryOutput,
    HistoryEventList,
    InvalidName,
    SensitiveCause,
    SensitiveError,
    StartExecutionOutput,
    StartSyncExecutionOutput,
    StateMachineType,
    SyncExecutionStatus,
    Timestamp,
    TraceHeader,
    VariableReferences,
)
from moto.stepfunctions.parser.asl.eval.evaluation_details import (
    AWSExecutionDetails,
    EvaluationDetails,
    ExecutionDetails,
    StateMachineDetails,
)
from moto.stepfunctions.parser.asl.eval.event.logging import (
    CloudWatchLoggingSession,
)
from moto.stepfunctions.parser.asl.eval.program_state import (
    ProgramEnded,
    ProgramError,
    ProgramState,
    ProgramStopped,
    ProgramTimedOut,
)
from moto.stepfunctions.parser.asl.static_analyser.variable_references_static_analyser import (
    VariableReferencesStaticAnalyser,
)
from moto.stepfunctions.parser.asl.utils.encoding import to_json_str
from moto.stepfunctions.parser.backend.activity import Activity
from moto.stepfunctions.parser.backend.execution_worker import (
    ExecutionWorker,
    SyncExecutionWorker,
)
from moto.stepfunctions.parser.backend.execution_worker_comm import (
    ExecutionWorkerCommunication,
)
from moto.stepfunctions.parser.backend.state_machine import (
    StateMachineInstance,
    StateMachineVersion,
)

LOG = logging.getLogger(__name__)


class BaseExecutionWorkerCommunication(ExecutionWorkerCommunication):
    execution: Execution

    def __init__(self, execution: Execution):
        self.execution = execution

    def _reflect_execution_status(self):
        exit_program_state: ProgramState = (
            self.execution.exec_worker.env.program_state()
        )
        self.execution.stop_date = datetime.datetime.now(tz=datetime.timezone.utc)
        if isinstance(exit_program_state, ProgramEnded):
            self.execution.exec_status = ExecutionStatus.SUCCEEDED
            self.execution.output = self.execution.exec_worker.env.states.get_input()
        elif isinstance(exit_program_state, ProgramStopped):
            self.execution.exec_status = ExecutionStatus.ABORTED
        elif isinstance(exit_program_state, ProgramError):
            self.execution.exec_status = ExecutionStatus.FAILED
            self.execution.error = exit_program_state.error.get("error")
            self.execution.cause = exit_program_state.error.get("cause")
        elif isinstance(exit_program_state, ProgramTimedOut):
            self.execution.exec_status = ExecutionStatus.TIMED_OUT
        else:
            raise RuntimeWarning(
                f"Execution ended with unsupported ProgramState type '{type(exit_program_state)}'."
            )

    def terminated(self) -> None:
        self._reflect_execution_status()


class Execution:
    name: str
    sm_type: StateMachineType
    role_arn: Arn
    exec_arn: Arn

    account_id: str
    region_name: str

    state_machine: StateMachineInstance
    start_date: Timestamp
    input_data: Optional[json]
    input_details: Optional[CloudWatchEventsExecutionDataDetails]
    trace_header: Optional[TraceHeader]
    _cloud_watch_logging_session: Optional[CloudWatchLoggingSession]

    exec_status: Optional[ExecutionStatus]
    stop_date: Optional[Timestamp]

    output: Optional[json]
    output_details: Optional[CloudWatchEventsExecutionDataDetails]

    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]

    exec_worker: Optional[ExecutionWorker]

    _activity_store: dict[Arn, Activity]

    def __init__(
        self,
        name: str,
        sm_type: StateMachineType,
        role_arn: Arn,
        exec_arn: Arn,
        account_id: str,
        region_name: str,
        state_machine: StateMachineInstance,
        start_date: Timestamp,
        cloud_watch_logging_session: Optional[CloudWatchLoggingSession],
        activity_store: dict[Arn, Activity],
        input_data: str,
        trace_header: Optional[TraceHeader] = None,
    ):
        self.name = name
        self.sm_type = sm_type
        self.role_arn = role_arn
        self.exec_arn = exec_arn
        self.execution_arn = exec_arn
        self.account_id = account_id
        self.region_name = region_name
        self.state_machine = state_machine
        self._cloud_watch_logging_session = cloud_watch_logging_session
        self.input_data = json.loads(input_data)
        self.input_details = CloudWatchEventsExecutionDataDetails(included=True)
        self.trace_header = trace_header
        self.exec_status = None
        self.stop_date = None
        self.output = None
        self.output_details = CloudWatchEventsExecutionDataDetails(included=True)
        self.exec_worker = None
        self.error = None
        self.cause = None
        self._activity_store = activity_store

        # Compatibility with mock SFN
        self.state_machine_arn = state_machine.arn
        self.start_date = start_date
        self.execution_input = input_data

    @property
    def status(self):
        return self.exec_status.value

    def to_start_output(self) -> StartExecutionOutput:
        return StartExecutionOutput(
            executionArn=self.exec_arn, startDate=self.start_date
        )

    def to_describe_output(self) -> DescribeExecutionOutput:
        describe_output = DescribeExecutionOutput(
            executionArn=self.exec_arn,
            stateMachineArn=self.state_machine.arn,
            name=self.name,
            status=self.exec_status,
            startDate=self.start_date,
            stopDate=self.stop_date,
            input=to_json_str(self.input_data, separators=(",", ":")),
            inputDetails=self.input_details,
            traceHeader=self.trace_header,
        )
        if describe_output["status"] == ExecutionStatus.SUCCEEDED:
            describe_output["output"] = to_json_str(self.output, separators=(",", ":"))
            describe_output["outputDetails"] = self.output_details
        if self.error is not None:
            describe_output["error"] = self.error
        if self.cause is not None:
            describe_output["cause"] = self.cause
        return describe_output

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
            loggingConfiguration=state_machine.logging_config,
        )
        revision_id = self.state_machine.revision_id
        if self.state_machine.revision_id:
            out["revisionId"] = revision_id
        variable_references: VariableReferences = (
            VariableReferencesStaticAnalyser.process_and_get(
                definition=self.state_machine.definition
            )
        )
        if variable_references:
            out["variableReferences"] = variable_references
        return out

    def to_execution_list_item(self) -> ExecutionListItem:
        if isinstance(self.state_machine, StateMachineVersion):
            state_machine_arn = self.state_machine.source_arn
            state_machine_version_arn = self.state_machine.arn
        else:
            state_machine_arn = self.state_machine.arn
            state_machine_version_arn = None

        item = ExecutionListItem(
            executionArn=self.exec_arn,
            stateMachineArn=state_machine_arn,
            name=self.name,
            status=self.exec_status,
            startDate=self.start_date,
            stopDate=self.stop_date,
        )
        if state_machine_version_arn is not None:
            item["stateMachineVersionArn"] = state_machine_version_arn
        return item

    def to_history_output(self) -> GetExecutionHistoryOutput:
        env = self.exec_worker.env
        event_history: HistoryEventList = []
        if env is not None:
            # The execution has not started yet.
            event_history: HistoryEventList = env.event_manager.get_event_history()
        return GetExecutionHistoryOutput(events=event_history)

    def _get_start_execution_worker_comm(self) -> BaseExecutionWorkerCommunication:
        return BaseExecutionWorkerCommunication(self)

    def _get_start_aws_execution_details(self) -> AWSExecutionDetails:
        return AWSExecutionDetails(
            account=self.account_id, region=self.region_name, role_arn=self.role_arn
        )

    def get_start_execution_details(self) -> ExecutionDetails:
        return ExecutionDetails(
            arn=self.exec_arn,
            name=self.name,
            role_arn=self.role_arn,
            inpt=self.input_data,
            start_time=self.start_date,
        )

    def get_start_state_machine_details(self) -> StateMachineDetails:
        return StateMachineDetails(
            arn=self.state_machine.arn,
            name=self.state_machine.name,
            typ=self.state_machine.sm_type,
            definition=self.state_machine.definition,
        )

    def _get_start_execution_worker(self) -> ExecutionWorker:
        return ExecutionWorker(
            evaluation_details=EvaluationDetails(
                aws_execution_details=self._get_start_aws_execution_details(),
                execution_details=self.get_start_execution_details(),
                state_machine_details=self.get_start_state_machine_details(),
            ),
            exec_comm=self._get_start_execution_worker_comm(),
            cloud_watch_logging_session=self._cloud_watch_logging_session,
            activity_store=self._activity_store,
        )

    def start(self) -> None:
        # TODO: checks exec_worker does not exists already?
        if self.exec_worker:
            raise InvalidName()  # TODO.
        self.exec_worker = self._get_start_execution_worker()
        self.exec_status = ExecutionStatus.RUNNING
        self.exec_worker.start()

    def stop(
        self, stop_date: datetime.datetime, error: Optional[str], cause: Optional[str]
    ):
        exec_worker: Optional[ExecutionWorker] = self.exec_worker
        if exec_worker:
            exec_worker.stop(stop_date=stop_date, cause=cause, error=error)


class SyncExecutionWorkerCommunication(BaseExecutionWorkerCommunication):
    execution: SyncExecution

    def _reflect_execution_status(self) -> None:
        super()._reflect_execution_status()
        exit_status: ExecutionStatus = self.execution.exec_status
        if exit_status == ExecutionStatus.SUCCEEDED:
            self.execution.sync_execution_status = SyncExecutionStatus.SUCCEEDED
        elif exit_status == ExecutionStatus.TIMED_OUT:
            self.execution.sync_execution_status = SyncExecutionStatus.TIMED_OUT
        else:
            self.execution.sync_execution_status = SyncExecutionStatus.FAILED


class SyncExecution(Execution):
    sync_execution_status: Optional[SyncExecutionStatus] = None

    def _get_start_execution_worker(self) -> SyncExecutionWorker:
        return SyncExecutionWorker(
            evaluation_details=EvaluationDetails(
                aws_execution_details=self._get_start_aws_execution_details(),
                execution_details=self.get_start_execution_details(),
                state_machine_details=self.get_start_state_machine_details(),
            ),
            exec_comm=self._get_start_execution_worker_comm(),
            cloud_watch_logging_session=self._cloud_watch_logging_session,
            activity_store=self._activity_store,
        )

    def _get_start_execution_worker_comm(self) -> BaseExecutionWorkerCommunication:
        return SyncExecutionWorkerCommunication(self)

    def to_start_sync_execution_output(self) -> StartSyncExecutionOutput:
        start_output = StartSyncExecutionOutput(
            executionArn=self.exec_arn,
            stateMachineArn=self.state_machine.arn,
            name=self.name,
            status=self.sync_execution_status,
            startDate=self.start_date,
            stopDate=self.stop_date,
            input=to_json_str(self.input_data, separators=(",", ":")),
            inputDetails=self.input_details,
            traceHeader=self.trace_header,
        )
        if self.sync_execution_status == SyncExecutionStatus.SUCCEEDED:
            start_output["output"] = to_json_str(self.output, separators=(",", ":"))
        if self.output_details:
            start_output["outputDetails"] = self.output_details
        if self.error is not None:
            start_output["error"] = self.error
        if self.cause is not None:
            start_output["cause"] = self.cause
        return start_output
