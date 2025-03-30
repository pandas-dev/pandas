import copy
import datetime
import json
from typing import Any, Dict, List, Optional

from moto.core.common_models import BackendDict
from moto.stepfunctions.models import StateMachine, StepFunctionBackend
from moto.stepfunctions.parser.api import (
    Definition,
    EncryptionConfiguration,
    ExecutionStatus,
    GetExecutionHistoryOutput,
    InvalidDefinition,
    InvalidExecutionInput,
    InvalidToken,
    LoggingConfiguration,
    MissingRequiredParameter,
    Name,
    ResourceNotFound,
    SendTaskFailureOutput,
    SendTaskHeartbeatOutput,
    SendTaskSuccessOutput,
    SensitiveCause,
    SensitiveData,
    SensitiveError,
    TaskDoesNotExist,
    TaskTimedOut,
    TaskToken,
    TraceHeader,
    TracingConfiguration,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.iteration.itemprocessor.map_run_record import (
    MapRunRecord,
)
from moto.stepfunctions.parser.asl.eval.callback.callback import (
    CallbackConsumerTimeout,
    CallbackNotifyConsumerError,
    CallbackOutcomeFailure,
    CallbackOutcomeSuccess,
)
from moto.stepfunctions.parser.asl.parse.asl_parser import (
    AmazonStateLanguageParser,
    ASLParserException,
)
from moto.stepfunctions.parser.backend.execution import Execution


class StepFunctionsParserBackend(StepFunctionBackend):
    def _get_executions(self, execution_status: Optional[ExecutionStatus] = None):
        executions = []
        for sm in self.state_machines:
            for execution in sm.executions:
                if execution_status is None or execution_status == execution.status:
                    executions.append(execution)
        return executions

    def _revision_by_name(self, name: str) -> Optional[StateMachine]:
        for state_machine in self.state_machines:
            if state_machine.name == name:
                return state_machine
        return None

    @staticmethod
    def _validate_definition(definition: str):
        # Validate
        # TODO: pass through static analyser.
        try:
            AmazonStateLanguageParser.parse(definition)
        except ASLParserException as asl_parser_exception:
            raise InvalidDefinition(message=repr(asl_parser_exception))
        except Exception as exception:
            exception_name = exception.__class__.__name__
            exception_args = list(exception.args)
            raise InvalidDefinition(
                message=f"Error={exception_name} Args={exception_args} in definition '{definition}'."
            )

    def create_state_machine(
        self,
        name: str,
        definition: str,
        roleArn: str,
        tags: Optional[List[Dict[str, str]]] = None,
        publish: Optional[bool] = None,
        loggingConfiguration: Optional[LoggingConfiguration] = None,
        tracingConfiguration: Optional[TracingConfiguration] = None,
        encryptionConfiguration: Optional[EncryptionConfiguration] = None,
        version_description: Optional[str] = None,
    ) -> StateMachine:
        StepFunctionsParserBackend._validate_definition(definition=definition)

        return super().create_state_machine(
            name=name,
            definition=definition,
            roleArn=roleArn,
            tags=tags,
            publish=publish,
            loggingConfiguration=loggingConfiguration,
            tracingConfiguration=tracingConfiguration,
            encryptionConfiguration=encryptionConfiguration,
            version_description=version_description,
        )

    def send_task_heartbeat(self, task_token: TaskToken) -> SendTaskHeartbeatOutput:
        running_executions = self._get_executions(ExecutionStatus.RUNNING)
        for execution in running_executions:
            try:
                if execution.exec_worker.env.callback_pool_manager.heartbeat(
                    callback_id=task_token
                ):
                    return
            except CallbackNotifyConsumerError as consumer_error:
                if isinstance(consumer_error, CallbackConsumerTimeout):
                    raise TaskTimedOut()
                else:
                    raise TaskDoesNotExist()
        raise InvalidToken()

    def send_task_success(
        self, task_token: TaskToken, outcome: str
    ) -> SendTaskSuccessOutput:
        outcome = CallbackOutcomeSuccess(callback_id=task_token, output=outcome)
        running_executions = self._get_executions(ExecutionStatus.RUNNING)
        for execution in running_executions:
            try:
                if execution.exec_worker.env.callback_pool_manager.notify(
                    callback_id=task_token, outcome=outcome
                ):
                    return
            except CallbackNotifyConsumerError as consumer_error:
                if isinstance(consumer_error, CallbackConsumerTimeout):
                    raise TaskTimedOut()
                else:
                    raise TaskDoesNotExist()
        raise InvalidToken()

    def send_task_failure(
        self,
        task_token: TaskToken,
        error: SensitiveError = None,
        cause: SensitiveCause = None,
    ) -> SendTaskFailureOutput:
        outcome = CallbackOutcomeFailure(
            callback_id=task_token, error=error, cause=cause
        )
        for execution in self._get_executions():
            try:
                if execution.exec_worker.env.callback_pool_manager.notify(
                    callback_id=task_token, outcome=outcome
                ):
                    return SendTaskFailureOutput()
            except CallbackNotifyConsumerError as consumer_error:
                if isinstance(consumer_error, CallbackConsumerTimeout):
                    raise TaskTimedOut()
                else:
                    raise TaskDoesNotExist()
        raise InvalidToken()

    def start_execution(
        self,
        state_machine_arn: str,
        name: Name = None,
        execution_input: SensitiveData = None,
        trace_header: TraceHeader = None,
    ) -> Execution:
        state_machine = self.describe_state_machine(state_machine_arn)

        # Update event change parameters about the state machine and should not affect those about this execution.
        state_machine_clone = copy.deepcopy(state_machine)

        if execution_input is None:
            input_data = dict()
        else:
            try:
                input_data = json.loads(execution_input)
            except Exception as ex:
                raise InvalidExecutionInput(
                    str(ex)
                )  # TODO: report parsing error like AWS.

        exec_name = name  # TODO: validate name format

        execution_arn = "arn:{}:states:{}:{}:execution:{}:{}"
        execution_arn = execution_arn.format(
            self.partition,
            self.region_name,
            self.account_id,
            state_machine.name,
            name,
        )

        execution = Execution(
            name=exec_name,
            sm_type=state_machine_clone.sm_type,
            role_arn=state_machine_clone.roleArn,
            exec_arn=execution_arn,
            account_id=self.account_id,
            region_name=self.region_name,
            state_machine=state_machine_clone,
            start_date=datetime.datetime.now(tz=datetime.timezone.utc),
            cloud_watch_logging_session=None,
            input_data=input_data,
            trace_header=trace_header,
            activity_store={},
        )
        state_machine.executions.append(execution)

        execution.start()
        return execution

    def update_state_machine(
        self,
        arn: str,
        definition: Definition = None,
        role_arn: str = None,
        logging_configuration: LoggingConfiguration = None,
        tracing_configuration: TracingConfiguration = None,
        encryption_configuration: EncryptionConfiguration = None,
        publish: Optional[bool] = None,
        version_description: str = None,
    ) -> StateMachine:
        if not any(
            [
                definition,
                role_arn,
                logging_configuration,
                tracing_configuration,
                encryption_configuration,
            ]
        ):
            raise MissingRequiredParameter(
                "Either the definition, the role ARN, the LoggingConfiguration, the EncryptionConfiguration or the TracingConfiguration must be specified"
            )

        if definition is not None:
            self._validate_definition(definition=definition)

        return super().update_state_machine(
            arn,
            definition,
            role_arn,
            logging_configuration=logging_configuration,
            tracing_configuration=tracing_configuration,
            encryption_configuration=encryption_configuration,
            publish=publish,
            version_description=version_description,
        )

    def describe_map_run(self, map_run_arn: str) -> Dict[str, Any]:
        for execution in self._get_executions():
            map_run_record: Optional[MapRunRecord] = (
                execution.exec_worker.env.map_run_record_pool_manager.get(map_run_arn)
            )
            if map_run_record is not None:
                return map_run_record.describe()
        raise ResourceNotFound()

    def list_map_runs(self, execution_arn: str) -> Dict[str, Any]:
        """
        Pagination is not yet implemented
        """
        execution = self.describe_execution(execution_arn=execution_arn)
        map_run_records: List[MapRunRecord] = (
            execution.exec_worker.env.map_run_record_pool_manager.get_all()
        )
        return dict(
            mapRuns=[map_run_record.list_item() for map_run_record in map_run_records]
        )

    def update_map_run(
        self,
        map_run_arn: str,
        max_concurrency: int,
        tolerated_failure_count: str,
        tolerated_failure_percentage: str,
    ) -> None:
        # TODO: investigate behaviour of empty requests.
        for execution in self._get_executions():
            map_run_record = execution.exec_worker.env.map_run_record_pool_manager.get(
                map_run_arn
            )
            if map_run_record is not None:
                map_run_record.update(
                    max_concurrency=max_concurrency,
                    tolerated_failure_count=tolerated_failure_count,
                    tolerated_failure_percentage=tolerated_failure_percentage,
                )
                return
        raise ResourceNotFound()

    def get_execution_history(self, execution_arn: str) -> GetExecutionHistoryOutput:
        execution = self.describe_execution(execution_arn=execution_arn)
        return execution.to_history_output()


stepfunctions_parser_backends = BackendDict(StepFunctionsParserBackend, "stepfunctions")
