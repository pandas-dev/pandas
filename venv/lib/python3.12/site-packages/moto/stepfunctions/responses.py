import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.config import default_user_config
from moto.core.responses import BaseResponse
from moto.core.utils import iso_8601_datetime_with_milliseconds

from .exceptions import ValidationException
from .models import StepFunctionBackend, stepfunctions_backends
from .parser.api import ExecutionStatus


def get_backend(account: str, region: str) -> StepFunctionBackend:
    sfn_config = default_user_config.get("stepfunctions", {})
    if sfn_config.get("execute_state_machine", False):
        from .parser.models import stepfunctions_parser_backends

        return stepfunctions_parser_backends[account][region]
    else:
        return stepfunctions_backends[account][region]


class StepFunctionResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="stepfunctions")

    @property
    def stepfunction_backend(self) -> StepFunctionBackend:
        return get_backend(self.current_account, self.region)

    def create_state_machine(self) -> TYPE_RESPONSE:
        name = self._get_param("name")
        definition = self._get_param("definition")
        roleArn = self._get_param("roleArn")
        tags = self._get_param("tags")
        publish = self._get_param("publish")
        encryptionConfiguration = self._get_param("encryptionConfiguration")
        loggingConfiguration = self._get_param("loggingConfiguration")
        tracingConfiguration = self._get_param("tracingConfiguration")
        version_description = self._get_param("versionDescription")

        if version_description and not publish:
            raise ValidationException(
                "Version description can only be set when publish is true"
            )

        state_machine = self.stepfunction_backend.create_state_machine(
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
        response = {
            "creationDate": state_machine.creation_date,
            "stateMachineArn": state_machine.arn,
        }
        if state_machine.latest_version:
            response["stateMachineVersionArn"] = state_machine.latest_version.arn
        return 200, {}, json.dumps(response)

    def list_state_machines(self) -> TYPE_RESPONSE:
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")
        results, next_token = self.stepfunction_backend.list_state_machines(
            max_results=max_results, next_token=next_token
        )
        state_machines = [
            {
                "creationDate": sm.creation_date,
                "name": sm.name,
                "stateMachineArn": sm.arn,
            }
            for sm in results
        ]
        response = {"stateMachines": state_machines, "nextToken": next_token}
        return 200, {}, json.dumps(response)

    def describe_state_machine(self) -> TYPE_RESPONSE:
        arn = self._get_param("stateMachineArn")
        return self._describe_state_machine(arn)

    def _describe_state_machine(self, state_machine_arn: str) -> TYPE_RESPONSE:
        state_machine = self.stepfunction_backend.describe_state_machine(
            state_machine_arn
        )
        response = {
            "creationDate": state_machine.creation_date,
            "stateMachineArn": state_machine.arn,
            "definition": state_machine.definition,
            "name": state_machine.name,
            "roleArn": state_machine.roleArn,
            "status": "ACTIVE",
            "type": state_machine.type,
            "encryptionConfiguration": state_machine.encryptionConfiguration,
            "tracingConfiguration": state_machine.tracingConfiguration,
            "loggingConfiguration": state_machine.loggingConfiguration,
        }
        if state_machine.description:
            response["description"] = state_machine.description
        return 200, {}, json.dumps(response)

    def delete_state_machine(self) -> TYPE_RESPONSE:
        arn = self._get_param("stateMachineArn")
        self.stepfunction_backend.delete_state_machine(arn)
        return 200, {}, json.dumps("{}")

    def update_state_machine(self) -> TYPE_RESPONSE:
        arn = self._get_param("stateMachineArn")
        definition = self._get_param("definition")
        role_arn = self._get_param("roleArn")
        tracing_config = self._get_param("tracingConfiguration")
        encryption_config = self._get_param("encryptionConfiguration")
        logging_config = self._get_param("loggingConfiguration")
        publish = self._get_param("publish")
        version_description = self._get_param("versionDescription")

        if version_description and not publish:
            raise ValidationException(
                "Version description can only be set when publish is true"
            )

        state_machine = self.stepfunction_backend.update_state_machine(
            arn=arn,
            definition=definition,
            role_arn=role_arn,
            tracing_configuration=tracing_config,
            encryption_configuration=encryption_config,
            logging_configuration=logging_config,
            publish=publish,
            version_description=version_description,
        )
        response = {"updateDate": state_machine.update_date}
        if publish:
            response["stateMachineVersionArn"] = state_machine.latest_version.arn  # type: ignore
        return 200, {}, json.dumps(response)

    def list_tags_for_resource(self) -> TYPE_RESPONSE:
        arn = self._get_param("resourceArn")
        response = self.stepfunction_backend.list_tags_for_resource(arn)
        return 200, {}, json.dumps(response)

    def tag_resource(self) -> TYPE_RESPONSE:
        arn = self._get_param("resourceArn")
        tags = self._get_param("tags", [])
        self.stepfunction_backend.tag_resource(arn, tags)
        return 200, {}, json.dumps({})

    def untag_resource(self) -> TYPE_RESPONSE:
        arn = self._get_param("resourceArn")
        tag_keys = self._get_param("tagKeys", [])
        self.stepfunction_backend.untag_resource(arn, tag_keys)
        return 200, {}, json.dumps({})

    def start_execution(self) -> TYPE_RESPONSE:
        arn = self._get_param("stateMachineArn")
        name = self._get_param("name")
        execution_input = self._get_param("input", if_none="{}")
        execution = self.stepfunction_backend.start_execution(
            arn, name, execution_input
        )
        response = {
            "executionArn": execution.execution_arn,
            "startDate": iso_8601_datetime_with_milliseconds(execution.start_date),
        }
        return 200, {}, json.dumps(response)

    def list_executions(self) -> TYPE_RESPONSE:
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")
        arn = self._get_param("stateMachineArn")
        status_filter = self._get_param("statusFilter")
        state_machine = self.stepfunction_backend.describe_state_machine(arn)
        results, next_token = self.stepfunction_backend.list_executions(
            arn,
            status_filter=status_filter,
            max_results=max_results,
            next_token=next_token,
        )
        executions = []
        for execution in results:
            result = {
                "executionArn": execution.execution_arn,
                "name": execution.name,
                "startDate": iso_8601_datetime_with_milliseconds(execution.start_date),
                "stateMachineArn": state_machine.arn,
                "status": execution.status,
            }
            if execution.status in [
                ExecutionStatus.SUCCEEDED,
                ExecutionStatus.FAILED,
                ExecutionStatus.ABORTED,
            ]:
                result["stopDate"] = iso_8601_datetime_with_milliseconds(
                    execution.stop_date
                )
            executions.append(result)
        response = {"executions": executions, "nextToken": next_token}
        return 200, {}, json.dumps(response)

    def describe_execution(self) -> TYPE_RESPONSE:
        arn = self._get_param("executionArn")
        execution = self.stepfunction_backend.describe_execution(arn)
        response = {
            "executionArn": arn,
            "input": json.dumps(execution.execution_input),
            "name": execution.name,
            "startDate": iso_8601_datetime_with_milliseconds(execution.start_date),
            "stateMachineArn": execution.state_machine_arn,
            "status": execution.status,
        }
        if execution.status in [
            ExecutionStatus.SUCCEEDED,
            ExecutionStatus.ABORTED,
            ExecutionStatus.FAILED,
        ]:
            response["stopDate"] = iso_8601_datetime_with_milliseconds(
                execution.stop_date
            )
        if execution.status in [
            ExecutionStatus.SUCCEEDED,
            ExecutionStatus.SUCCEEDED.value,
        ]:
            if isinstance(execution.output, str):
                response["output"] = execution.output
            elif execution.output is not None:
                response["output"] = json.dumps(execution.output)
            response["outputDetails"] = execution.output_details
        if execution.error is not None:
            response["error"] = execution.error
        if execution.cause is not None:
            response["cause"] = execution.cause
        return 200, {}, json.dumps(response)

    def describe_state_machine_for_execution(self) -> TYPE_RESPONSE:
        arn = self._get_param("executionArn")
        sm = self.stepfunction_backend.describe_state_machine_for_execution(arn)
        return self._describe_state_machine(sm.arn)

    def stop_execution(self) -> TYPE_RESPONSE:
        arn = self._get_param("executionArn")
        execution = self.stepfunction_backend.stop_execution(arn)
        response = {
            "stopDate": iso_8601_datetime_with_milliseconds(execution.stop_date)
        }
        return 200, {}, json.dumps(response)

    def get_execution_history(self) -> TYPE_RESPONSE:
        execution_arn = self._get_param("executionArn")
        execution_history = self.stepfunction_backend.get_execution_history(
            execution_arn
        )
        return 200, {}, json.dumps(execution_history)

    def send_task_failure(self) -> TYPE_RESPONSE:
        task_token = self._get_param("taskToken")
        error = self._get_param("error")
        self.stepfunction_backend.send_task_failure(task_token, error=error)
        return 200, {}, "{}"

    def send_task_heartbeat(self) -> TYPE_RESPONSE:
        task_token = self._get_param("taskToken")
        self.stepfunction_backend.send_task_heartbeat(task_token)
        return 200, {}, "{}"

    def send_task_success(self) -> TYPE_RESPONSE:
        task_token = self._get_param("taskToken")
        output = self._get_param("output")
        self.stepfunction_backend.send_task_success(task_token, output)
        return 200, {}, "{}"

    def list_map_runs(self) -> TYPE_RESPONSE:
        execution_arn = self._get_param("executionArn")
        runs = self.stepfunction_backend.list_map_runs(execution_arn)
        return 200, {}, json.dumps(runs)

    def describe_map_run(self) -> TYPE_RESPONSE:
        map_run_arn = self._get_param("mapRunArn")
        run = self.stepfunction_backend.describe_map_run(map_run_arn)
        return 200, {}, json.dumps(run)

    def update_map_run(self) -> TYPE_RESPONSE:
        map_run_arn = self._get_param("mapRunArn")
        max_concurrency = self._get_param("maxConcurrency")
        tolerated_failure_count = self._get_param("toleratedFailureCount")
        tolerated_failure_percentage = self._get_param("toleratedFailurePercentage")
        self.stepfunction_backend.update_map_run(
            map_run_arn,
            max_concurrency,
            tolerated_failure_count=tolerated_failure_count,
            tolerated_failure_percentage=tolerated_failure_percentage,
        )
        return 200, {}, "{}"

    def create_activity(self) -> TYPE_RESPONSE:
        name = self._get_param("name")
        tags = self._get_param("tags")
        encryption_configuration = self._get_param("encryptionConfiguration")
        activity = self.stepfunction_backend.create_activity(
            name=name, tags=tags, encryption_configuration=encryption_configuration
        )
        response = {
            "creationDate": activity.creation_date,
            "activityArn": activity.arn,
        }
        return 200, {}, json.dumps(response)

    def describe_activity(self) -> TYPE_RESPONSE:
        activity_arn = self._get_param("activityArn")
        activity = self.stepfunction_backend.describe_activity(activity_arn)

        response = {
            "activityArn": activity.arn,
            "name": activity.name,
            "creationDate": activity.creation_date,
            "encryptionConfiguration": activity.encryption_configuration,
        }
        return 200, {}, json.dumps(response)

    def delete_activity(self) -> TYPE_RESPONSE:
        activity_arn = self._get_param("activityArn")
        self.stepfunction_backend.delete_activity(activity_arn)
        return 200, {}, "{}"

    def list_activities(self) -> TYPE_RESPONSE:
        max_results = self._get_param("maxResults")
        next_token = self._get_param("nextToken")
        results, next_token = self.stepfunction_backend.list_activities(
            max_results=max_results, next_token=next_token
        )

        activities = [
            {
                "activityArn": activity.arn,
                "name": activity.name,
                "creationDate": activity.creation_date,
            }
            for activity in results
        ]
        response = {"activities": activities, "nextToken": next_token}
        return 200, {}, json.dumps(response)
