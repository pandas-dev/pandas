import json
import re
from datetime import datetime, timezone
from re import Pattern
from typing import Any, Optional

from moto import settings
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import CloudFormationModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from .exceptions import (
    ActivityAlreadyExists,
    ActivityDoesNotExist,
    ExecutionAlreadyExists,
    ExecutionDoesNotExist,
    InvalidArn,
    InvalidEncryptionConfiguration,
    InvalidExecutionInput,
    InvalidName,
    NameTooLongException,
    ResourceNotFound,
    StateMachineDoesNotExist,
)
from .parser.api import EncryptionType
from .utils import PAGINATION_MODEL, api_to_cfn_tags, cfn_to_api_tags


class StateMachineInstance:
    def __init__(
        self,
        arn: str,
        name: str,
        definition: str,
        roleArn: str,
        encryptionConfiguration: Optional[dict[str, Any]] = None,
        loggingConfiguration: Optional[dict[str, Any]] = None,
        tracingConfiguration: Optional[dict[str, Any]] = None,
    ):
        self.creation_date = iso_8601_datetime_with_milliseconds()
        self.update_date = self.creation_date
        self.arn = arn
        self.name = name
        self.definition = definition
        self.roleArn = roleArn
        self.executions: list[Execution] = []
        self.type = "STANDARD"
        self.encryptionConfiguration = encryptionConfiguration or {
            "type": "AWS_OWNED_KEY"
        }
        self.loggingConfiguration = loggingConfiguration or {"level": "OFF"}
        self.tracingConfiguration = tracingConfiguration or {"enabled": False}
        self.sm_type = "STANDARD"  # or express
        self.description: Optional[str] = None


class StateMachineVersion(StateMachineInstance, CloudFormationModel):
    def __init__(
        self, source: StateMachineInstance, version: int, description: Optional[str]
    ):
        version_arn = f"{source.arn}:{version}"
        StateMachineInstance.__init__(
            self,
            arn=version_arn,
            name=source.name,
            definition=source.definition,
            roleArn=source.roleArn,
            encryptionConfiguration=source.encryptionConfiguration,
            loggingConfiguration=source.loggingConfiguration,
            tracingConfiguration=source.tracingConfiguration,
        )
        self.source_arn = source.arn
        self.version = version
        self.description = description


class StateMachine(StateMachineInstance, CloudFormationModel):
    def __init__(
        self,
        arn: str,
        name: str,
        definition: str,
        roleArn: str,
        backend: "StepFunctionBackend",
        encryptionConfiguration: Optional[dict[str, Any]] = None,
        loggingConfiguration: Optional[dict[str, Any]] = None,
        tracingConfiguration: Optional[dict[str, Any]] = None,
    ):
        StateMachineInstance.__init__(
            self,
            arn=arn,
            name=name,
            definition=definition,
            roleArn=roleArn,
            encryptionConfiguration=encryptionConfiguration,
            loggingConfiguration=loggingConfiguration,
            tracingConfiguration=tracingConfiguration,
        )
        self.latest_version_number = 0
        self.versions: dict[int, StateMachineVersion] = {}
        self.latest_version: Optional[StateMachineVersion] = None
        self.backend = backend

    def publish(self, description: Optional[str]) -> None:
        new_version_number = self.latest_version_number + 1
        new_version = StateMachineVersion(
            source=self, version=new_version_number, description=description
        )
        self.versions[new_version_number] = new_version
        self.latest_version = new_version
        self.latest_version_number = new_version_number

    def start_execution(
        self,
        region_name: str,
        account_id: str,
        execution_name: str,
        execution_input: str,
    ) -> "Execution":
        self._validate_execution_input(execution_input)
        existing_execution = self._handle_name_input_idempotency(
            execution_name, execution_input
        )
        if existing_execution is not None:
            # If we found a match for the name and input, return the existing execution.
            return existing_execution

        execution = Execution(
            region_name=region_name,
            account_id=account_id,
            state_machine_name=self.name,
            execution_name=execution_name,
            state_machine_arn=self.arn,
            execution_input=execution_input,
        )
        self.executions.append(execution)
        return execution

    def stop_execution(self, execution_arn: str) -> "Execution":
        execution = next(
            (x for x in self.executions if x.execution_arn == execution_arn), None
        )
        if not execution:
            raise ExecutionDoesNotExist(
                "Execution Does Not Exist: '" + execution_arn + "'"
            )
        execution.stop(stop_date=datetime.now(), error="", cause="")
        return execution

    def _handle_name_input_idempotency(
        self, name: str, execution_input: str
    ) -> Optional["Execution"]:
        for execution in self.executions:
            if execution.name == name:
                # Executions with the same name and input are considered idempotent
                if (
                    execution_input == execution.execution_input
                    and execution.status == "RUNNING"
                ):
                    return execution

                # If the inputs are different _or_ the execution already finished, raise
                raise ExecutionAlreadyExists(
                    "Execution Already Exists: '" + execution.execution_arn + "'"
                )
        return None

    def _validate_execution_input(self, execution_input: str) -> None:
        try:
            json.loads(execution_input)
        except Exception as ex:
            raise InvalidExecutionInput(
                "Invalid State Machine Execution Input: '" + str(ex) + "'"
            )

    def update(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            if value is not None:
                setattr(self, key, value)
        self.update_date = iso_8601_datetime_with_milliseconds()

    @property
    def physical_resource_id(self) -> str:
        return self.arn

    def get_cfn_properties(self, prop_overrides: dict[str, Any]) -> dict[str, Any]:
        property_names = [
            "DefinitionString",
            "RoleArn",
            "StateMachineName",
        ]
        properties = {}
        for prop in property_names:
            properties[prop] = prop_overrides.get(prop, self.get_cfn_attribute(prop))
        # Special handling for Tags
        overridden_keys = [tag["Key"] for tag in prop_overrides.get("Tags", [])]
        original_tags_to_include = [
            tag
            for tag in self.get_cfn_attribute("Tags")
            if tag["Key"] not in overridden_keys
        ]
        properties["Tags"] = original_tags_to_include + prop_overrides.get("Tags", [])
        return properties

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in [
            "Name",
            "DefinitionString",
            "RoleArn",
            "StateMachineName",
            "Tags",
        ]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Name":
            return self.name
        elif attribute_name == "DefinitionString":
            return self.definition
        elif attribute_name == "RoleArn":
            return self.roleArn
        elif attribute_name == "StateMachineName":
            return self.name
        elif attribute_name == "Tags":
            return api_to_cfn_tags(
                self.backend.get_tags_list_for_state_machine(self.arn)
            )

        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return "StateMachine"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::StepFunctions::StateMachine"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "StateMachine":
        properties = cloudformation_json["Properties"]
        name = properties.get("StateMachineName", resource_name)
        definition = properties.get("DefinitionString", "")
        role_arn = properties.get("RoleArn", "")
        tags = cfn_to_api_tags(properties.get("Tags", []))
        sf_backend = stepfunctions_backends[account_id][region_name]
        return sf_backend.create_state_machine(name, definition, role_arn, tags=tags)

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        sf_backend = stepfunctions_backends[account_id][region_name]
        sf_backend.delete_state_machine(resource_name)

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "StateMachine":
        properties = cloudformation_json.get("Properties", {})
        name = properties.get("StateMachineName", original_resource.name)

        if name != original_resource.name:
            # Replacement
            new_properties = original_resource.get_cfn_properties(properties)
            cloudformation_json["Properties"] = new_properties
            new_resource = cls.create_from_cloudformation_json(
                name, cloudformation_json, account_id, region_name
            )
            cls.delete_from_cloudformation_json(
                original_resource.arn, cloudformation_json, account_id, region_name
            )
            return new_resource

        else:
            # No Interruption
            definition = properties.get("DefinitionString")
            role_arn = properties.get("RoleArn")
            tags = cfn_to_api_tags(properties.get("Tags", []))
            sf_backend = stepfunctions_backends[account_id][region_name]
            state_machine = sf_backend.update_state_machine(
                original_resource.arn, definition=definition, role_arn=role_arn
            )
            sf_backend.tag_resource(state_machine.arn, tags)
            return state_machine


class Execution:
    def __init__(
        self,
        region_name: str,
        account_id: str,
        state_machine_name: str,
        execution_name: str,
        state_machine_arn: str,
        execution_input: str,
    ):
        execution_arn = "arn:{}:states:{}:{}:execution:{}:{}"
        execution_arn = execution_arn.format(
            get_partition(region_name),
            region_name,
            account_id,
            state_machine_name,
            execution_name,
        )
        self.execution_arn = execution_arn
        self.name = execution_name
        self.start_date = datetime.now()
        self.state_machine_arn = state_machine_arn
        self.execution_input = execution_input
        self.status = (
            "RUNNING"
            if settings.get_sf_execution_history_type() == "SUCCESS"
            else "FAILED"
        )
        self.stop_date: Optional[datetime] = None
        self.account_id = account_id
        self.region_name = region_name
        self.output: Optional[str] = None
        self.output_details: Optional[str] = None
        self.cause: Optional[str] = None
        self.error: Optional[str] = None

    def get_execution_history(self, roleArn: str) -> list[dict[str, Any]]:
        sf_execution_history_type = settings.get_sf_execution_history_type()
        tzlocal = datetime.now(timezone.utc).astimezone().tzinfo
        if sf_execution_history_type == "SUCCESS":
            return [
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 0, tzinfo=tzlocal)
                    ),
                    "type": "ExecutionStarted",
                    "id": 1,
                    "previousEventId": 0,
                    "executionStartedEventDetails": {
                        "input": "{}",
                        "inputDetails": {"truncated": False},
                        "roleArn": roleArn,
                    },
                },
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 10, tzinfo=tzlocal)
                    ),
                    "type": "PassStateEntered",
                    "id": 2,
                    "previousEventId": 0,
                    "stateEnteredEventDetails": {
                        "name": "A State",
                        "input": "{}",
                        "inputDetails": {"truncated": False},
                    },
                },
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 10, tzinfo=tzlocal)
                    ),
                    "type": "PassStateExited",
                    "id": 3,
                    "previousEventId": 2,
                    "stateExitedEventDetails": {
                        "name": "A State",
                        "output": "An output",
                        "outputDetails": {"truncated": False},
                    },
                },
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 20, tzinfo=tzlocal)
                    ),
                    "type": "ExecutionSucceeded",
                    "id": 4,
                    "previousEventId": 3,
                    "executionSucceededEventDetails": {
                        "output": "An output",
                        "outputDetails": {"truncated": False},
                    },
                },
            ]
        elif sf_execution_history_type == "FAILURE":
            return [
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 0, tzinfo=tzlocal)
                    ),
                    "type": "ExecutionStarted",
                    "id": 1,
                    "previousEventId": 0,
                    "executionStartedEventDetails": {
                        "input": "{}",
                        "inputDetails": {"truncated": False},
                        "roleArn": roleArn,
                    },
                },
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 10, tzinfo=tzlocal)
                    ),
                    "type": "FailStateEntered",
                    "id": 2,
                    "previousEventId": 0,
                    "stateEnteredEventDetails": {
                        "name": "A State",
                        "input": "{}",
                        "inputDetails": {"truncated": False},
                    },
                },
                {
                    "timestamp": iso_8601_datetime_with_milliseconds(
                        datetime(2020, 1, 1, 0, 0, 10, tzinfo=tzlocal)
                    ),
                    "type": "ExecutionFailed",
                    "id": 3,
                    "previousEventId": 2,
                    "executionFailedEventDetails": {
                        "error": "AnError",
                        "cause": "An error occurred!",
                    },
                },
            ]
        return []

    def stop(self, *args: Any, **kwargs: Any) -> None:
        self.status = "ABORTED"
        self.stop_date = datetime.now()


class Activity:
    def __init__(
        self,
        arn: str,
        name: str,
        encryption_configuration: Optional[dict[str, Any]] = None,
    ):
        self.arn = arn
        self.name = name
        self.encryption_configuration = encryption_configuration

        self.creation_date = iso_8601_datetime_with_milliseconds()
        self.update_date = self.creation_date


class StepFunctionBackend(BaseBackend):
    """
    Configure Moto to explicitly parse and execute the StateMachine:

    .. sourcecode:: python

        @mock_aws(config={"stepfunctions": {"execute_state_machine": True}})

    By default, executing a StateMachine does nothing, and calling `describe_state_machine` will return static data.

    Set the following environment variable if you want to get the static data to have a FAILED status:

    .. sourcecode:: bash

        SF_EXECUTION_HISTORY_TYPE=FAILURE

    """

    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/stepfunctions.html#SFN.Client.create_state_machine
    # A name must not contain:
    # whitespace
    # brackets < > { } [ ]
    # wildcard characters ? *
    # special characters " # % \ ^ | ~ ` $ & , ; : /
    invalid_chars_for_name = [
        " ",
        "{",
        "}",
        "[",
        "]",
        "<",
        ">",
        "?",
        "*",
        '"',
        "#",
        "%",
        "\\",
        "^",
        "|",
        "~",
        "`",
        "$",
        "&",
        ",",
        ";",
        ":",
        "/",
    ]
    # control characters (U+0000-001F , U+007F-009F )
    invalid_unicodes_for_name = [
        "\u0000",
        "\u0001",
        "\u0002",
        "\u0003",
        "\u0004",
        "\u0005",
        "\u0006",
        "\u0007",
        "\u0008",
        "\u0009",
        "\u000a",
        "\u000b",
        "\u000c",
        "\u000d",
        "\u000e",
        "\u000f",
        "\u0010",
        "\u0011",
        "\u0012",
        "\u0013",
        "\u0014",
        "\u0015",
        "\u0016",
        "\u0017",
        "\u0018",
        "\u0019",
        "\u001a",
        "\u001b",
        "\u001c",
        "\u001d",
        "\u001e",
        "\u001f",
        "\u007f",
        "\u0080",
        "\u0081",
        "\u0082",
        "\u0083",
        "\u0084",
        "\u0085",
        "\u0086",
        "\u0087",
        "\u0088",
        "\u0089",
        "\u008a",
        "\u008b",
        "\u008c",
        "\u008d",
        "\u008e",
        "\u008f",
        "\u0090",
        "\u0091",
        "\u0092",
        "\u0093",
        "\u0094",
        "\u0095",
        "\u0096",
        "\u0097",
        "\u0098",
        "\u0099",
        "\u009a",
        "\u009b",
        "\u009c",
        "\u009d",
        "\u009e",
        "\u009f",
    ]
    accepted_role_arn_format = re.compile(
        ARN_PARTITION_REGEX + r":iam::(?P<account_id>[0-9]{12}):role/.+"
    )
    accepted_mchn_arn_format = re.compile(
        ARN_PARTITION_REGEX
        + r":states:[-0-9a-zA-Z]+:(?P<account_id>[0-9]{12}):stateMachine:.+"
    )
    accepted_exec_arn_format = re.compile(
        ARN_PARTITION_REGEX
        + r":states:[-0-9a-zA-Z]+:(?P<account_id>[0-9]{12}):execution:.+"
    )
    accepted_activity_arn_format = re.compile(
        ARN_PARTITION_REGEX
        + r":states:[-0-9a-zA-Z]+:(?P<account_id>[0-9]{12}):activity:.+"
    )

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.tagger = TaggingService(
            tag_name="tags", key_name="key", value_name="value"
        )

        self.state_machines: list[StateMachine] = []
        self.activities: dict[str, Activity] = {}
        self._account_id = None

    def create_state_machine(
        self,
        name: str,
        definition: str,
        roleArn: str,
        tags: Optional[list[dict[str, str]]] = None,
        publish: Optional[bool] = None,
        loggingConfiguration: Optional[dict[str, Any]] = None,
        tracingConfiguration: Optional[dict[str, Any]] = None,
        encryptionConfiguration: Optional[dict[str, Any]] = None,
        version_description: Optional[str] = None,
    ) -> StateMachine:
        self._validate_name(name)
        self._validate_role_arn(roleArn)
        arn = f"arn:{get_partition(self.region_name)}:states:{self.region_name}:{self.account_id}:stateMachine:{name}"
        try:
            return self.describe_state_machine(arn)
        except StateMachineDoesNotExist:
            state_machine = StateMachine(
                arn,
                name,
                definition,
                roleArn,
                self,
                encryptionConfiguration,
                loggingConfiguration,
                tracingConfiguration,
            )
            if publish:
                state_machine.publish(description=version_description)

            if tags:
                self.tag_resource(arn, tags)

            self.state_machines.append(state_machine)
            return state_machine

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_state_machines(self) -> list[StateMachine]:
        return sorted(self.state_machines, key=lambda x: x.creation_date)

    def describe_state_machine(self, arn: str) -> StateMachine:
        self._validate_machine_arn(arn)
        sm = next((x for x in self.state_machines if x.arn == arn), None)
        if not sm:
            if (
                (arn_parts := arn.split(":"))
                and len(arn_parts) > 7
                and arn_parts[-1].isnumeric()
            ):
                # we might have a versioned arn, ending in :stateMachine:name:version_nr
                source_arn = ":".join(arn_parts[:-1])
                source_sm = next(
                    (x for x in self.state_machines if x.arn == source_arn), None
                )
                if source_sm:
                    sm = source_sm.versions.get(int(arn_parts[-1]))  # type: ignore[assignment]
        if not sm:
            raise StateMachineDoesNotExist(f"State Machine Does Not Exist: '{arn}'")
        return sm  # type: ignore[return-value]

    def delete_state_machine(self, arn: str) -> None:
        self._validate_machine_arn(arn)
        sm = next((x for x in self.state_machines if x.arn == arn), None)
        if sm:
            self.state_machines.remove(sm)

    def update_state_machine(
        self,
        arn: str,
        definition: Optional[str] = None,
        role_arn: Optional[str] = None,
        logging_configuration: Optional[dict[str, bool]] = None,
        tracing_configuration: Optional[dict[str, bool]] = None,
        encryption_configuration: Optional[dict[str, Any]] = None,
        publish: Optional[bool] = None,
        version_description: Optional[str] = None,
    ) -> StateMachine:
        sm = self.describe_state_machine(arn)
        updates: dict[str, Any] = {
            "definition": definition,
            "roleArn": role_arn,
        }
        if encryption_configuration:
            updates["encryptionConfiguration"] = encryption_configuration
        if logging_configuration:
            updates["loggingConfiguration"] = logging_configuration
        if tracing_configuration:
            updates["tracingConfiguration"] = tracing_configuration
        sm.update(**updates)
        if publish:
            sm.publish(version_description)
        return sm

    def start_execution(
        self, state_machine_arn: str, name: str, execution_input: str
    ) -> Execution:
        if name:
            self._validate_name(name)
        state_machine = self.describe_state_machine(state_machine_arn)
        return state_machine.start_execution(
            region_name=self.region_name,
            account_id=self.account_id,
            execution_name=name or str(mock_random.uuid4()),
            execution_input=execution_input,
        )

    def stop_execution(self, execution_arn: str) -> Execution:
        self._validate_execution_arn(execution_arn)
        state_machine = self._get_state_machine_for_execution(execution_arn)
        return state_machine.stop_execution(execution_arn)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_executions(
        self, state_machine_arn: str, status_filter: Optional[str] = None
    ) -> list[Execution]:
        executions = self.describe_state_machine(state_machine_arn).executions

        if status_filter:
            executions = list(filter(lambda e: e.status == status_filter, executions))

        return sorted(executions, key=lambda x: x.start_date, reverse=True)

    def describe_execution(self, execution_arn: str) -> Execution:
        self._validate_execution_arn(execution_arn)
        state_machine = self._get_state_machine_for_execution(execution_arn)
        exctn = next(
            (x for x in state_machine.executions if x.execution_arn == execution_arn),
            None,
        )
        if not exctn:
            raise ExecutionDoesNotExist(
                "Execution Does Not Exist: '" + execution_arn + "'"
            )
        return exctn

    def get_execution_history(self, execution_arn: str) -> dict[str, Any]:
        self._validate_execution_arn(execution_arn)
        state_machine = self._get_state_machine_for_execution(execution_arn)
        execution = next(
            (x for x in state_machine.executions if x.execution_arn == execution_arn),
            None,
        )
        if not execution:
            raise ExecutionDoesNotExist(
                "Execution Does Not Exist: '" + execution_arn + "'"
            )
        return {"events": execution.get_execution_history(state_machine.roleArn)}

    def describe_state_machine_for_execution(self, execution_arn: str) -> StateMachine:
        for sm in self.state_machines:
            for exc in sm.executions:
                if exc.execution_arn == execution_arn:
                    return sm
        raise ResourceNotFound(execution_arn)

    def list_tags_for_resource(self, arn: str) -> dict[str, list[dict[str, str]]]:
        return self.tagger.list_tags_for_resource(arn)

    def tag_resource(self, resource_arn: str, tags: list[dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def get_tags_list_for_state_machine(self, arn: str) -> list[dict[str, str]]:
        return self.list_tags_for_resource(arn)[self.tagger.tag_name]

    def send_task_failure(self, task_token: str, error: Optional[str] = None) -> None:
        pass

    def send_task_heartbeat(self, task_token: str) -> None:
        pass

    def send_task_success(self, task_token: str, outcome: str) -> None:
        pass

    def describe_map_run(self, map_run_arn: str) -> dict[str, Any]:
        return {}

    def list_map_runs(self, execution_arn: str) -> Any:
        return []

    def update_map_run(
        self,
        map_run_arn: str,
        max_concurrency: int,
        tolerated_failure_count: str,
        tolerated_failure_percentage: str,
    ) -> None:
        pass

    def _validate_name(self, name: str) -> None:
        if any(invalid_char in name for invalid_char in self.invalid_chars_for_name):
            raise InvalidName("Invalid Name: '" + name + "'")

        if any(name.find(char) >= 0 for char in self.invalid_unicodes_for_name):
            raise InvalidName("Invalid Name: '" + name + "'")

        if len(name) > 80:
            raise NameTooLongException(name)

    def _validate_role_arn(self, role_arn: str) -> None:
        self._validate_arn(
            arn=role_arn,
            regex=self.accepted_role_arn_format,
            invalid_msg="Invalid Role Arn: '" + role_arn + "'",
        )

    def _validate_machine_arn(self, machine_arn: str) -> None:
        self._validate_arn(
            arn=machine_arn,
            regex=self.accepted_mchn_arn_format,
            invalid_msg="Invalid State Machine Arn: '" + machine_arn + "'",
        )

    def _validate_execution_arn(self, execution_arn: str) -> None:
        self._validate_arn(
            arn=execution_arn,
            regex=self.accepted_exec_arn_format,
            invalid_msg="Execution Does Not Exist: '" + execution_arn + "'",
        )

    def _validate_activity_arn(self, activity_arn: str) -> None:
        self._validate_arn(
            arn=activity_arn,
            regex=self.accepted_activity_arn_format,
            invalid_msg="Invalid Activity Arn: '" + activity_arn + "'",
        )

    def _validate_arn(self, arn: str, regex: Pattern[str], invalid_msg: str) -> None:
        match = regex.match(arn)
        if not arn or not match:
            raise InvalidArn(invalid_msg)

    def _validate_encryption_configuration(
        self, encryption_configuration: dict[str, Any]
    ) -> None:
        encryption_type = encryption_configuration.get("type")
        if not encryption_type:
            raise InvalidEncryptionConfiguration(
                "Invalid Encryption Configuration: 'type' is required."
            )

        if (
            encryption_type == EncryptionType.CUSTOMER_MANAGED_KMS_KEY
            and not encryption_configuration.get("kmsKeyId")
        ):
            raise InvalidEncryptionConfiguration(
                "Invalid Encryption Configuration: 'kmsKeyId' is required when 'type' is 'CUSTOMER_MANAGED_KMS_KEY'"
            )

    def _get_state_machine_for_execution(self, execution_arn: str) -> StateMachine:
        state_machine_name = execution_arn.split(":")[6]
        state_machine_arn = next(
            (x.arn for x in self.state_machines if x.name == state_machine_name), None
        )
        if not state_machine_arn:
            # Assume that if the state machine arn is not present, then neither will the
            # execution
            raise ExecutionDoesNotExist(
                "Execution Does Not Exist: '" + execution_arn + "'"
            )
        return self.describe_state_machine(state_machine_arn)

    def create_activity(
        self,
        name: str,
        tags: Optional[list[dict[str, str]]] = None,
        encryption_configuration: Optional[dict[str, Any]] = None,
    ) -> Activity:
        self._validate_name(name)

        if encryption_configuration is not None:
            self._validate_encryption_configuration(encryption_configuration)

        arn = f"arn:{get_partition(self.region_name)}:states:{self.region_name}:{self.account_id}:activity:{name}"

        if arn in self.activities:
            raise ActivityAlreadyExists("Activity already exists.")

        activity = Activity(
            arn=arn,
            name=name,
            encryption_configuration=encryption_configuration,
        )
        self.activities[arn] = activity

        if tags:
            self.tag_resource(arn, tags)

        return activity

    def describe_activity(self, activity_arn: str) -> Activity:
        self._validate_activity_arn(activity_arn)
        if activity_arn in self.activities:
            return self.activities[activity_arn]
        else:
            raise ActivityDoesNotExist(activity_arn)

    def delete_activity(self, activity_arn: str) -> None:
        self._validate_activity_arn(activity_arn)
        if activity_arn in self.activities:
            del self.activities[activity_arn]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_activities(self) -> list[Activity]:
        return sorted(self.activities.values(), key=lambda x: x.creation_date)


stepfunctions_backends = BackendDict(StepFunctionBackend, "stepfunctions")
