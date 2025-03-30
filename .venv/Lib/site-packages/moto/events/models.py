import copy
import json
import re
import sys
import warnings
from collections import OrderedDict
from enum import Enum, unique
from json import JSONDecodeError
from operator import eq, ge, gt, le, lt
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import requests

from moto import settings
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.exceptions import JsonRESTError
from moto.core.utils import (
    iso_8601_datetime_without_milliseconds,
    unix_time,
    unix_time_millis,
    utcfromtimestamp,
)
from moto.events.exceptions import (
    IllegalStatusException,
    InvalidEventPatternException,
    ResourceAlreadyExistsException,
    ResourceNotFoundException,
    ValidationException,
)
from moto.moto_api._internal import mock_random as random
from moto.utilities.arns import parse_arn
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import (
    ARN_PARTITION_REGEX,
    CamelToUnderscoresWalker,
    get_partition,
)

from .utils import (
    _BASE_EVENT_MESSAGE,
    PAGINATION_MODEL,
    EventMessageType,
    EventTemplateParser,
)

if TYPE_CHECKING:
    from moto.secretsmanager.models import SecretsManagerBackend

# Sentinel to signal the absence of a field for `Exists` pattern matching
UNDEFINED = object()


def get_secrets_manager_backend(
    account_id: str, region: str
) -> "SecretsManagerBackend":
    from moto.secretsmanager import secretsmanager_backends

    return secretsmanager_backends[account_id][region]


class Rule(CloudFormationModel):
    def __init__(
        self,
        name: str,
        account_id: str,
        region_name: str,
        description: Optional[str],
        event_pattern: Optional[str],
        schedule_exp: Optional[str],
        role_arn: Optional[str],
        event_bus_name: str,
        state: Optional[str],
        managed_by: Optional[str] = None,
        targets: Optional[List[Dict[str, Any]]] = None,
    ):
        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        self.description = description
        self.event_pattern = EventPattern.load(event_pattern)
        self.scheduled_expression = schedule_exp
        self.role_arn = role_arn
        self.event_bus_name = event_bus_name
        self.state = state or "ENABLED"
        self.managed_by = managed_by  # can only be set by AWS services
        self.created_by = account_id
        self.targets = targets or []

    @property
    def arn(self) -> str:
        event_bus_name = (
            "" if self.event_bus_name == "default" else f"{self.event_bus_name}/"
        )

        return f"arn:{get_partition(self.region_name)}:events:{self.region_name}:{self.account_id}:rule/{event_bus_name}{self.name}"

    @property
    def physical_resource_id(self) -> str:
        return self.name

    # This song and dance for targets is because we need order for Limits and NextTokens, but can't use OrderedDicts
    # with Python 2.6, so tracking it with an array it is.
    def _check_target_exists(self, target_id: str) -> Optional[int]:
        for i in range(0, len(self.targets)):
            if target_id == self.targets[i]["Id"]:
                return i
        return None

    def enable(self) -> None:
        self.state = "ENABLED"

    def disable(self) -> None:
        self.state = "DISABLED"

    def delete(self, account_id: str, region_name: str) -> None:
        event_backend: EventsBackend = events_backends[account_id][region_name]
        self.remove_targets([t["Id"] for t in self.targets])
        event_backend.delete_rule(name=self.name, event_bus_arn=self.event_bus_name)

    def put_targets(self, targets: List[Dict[str, Any]]) -> None:
        # Not testing for valid ARNs.
        for target in targets:
            index = self._check_target_exists(target["Id"])
            if index is not None:
                self.targets[index] = target
            else:
                self.targets.append(target)

    def remove_targets(self, ids: List[str]) -> None:
        for target_id in ids:
            index = self._check_target_exists(target_id)
            if index is not None:
                self.targets.pop(index)

    def send_to_targets(
        self, original_event: EventMessageType, transform_input: bool = True
    ) -> None:
        if not self.event_pattern.matches_event(original_event):
            return

        # supported targets
        # - CloudWatch Log Group
        # - EventBridge Archive
        # - SQS Queue + FIFO Queue
        # - Cross-region/account EventBus
        for target in self.targets:
            arn = parse_arn(target["Arn"])

            if transform_input:
                input_transformer = target.get("InputTransformer", {})
                event = EventTemplateParser.parse(
                    input_template=input_transformer.get("InputTemplate"),
                    input_paths_map=input_transformer.get("InputPathsMap", {}),
                    event=original_event,
                )
            else:
                event = original_event.copy()  # type: ignore[assignment]

            if arn.service == "logs" and arn.resource_type == "log-group":
                self._send_to_cw_log_group(arn.resource_id, event)
            elif arn.service == "events" and not arn.resource_type:
                archive_arn = parse_arn(event["archive-arn"])

                self._send_to_events_archive(archive_arn.resource_id, original_event)
            elif arn.service == "sqs":
                group_id = target.get("SqsParameters", {}).get("MessageGroupId")
                self._send_to_sqs_queue(arn.resource_id, event, group_id)
            elif arn.service == "events" and arn.resource_type == "event-bus":
                cross_account_backend: EventsBackend = events_backends[arn.account][
                    arn.region
                ]
                new_event = {
                    "Source": event["source"],
                    "DetailType": event["detail-type"],
                    "Detail": json.dumps(event["detail"]),
                    "EventBusName": arn.resource_id,
                }
                cross_account_backend.put_events([new_event])
            elif arn.service == "events" and arn.resource_type == "api-destination":
                if settings.events_invoke_http():
                    api_destination = self._find_api_destination(arn.resource_id)
                    request_parameters = target.get("HttpParameters", {})
                    headers = request_parameters.get("HeaderParameters", {})
                    qs_params = request_parameters.get("QueryStringParameters", {})
                    query_string = "&".join(
                        [f"{key}={val}" for key, val in qs_params.items()]
                    )
                    url = api_destination.invocation_endpoint + (
                        f"?{query_string}" if query_string else ""
                    )
                    requests.request(
                        method=api_destination.http_method,
                        url=url,
                        headers=headers,
                    )
            else:
                raise NotImplementedError(f"Expr not defined for {type(self)}")

    def _send_to_cw_log_group(self, name: str, event: Dict[str, Any]) -> None:
        from moto.logs import logs_backends

        event["time"] = iso_8601_datetime_without_milliseconds(
            utcfromtimestamp(event["time"])  # type: ignore[arg-type]
        )

        log_stream_name = str(random.uuid4())
        log_events = [{"timestamp": unix_time_millis(), "message": json.dumps(event)}]

        log_backend = logs_backends[self.account_id][self.region_name]
        log_backend.create_log_stream(name, log_stream_name)
        log_backend.put_log_events(name, log_stream_name, log_events)

    def _send_to_events_archive(
        self, resource_id: str, event: EventMessageType
    ) -> None:
        archive_name, archive_uuid = resource_id.split(":")
        archive = events_backends[self.account_id][self.region_name].archives.get(
            archive_name
        )
        if archive.uuid == archive_uuid:  # type: ignore[union-attr]
            archive.events.append(event)  # type: ignore[union-attr]

    def _find_api_destination(self, resource_id: str) -> "Destination":
        backend: "EventsBackend" = events_backends[self.account_id][self.region_name]
        destination_name = resource_id.split("/")[0]
        return backend.destinations[destination_name]

    def _send_to_sqs_queue(
        self, resource_id: str, event: Dict[str, Any], group_id: Optional[str] = None
    ) -> None:
        from moto.sqs import sqs_backends

        event["time"] = iso_8601_datetime_without_milliseconds(
            utcfromtimestamp(float(event["time"]))  # type: ignore[arg-type]
        )

        if group_id:
            queue_attr = sqs_backends[self.account_id][
                self.region_name
            ].get_queue_attributes(
                queue_name=resource_id, attribute_names=["ContentBasedDeduplication"]
            )
            if queue_attr["ContentBasedDeduplication"] == "false":
                warnings.warn(
                    "To let EventBridge send messages to your SQS FIFO queue, "
                    "you must enable content-based deduplication."
                )
                return

        sqs_backends[self.account_id][self.region_name].send_message(
            queue_name=resource_id,
            message_body=json.dumps(event),
            group_id=group_id,
        )

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn"]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn

        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-events-rule.html
        return "AWS::Events::Rule"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Rule":
        properties = cloudformation_json["Properties"]
        properties.setdefault("EventBusName", "default")

        if "EventPattern" in properties:
            properties["EventPattern"] = json.dumps(properties["EventPattern"])

        event_name = resource_name

        event_pattern = properties.get("EventPattern")
        scheduled_expression = properties.get("ScheduleExpression")
        state = properties.get("State")
        desc = properties.get("Description")
        role_arn = properties.get("RoleArn")
        event_bus_arn = properties.get("EventBusName")
        tags = properties.get("Tags")

        backend: "EventsBackend" = events_backends[account_id][region_name]
        rule = backend.put_rule(
            event_name,
            scheduled_expression=scheduled_expression,
            event_pattern=event_pattern,
            state=state,
            description=desc,
            role_arn=role_arn,
            event_bus_arn=event_bus_arn,
            tags=tags,
        )

        targets = properties.get("Targets", [])
        if targets:
            backend.put_targets(
                name=rule.name,
                event_bus_arn=event_bus_arn,
                targets=targets,
            )

        return rule

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Rule":
        original_resource.delete(account_id, region_name)
        return cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        event_backend: EventsBackend = events_backends[account_id][region_name]
        properties = cloudformation_json["Properties"]
        event_bus_arn = properties.get("EventBusName")

        event_backend.delete_rule(resource_name, event_bus_arn)

    def describe(self) -> Dict[str, Any]:
        attributes = {
            "Arn": self.arn,
            "CreatedBy": self.created_by,
            "Description": self.description,
            "EventBusName": self.event_bus_name,
            "EventPattern": self.event_pattern.dump(),
            "ManagedBy": self.managed_by,
            "Name": self.name,
            "RoleArn": self.role_arn,
            "ScheduleExpression": self.scheduled_expression,
            "State": self.state,
        }
        attributes = {
            attr: value for attr, value in attributes.items() if value is not None
        }
        return attributes


class EventBus(CloudFormationModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.account_id = account_id
        self.region = region_name
        self.name = name
        self.arn = f"arn:{get_partition(self.region)}:events:{self.region}:{account_id}:event-bus/{name}"
        self.tags = tags or []

        self._statements: Dict[str, EventBusPolicyStatement] = {}
        self.rules: Dict[str, Rule] = OrderedDict()

    @property
    def policy(self) -> Optional[str]:
        if self._statements:
            policy = {
                "Version": "2012-10-17",
                "Statement": [stmt.describe() for stmt in self._statements.values()],
            }
            return json.dumps(policy)
        return None

    def has_permissions(self) -> bool:
        return len(self._statements) > 0

    def delete(self, account_id: str, region_name: str) -> None:
        event_backend = events_backends[account_id][region_name]
        for rule in self.rules.values():
            rule.delete(account_id, region_name)
        event_backend.delete_event_bus(name=self.name)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn", "Name", "Policy"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "Arn":
            return self.arn
        elif attribute_name == "Name":
            return self.name
        elif attribute_name == "Policy":
            return self.policy

        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-events-eventbus.html
        return "AWS::Events::EventBus"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "EventBus":
        properties = cloudformation_json["Properties"]
        event_backend = events_backends[account_id][region_name]
        event_name = resource_name
        event_source_name = properties.get("EventSourceName")
        return event_backend.create_event_bus(
            name=event_name, event_source_name=event_source_name
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "EventBus":
        original_resource.delete(account_id, region_name)
        return cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        event_backend = events_backends[account_id][region_name]
        event_bus_name = resource_name
        event_backend.delete_event_bus(event_bus_name)

    def _remove_principals_statements(self, *principals: Any) -> None:
        statements_to_delete = set()

        for principal in principals:
            for sid, statement in self._statements.items():
                if statement.principal == principal:
                    statements_to_delete.add(sid)

        # This is done separately to avoid:
        # RuntimeError: dictionary changed size during iteration
        for sid in statements_to_delete:
            del self._statements[sid]

    def add_permission(
        self,
        statement_id: str,
        action: str,
        principal: Dict[str, str],
        condition: Optional[Dict[str, Any]],
    ) -> None:
        self._remove_principals_statements(principal)
        statement = EventBusPolicyStatement(
            sid=statement_id,
            action=action,
            principal=principal,
            condition=condition,
            resource=self.arn,
        )
        self._statements[statement_id] = statement

    def add_policy(self, policy: Dict[str, Any]) -> None:
        policy_statements = policy["Statement"]

        principals = [stmt["Principal"] for stmt in policy_statements]
        self._remove_principals_statements(*principals)

        for new_statement in policy_statements:
            sid = new_statement["Sid"]
            self._statements[sid] = EventBusPolicyStatement.from_dict(new_statement)

    def remove_statement(self, sid: str) -> Optional["EventBusPolicyStatement"]:
        return self._statements.pop(sid, None)

    def remove_statements(self) -> None:
        self._statements.clear()


class EventBusPolicyStatement:
    def __init__(
        self,
        sid: str,
        principal: Dict[str, str],
        action: str,
        resource: str,
        effect: str = "Allow",
        condition: Optional[Dict[str, Any]] = None,
    ):
        self.sid = sid
        self.principal = principal
        self.action = action
        self.resource = resource
        self.effect = effect
        self.condition = condition

    def describe(self) -> Dict[str, Any]:
        statement: Dict[str, Any] = dict(
            Sid=self.sid,
            Effect=self.effect,
            Principal=self.principal,
            Action=self.action,
            Resource=self.resource,
        )

        if self.condition:
            statement["Condition"] = self.condition
        return statement

    @classmethod
    def from_dict(cls, statement_dict: Dict[str, Any]) -> "EventBusPolicyStatement":  # type: ignore[misc]
        params = dict(
            sid=statement_dict["Sid"],
            effect=statement_dict["Effect"],
            principal=statement_dict["Principal"],
            action=statement_dict["Action"],
            resource=statement_dict["Resource"],
        )
        condition = statement_dict.get("Condition")
        if condition:
            params["condition"] = condition

        return cls(**params)


class Archive(CloudFormationModel):
    # https://docs.aws.amazon.com/eventbridge/latest/APIReference/API_ListArchives.html#API_ListArchives_RequestParameters
    VALID_STATES = [
        "ENABLED",
        "DISABLED",
        "CREATING",
        "UPDATING",
        "CREATE_FAILED",
        "UPDATE_FAILED",
    ]

    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        source_arn: str,
        description: str,
        event_pattern: str,
        retention: str,
    ):
        self.region = region_name
        self.name = name
        self.source_arn = source_arn
        self.description = description
        self.event_pattern = EventPattern.load(event_pattern)
        self.retention = retention if retention else 0

        self.arn = f"arn:{get_partition(region_name)}:events:{region_name}:{account_id}:archive/{name}"
        self.creation_time = unix_time()
        self.state = "ENABLED"
        self.uuid = str(random.uuid4())

        self.events: List[EventMessageType] = []
        self.event_bus_name = source_arn.split("/")[-1]
        self.rule: Optional[Rule] = None

    def describe_short(self) -> Dict[str, Any]:
        return {
            "ArchiveName": self.name,
            "EventSourceArn": self.source_arn,
            "State": self.state,
            "RetentionDays": self.retention,
            "SizeBytes": sys.getsizeof(self.events) if len(self.events) > 0 else 0,
            "EventCount": len(self.events),
            "CreationTime": self.creation_time,
        }

    def describe(self) -> Dict[str, Any]:
        result = {
            "ArchiveArn": self.arn,
            "Description": self.description,
            "EventPattern": self.event_pattern.dump(),
        }
        result.update(self.describe_short())

        return result

    def update(
        self,
        description: Optional[str],
        event_pattern: Optional[str],
        retention: Optional[str],
    ) -> None:
        if description:
            self.description = description
        if event_pattern:
            self.event_pattern = EventPattern.load(event_pattern)
        if retention:
            self.retention = retention

    def delete(self, account_id: str, region_name: str) -> None:
        event_backend: EventsBackend = events_backends[account_id][region_name]

        if self.rule:
            self.rule.delete(account_id, region_name)

        event_backend.archives.pop(self.name)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["Arn", "ArchiveName"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "ArchiveName":
            return self.name
        elif attribute_name == "Arn":
            return self.arn

        raise UnformattedGetAttTemplateException()

    @staticmethod
    def cloudformation_name_type() -> str:
        return "ArchiveName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-events-archive.html
        return "AWS::Events::Archive"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Archive":
        properties = cloudformation_json["Properties"]
        event_backend = events_backends[account_id][region_name]

        source_arn = properties.get("SourceArn")
        description = properties.get("Description")
        event_pattern = properties.get("EventPattern")
        retention = properties.get("RetentionDays")

        return event_backend.create_archive(
            resource_name, source_arn, description, event_pattern, retention
        )

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "Archive":
        if new_resource_name == original_resource.name:
            properties = cloudformation_json["Properties"]

            original_resource.update(
                properties.get("Description"),
                properties.get("EventPattern"),
                properties.get("Retention"),
            )

            return original_resource
        else:
            original_resource.delete(account_id, region_name)
            return cls.create_from_cloudformation_json(
                new_resource_name, cloudformation_json, account_id, region_name
            )


@unique
class ReplayState(Enum):
    # https://docs.aws.amazon.com/eventbridge/latest/APIReference/API_ListReplays.html#API_ListReplays_RequestParameters
    STARTING = "STARTING"
    RUNNING = "RUNNING"
    CANCELLING = "CANCELLING"
    COMPLETED = "COMPLETED"
    CANCELLED = "CANCELLED"
    FAILED = "FAILED"


class Replay(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        description: str,
        source_arn: str,
        start_time: str,
        end_time: str,
        destination: Dict[str, Any],
    ):
        self.account_id = account_id
        self.region = region_name
        self.name = name
        self.description = description
        self.source_arn = source_arn
        self.event_start_time = start_time
        self.event_end_time = end_time
        self.destination = destination

        self.arn = f"arn:{get_partition(region_name)}:events:{region_name}:{account_id}:replay/{name}"
        self.state = ReplayState.STARTING
        self.start_time = unix_time()
        self.end_time: Optional[float] = None

    def describe_short(self) -> Dict[str, Any]:
        return {
            "ReplayName": self.name,
            "EventSourceArn": self.source_arn,
            "State": self.state.value,
            "EventStartTime": self.event_start_time,
            "EventEndTime": self.event_end_time,
            "ReplayStartTime": self.start_time,
            "ReplayEndTime": self.end_time,
        }

    def describe(self) -> Dict[str, Any]:
        result = {
            "ReplayArn": self.arn,
            "Description": self.description,
            "Destination": self.destination,
        }

        result.update(self.describe_short())

        return result

    def replay_events(self, archive: Archive) -> None:
        event_bus_name = self.destination["Arn"].split("/")[-1]

        for event in archive.events:
            event_backend = events_backends[self.account_id][self.region]
            event_bus = event_backend.describe_event_bus(event_bus_name)
            for rule in event_bus.rules.values():
                rule.send_to_targets(
                    dict(
                        event,
                        **{"id": str(random.uuid4()), "replay-name": self.name},  # type: ignore
                    ),
                    transform_input=False,
                )

        self.state = ReplayState.COMPLETED
        self.end_time = unix_time()


class Connection(BaseModel):
    def __init__(
        self,
        name: str,
        account_id: str,
        region_name: str,
        description: str,
        authorization_type: str,
        auth_parameters: Dict[str, Any],
    ):
        self.uuid = random.uuid4()
        self.name = name
        self.region = region_name
        self.description = description
        self.authorization_type = authorization_type
        self.auth_parameters = auth_parameters
        self.creation_time = unix_time()
        self.state = "AUTHORIZED"

        connection_id = f"{self.name}/{self.uuid}"
        secretsmanager_backend = get_secrets_manager_backend(account_id, region_name)
        secret_value = {}
        for key, value in self.auth_parameters.items():
            if key == "InvocationHttpParameters":
                secret_value.update({"InvocationHttpParameters": value})
            else:
                secret_value.update(value)
        secret_value = CamelToUnderscoresWalker.parse(secret_value)

        secret = secretsmanager_backend.create_secret(
            name=f"events!connection/{connection_id}/auth",
            secret_string=json.dumps(secret_value),
            replica_regions=[],
            force_overwrite=False,
            secret_binary=None,
            description=f"Auth parameters for Eventbridge connection {connection_id}",
            tags=[{"Key": "aws:secretsmanager:owningService", "Value": "events"}],
            kms_key_id=None,
            client_request_token=None,
        )
        self.secret_arn = json.loads(secret).get("ARN")
        self.arn = f"arn:{get_partition(region_name)}:events:{region_name}:{account_id}:connection/{connection_id}"

    def describe_short(self) -> Dict[str, Any]:
        """
        Create the short description for the Connection object.

        Taken our from the Response Syntax of this API doc:
            - https://docs.aws.amazon.com/eventbridge/latest/APIReference/API_DeleteConnection.html

        Something to consider:
            - The original response also has
                - LastAuthorizedTime (number)
                - LastModifiedTime (number)
            - At the time of implementing this, there was no place where to set/get
            those attributes. That is why they are not in the response.

        Returns:
            dict
        """
        return {
            "ConnectionArn": self.arn,
            "ConnectionState": self.state,
            "CreationTime": self.creation_time,
        }

    def describe(self) -> Dict[str, Any]:
        """
        Create a complete description for the Connection object.

        Taken our from the Response Syntax of this API doc:
            - https://docs.aws.amazon.com/eventbridge/latest/APIReference/API_DescribeConnection.html

        Something to consider:
            - The original response also has:
                - LastAuthorizedTime (number)
                - LastModifiedTime (number)
                - StateReason (string)
            - At the time of implementing this, there was no place where to set/get
            those attributes. That is why they are not in the response.

        Returns:
            dict
        """
        return {
            "AuthorizationType": self.authorization_type,
            "AuthParameters": self.auth_parameters,
            "ConnectionArn": self.arn,
            "ConnectionState": self.state,
            "CreationTime": self.creation_time,
            "Description": self.description,
            "Name": self.name,
            "SecretArn": self.secret_arn,
        }


class Destination(BaseModel):
    def __init__(
        self,
        name: str,
        account_id: str,
        region_name: str,
        description: str,
        connection_arn: str,
        invocation_endpoint: str,
        invocation_rate_limit_per_second: str,
        http_method: str,
    ):
        self.uuid = random.uuid4()
        self.name = name
        self.region = region_name
        self.description = description
        self.connection_arn = connection_arn
        self.invocation_endpoint = invocation_endpoint
        self.invocation_rate_limit_per_second = invocation_rate_limit_per_second
        self.creation_time = unix_time()
        self.http_method = http_method
        self.state = "ACTIVE"
        self.arn = f"arn:{get_partition(region_name)}:events:{region_name}:{account_id}:api-destination/{name}/{self.uuid}"

    def describe(self) -> Dict[str, Any]:
        return {
            "ApiDestinationArn": self.arn,
            "ApiDestinationState": self.state,
            "ConnectionArn": self.connection_arn,
            "CreationTime": self.creation_time,
            "Description": self.description,
            "HttpMethod": self.http_method,
            "InvocationEndpoint": self.invocation_endpoint,
            "InvocationRateLimitPerSecond": self.invocation_rate_limit_per_second,
            "LastModifiedTime": self.creation_time,
            "Name": self.name,
        }

    def describe_short(self) -> Dict[str, Any]:
        return {
            "ApiDestinationArn": self.arn,
            "ApiDestinationState": self.state,
            "CreationTime": self.creation_time,
            "LastModifiedTime": self.creation_time,
        }


class EventPattern:
    def __init__(self, raw_pattern: Optional[str], pattern: Dict[str, Any]):
        self._raw_pattern = raw_pattern
        self._pattern = pattern

    def get_pattern(self) -> Dict[str, Any]:
        return self._pattern

    def matches_event(self, event: EventMessageType) -> bool:
        if not self._pattern:
            return True
        event = json.loads(json.dumps(event))
        return self._does_event_match(event, self._pattern)

    def _does_event_match(
        self, event: EventMessageType, pattern: Dict[str, Any]
    ) -> bool:
        items_and_filters = [(event.get(k, UNDEFINED), v) for k, v in pattern.items()]
        nested_filter_matches = [
            self._does_event_match(item, nested_filter)  # type: ignore
            for item, nested_filter in items_and_filters
            if isinstance(nested_filter, dict)
        ]
        filter_list_matches = [
            self._does_item_match_filters(item, filter_list)
            for item, filter_list in items_and_filters
            if isinstance(filter_list, list)
        ]
        return all(nested_filter_matches + filter_list_matches)

    def _does_item_match_filters(self, item: Any, filters: Any) -> bool:
        allowed_values = [value for value in filters if isinstance(value, str)]
        allowed_values_match = item in allowed_values if allowed_values else True
        full_match = isinstance(item, list) and item == allowed_values
        named_filter_matches = [
            self._does_item_match_named_filter(item, pattern)
            for pattern in filters
            if isinstance(pattern, dict)
        ]
        return (full_match or allowed_values_match) and all(named_filter_matches)

    @staticmethod
    def _does_item_match_named_filter(item: Any, pattern: Any) -> bool:  # type: ignore[misc]
        filter_name, filter_value = list(pattern.items())[0]
        if filter_name == "exists":
            is_leaf_node = not isinstance(item, dict)
            leaf_exists = is_leaf_node and item is not UNDEFINED
            should_exist = filter_value
            return leaf_exists if should_exist else not leaf_exists
        if filter_name == "prefix":
            prefix = filter_value
            return item.startswith(prefix)
        if filter_name == "numeric":
            as_function = {"<": lt, "<=": le, "=": eq, ">=": ge, ">": gt}
            operators_and_values = zip(filter_value[::2], filter_value[1::2])
            numeric_matches = [
                as_function[operator](item, value)
                for operator, value in operators_and_values
            ]
            return all(numeric_matches)
        else:
            warnings.warn(
                f"'{filter_name}' filter logic unimplemented. defaulting to True"
            )
            return True

    @classmethod
    def load(cls, raw_pattern: Optional[str]) -> "EventPattern":
        parser = EventPatternParser(raw_pattern)
        pattern = parser.parse()
        return cls(raw_pattern, pattern)

    def dump(self) -> Optional[str]:
        return self._raw_pattern


class EventPatternParser:
    def __init__(self, pattern: Optional[str]):
        self.pattern = pattern

    def _validate_event_pattern(self, pattern: Dict[str, Any]) -> None:
        # values in the event pattern have to be either a dict or an array
        for attr, value in pattern.items():
            if isinstance(value, dict):
                self._validate_event_pattern(value)
            elif isinstance(value, list):
                if len(value) == 0:
                    raise InvalidEventPatternException(
                        reason="Empty arrays are not allowed"
                    )
            else:
                raise InvalidEventPatternException(
                    reason=f"'{attr}' must be an object or an array"
                )

    def parse(self) -> Dict[str, Any]:
        try:
            parsed_pattern = json.loads(self.pattern) if self.pattern else dict()
            self._validate_event_pattern(parsed_pattern)
            return parsed_pattern
        except JSONDecodeError:
            raise InvalidEventPatternException(reason="Invalid JSON")


class PartnerEventSource(BaseModel):
    def __init__(self, region: str, name: str):
        self.name = name
        self.arn = f"arn:{get_partition(region)}:events:{region}::event-source/aws.partner/{name}"
        self.created_on = unix_time()
        self.accounts: List[str] = []
        self.state = "ACTIVE"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "CreatedBy": self.name.split("/")[0],
            "CreatedOn": self.created_on,
            "Name": self.name,
            "State": self.state,
        }


class EventsBackend(BaseBackend):
    """
    Some Moto services are configured to generate events and send them to EventBridge. See the AWS documentation here:
    https://docs.aws.amazon.com/eventbridge/latest/userguide/eb-service-event.html

    Events that currently supported

     - S3:CreateBucket

    Targets that are currently supported

     - AWSLambda functions

    Please let us know if you want support for an event/target that is not yet listed here.
    """

    ACCOUNT_ID = re.compile(r"^(\d{1,12}|\*)$")
    STATEMENT_ID = re.compile(r"^[a-zA-Z0-9-_]{1,64}$")
    _CRON_REGEX = re.compile(r"^cron\(.*\)")
    _RATE_REGEX = re.compile(r"^rate\(\d*\s(minute|minutes|hour|hours|day|days)\)")

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.next_tokens: Dict[str, int] = {}
        self.event_buses: Dict[str, EventBus] = {}
        self.event_sources: Dict[str, PartnerEventSource] = {}
        self.archives: Dict[str, Archive] = {}
        self.replays: Dict[str, Replay] = {}
        self.tagger = TaggingService()

        self._add_default_event_bus()
        self.connections: Dict[str, Connection] = {}
        self.destinations: Dict[str, Destination] = {}
        self.partner_event_sources: Dict[str, PartnerEventSource] = {}
        self.approved_parent_event_bus_names: List[str] = []

    def _add_default_event_bus(self) -> None:
        self.event_buses["default"] = EventBus(
            self.account_id, self.region_name, "default"
        )

    def _get_event_bus(self, name: str) -> EventBus:
        event_bus_name = name.split(f"{self.account_id}:event-bus/")[-1]

        event_bus = self.event_buses.get(event_bus_name)
        if not event_bus:
            raise ResourceNotFoundException(
                f"Event bus {event_bus_name} does not exist."
            )

        return event_bus

    def _get_replay(self, name: str) -> Replay:
        replay = self.replays.get(name)
        if not replay:
            raise ResourceNotFoundException(f"Replay {name} does not exist.")

        return replay

    def put_rule(
        self,
        name: str,
        description: Optional[str] = None,
        event_bus_arn: Optional[str] = None,
        event_pattern: Optional[str] = None,
        role_arn: Optional[str] = None,
        scheduled_expression: Optional[str] = None,
        state: Optional[str] = None,
        managed_by: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> Rule:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)

        if not event_pattern and not scheduled_expression:
            raise JsonRESTError(
                "ValidationException",
                "Parameter(s) EventPattern or ScheduleExpression must be specified.",
            )

        if scheduled_expression:
            if event_bus_name != "default":
                raise ValidationException(
                    "ScheduleExpression is supported only on the default event bus."
                )

            if not (
                self._CRON_REGEX.match(scheduled_expression)
                or self._RATE_REGEX.match(scheduled_expression)
            ):
                raise ValidationException("Parameter ScheduleExpression is not valid.")

        event_bus = self._get_event_bus(event_bus_name)
        existing_rule = event_bus.rules.get(name)
        targets = existing_rule.targets if existing_rule else list()
        rule = Rule(
            name,
            self.account_id,
            self.region_name,
            description,
            event_pattern,
            scheduled_expression,
            role_arn,
            event_bus_name,
            state,
            managed_by,
            targets=targets,
        )
        event_bus.rules[name] = rule

        if tags:
            self.tagger.tag_resource(rule.arn, tags)

        return rule

    def _normalize_event_bus_arn(self, event_bus_arn: Optional[str]) -> str:
        if event_bus_arn is None:
            return "default"
        return event_bus_arn.split("/")[-1]

    def delete_rule(self, name: str, event_bus_arn: Optional[str]) -> None:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        try:
            event_bus = self._get_event_bus(event_bus_name)
        except ResourceNotFoundException:
            # If the EventBus is deleted, the Rule is also gone
            return
        rule = event_bus.rules.get(name)
        if not rule:
            return
        if len(rule.targets) > 0:
            raise ValidationException("Rule can't be deleted since it has targets.")

        arn = rule.arn
        if self.tagger.has_tags(arn):
            self.tagger.delete_all_tags_for_resource(arn)
        event_bus.rules.pop(name)

    def describe_rule(self, name: str, event_bus_arn: Optional[str]) -> Rule:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        rule = event_bus.rules.get(name)
        if not rule:
            raise ResourceNotFoundException(f"Rule {name} does not exist.")
        return rule

    def disable_rule(self, name: str, event_bus_arn: Optional[str]) -> bool:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        if name in event_bus.rules:
            event_bus.rules[name].disable()
            return True

        return False

    def enable_rule(self, name: str, event_bus_arn: Optional[str]) -> bool:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        if name in event_bus.rules:
            event_bus.rules[name].enable()
            return True

        return False

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_rule_names_by_target(
        self, target_arn: str, event_bus_arn: Optional[str]
    ) -> List[Rule]:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        matching_rules = []

        for _, rule in event_bus.rules.items():
            for target in rule.targets:
                if target["Arn"] == target_arn:
                    matching_rules.append(rule)

        return matching_rules

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_rules(
        self, prefix: Optional[str] = None, event_bus_arn: Optional[str] = None
    ) -> List[Rule]:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        match_string = ".*"
        if prefix is not None:
            match_string = "^" + prefix + match_string

        match_regex = re.compile(match_string)

        matching_rules = []

        for name, rule in event_bus.rules.items():
            if match_regex.match(name):
                matching_rules.append(rule)

        return matching_rules

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_targets_by_rule(  # type: ignore[misc]
        self, rule_id: str, event_bus_arn: Optional[str]
    ) -> List[Dict[str, Any]]:
        # We'll let a KeyError exception be thrown for response to handle if
        # rule doesn't exist.
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        rule = event_bus.rules[rule_id]
        return rule.targets

    def put_targets(
        self, name: str, event_bus_arn: Optional[str], targets: List[Dict[str, Any]]
    ) -> None:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        # super simple ARN check
        invalid_arn = next(
            (
                target["Arn"]
                for target in targets
                if not re.match(r"arn:[\d\w:\-/]*", target["Arn"])
            ),
            None,
        )
        if invalid_arn:
            raise ValidationException(
                f"Parameter {invalid_arn} is not valid. Reason: Provided Arn is not in correct format."
            )

        for target in targets:
            arn = target["Arn"]

            if (
                ":sqs:" in arn
                and arn.endswith(".fifo")
                and not target.get("SqsParameters")
            ):
                raise ValidationException(
                    f"Parameter(s) SqsParameters must be specified for target: {target['Id']}."
                )

        rule = event_bus.rules.get(name)

        if not rule:
            raise ResourceNotFoundException(
                f"Rule {name} does not exist on EventBus {event_bus_name}."
            )

        rule.put_targets(targets)

    def put_events(self, events: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        The following targets are supported at the moment:

         - CloudWatch Log Group
         - EventBridge Archive
         - SQS Queue + FIFO Queue
         - Cross-region/account EventBus
         - HTTP requests (only enabled when MOTO_EVENTS_INVOKE_HTTP=true)
        """
        num_events = len(events)

        if num_events > 10:
            # the exact error text is longer, the Value list consists of all the put events
            raise ValidationException(
                "1 validation error detected: "
                "Value '[PutEventsRequestEntry]' at 'entries' failed to satisfy constraint: "
                "Member must have length less than or equal to 10"
            )

        entries = []
        for event in events:
            if not event.get("Source"):
                entries.append(
                    {
                        "ErrorCode": "InvalidArgument",
                        "ErrorMessage": "Parameter Source is not valid. Reason: Source is a required argument.",
                    }
                )
            elif not event.get("DetailType"):
                entries.append(
                    {
                        "ErrorCode": "InvalidArgument",
                        "ErrorMessage": "Parameter DetailType is not valid. Reason: DetailType is a required argument.",
                    }
                )
            elif not event.get("Detail"):
                entries.append(
                    {
                        "ErrorCode": "InvalidArgument",
                        "ErrorMessage": "Parameter Detail is not valid. Reason: Detail is a required argument.",
                    }
                )
            else:
                try:
                    detail = json.loads(event["Detail"])
                    if not isinstance(detail, dict):
                        warnings.warn(
                            f"EventDetail should be of type dict - types such as {type(detail)} are ignored by AWS"
                        )
                except ValueError:  # json.JSONDecodeError exists since Python 3.5
                    entries.append(
                        {
                            "ErrorCode": "MalformedDetail",
                            "ErrorMessage": "Detail is malformed.",
                        }
                    )
                    continue

                event_id = str(random.uuid4())
                entries.append({"EventId": event_id})

                # if 'EventBusName' is not especially set, it will be sent to the default one
                event_bus_name = self._normalize_event_bus_arn(
                    event.get("EventBusName")
                )

                event_bus = self.describe_event_bus(event_bus_name)
                for rule in event_bus.rules.values():
                    event_msg = copy.deepcopy(_BASE_EVENT_MESSAGE)
                    event_msg["id"] = event_id
                    event_msg["detail-type"] = event["DetailType"]
                    event_msg["source"] = event["Source"]
                    event_msg["account"] = self.account_id
                    event_msg["time"] = event.get("Time", unix_time())
                    event_msg["region"] = self.region_name
                    event_msg["resources"] = event.get("Resources", [])
                    event_msg["detail"] = json.loads(event["Detail"])
                    event_msg["ingestion-time"] = unix_time()
                    rule.send_to_targets(event_msg)

        return entries

    def remove_targets(
        self, name: str, event_bus_arn: Optional[str], ids: List[str]
    ) -> None:
        event_bus_name = self._normalize_event_bus_arn(event_bus_arn)
        event_bus = self._get_event_bus(event_bus_name)
        rule = event_bus.rules.get(name)

        if not rule:
            raise ResourceNotFoundException(
                f"Rule {name} does not exist on EventBus {event_bus_name}."
            )

        rule.remove_targets(ids)

    def test_event_pattern(self) -> None:
        raise NotImplementedError()

    @staticmethod
    def _put_permission_from_policy(event_bus: EventBus, policy: str) -> None:
        try:
            policy_doc = json.loads(policy)
            event_bus.add_policy(policy_doc)
        except JSONDecodeError:
            raise JsonRESTError(
                "ValidationException", "This policy contains invalid Json"
            )

    @staticmethod
    def _condition_param_to_stmt_condition(  # type: ignore[misc]
        condition: Optional[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        if condition:
            key = condition["Key"]
            value = condition["Value"]
            condition_type = condition["Type"]
            return {condition_type: {key: value}}
        return None

    def _put_permission_from_params(
        self,
        event_bus: EventBus,
        action: Optional[str],
        principal: str,
        statement_id: str,
        condition: Dict[str, str],
    ) -> None:
        if principal is None:
            raise JsonRESTError(
                "ValidationException", "Parameter Principal must be specified."
            )

        if condition and principal != "*":
            raise JsonRESTError(
                "InvalidParameterValue",
                "Value of the parameter 'principal' must be '*' when the parameter 'condition' is set.",
            )

        if not condition and self.ACCOUNT_ID.match(principal) is None:
            raise JsonRESTError(
                "InvalidParameterValue",
                f"Value {principal} at 'principal' failed to satisfy constraint: "
                r"Member must satisfy regular expression pattern: (\d{12}|\*)",
            )

        if action is None or action != "events:PutEvents":
            raise JsonRESTError(
                "ValidationException",
                "Provided value in parameter 'action' is not supported.",
            )

        if statement_id is None or self.STATEMENT_ID.match(statement_id) is None:
            raise JsonRESTError(
                "InvalidParameterValue", r"StatementId must match ^[a-zA-Z0-9-_]{1,64}$"
            )

        principal_arn = {
            "AWS": f"arn:{get_partition(self.region_name)}:iam::{principal}:root"
        }
        stmt_condition = self._condition_param_to_stmt_condition(condition)
        event_bus.add_permission(statement_id, action, principal_arn, stmt_condition)

    def put_permission(
        self,
        event_bus_name: str,
        action: str,
        principal: str,
        statement_id: str,
        condition: Dict[str, str],
        policy: str,
    ) -> None:
        if not event_bus_name:
            event_bus_name = "default"

        event_bus = self.describe_event_bus(event_bus_name)

        if policy:
            self._put_permission_from_policy(event_bus, policy)
        else:
            self._put_permission_from_params(
                event_bus, action, principal, statement_id, condition
            )

    def remove_permission(
        self,
        event_bus_name: Optional[str],
        statement_id: str,
        remove_all_permissions: bool,
    ) -> None:
        if not event_bus_name:
            event_bus_name = "default"

        event_bus = self.describe_event_bus(event_bus_name)

        if remove_all_permissions:
            event_bus.remove_statements()
        else:
            if not event_bus.has_permissions():
                raise JsonRESTError(
                    "ResourceNotFoundException", "EventBus does not have a policy."
                )

            statement = event_bus.remove_statement(statement_id)
            if not statement:
                raise JsonRESTError(
                    "ResourceNotFoundException",
                    "Statement with the provided id does not exist.",
                )

    def describe_event_bus(self, name: str) -> EventBus:
        if not name:
            name = "default"

        event_bus = self._get_event_bus(name)

        return event_bus

    def create_event_bus(
        self,
        name: str,
        event_source_name: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ) -> EventBus:
        if name in self.event_buses:
            raise JsonRESTError(
                "ResourceAlreadyExistsException", f"Event bus {name} already exists."
            )

        if not event_source_name and "/" in name:
            raise JsonRESTError(
                "ValidationException", "Event bus name must not contain '/'."
            )

        if event_source_name and event_source_name not in self.event_sources:
            raise JsonRESTError(
                "ResourceNotFoundException",
                f"Event source {event_source_name} does not exist.",
            )

        event_bus = EventBus(self.account_id, self.region_name, name, tags=tags)
        self.event_buses[name] = event_bus
        if tags:
            self.tagger.tag_resource(event_bus.arn, tags)

        return self.event_buses[name]

    def list_event_buses(self, name_prefix: Optional[str]) -> List[EventBus]:
        if name_prefix:
            return [
                event_bus
                for event_bus in self.event_buses.values()
                if event_bus.name.startswith(name_prefix)
            ]

        return list(self.event_buses.values())

    def delete_event_bus(self, name: str) -> None:
        if name == "default":
            raise JsonRESTError(
                "ValidationException", "Cannot delete event bus default."
            )
        event_bus = self.event_buses.pop(name, None)
        if event_bus:
            self.tagger.delete_all_tags_for_resource(event_bus.arn)

    def list_tags_for_resource(self, arn: str) -> Dict[str, List[Dict[str, str]]]:
        name = arn.split("/")[-1]
        rules = [bus.rules for bus in self.event_buses.values()]
        for registry in rules + [self.event_buses]:
            if name in registry:
                return self.tagger.list_tags_for_resource(registry[name].arn)
        raise ResourceNotFoundException(
            f"Rule {name} does not exist on EventBus default."
        )

    def tag_resource(self, arn: str, tags: List[Dict[str, str]]) -> None:
        name = arn.split("/")[-1]
        rules = [bus.rules for bus in self.event_buses.values()]
        for registry in rules + [self.event_buses]:
            if name in registry:
                self.tagger.tag_resource(registry[name].arn, tags)
                return
        raise ResourceNotFoundException(
            f"Rule {name} does not exist on EventBus default."
        )

    def untag_resource(self, arn: str, tag_names: List[str]) -> None:
        name = arn.split("/")[-1]
        rules = [bus.rules for bus in self.event_buses.values()]
        for registry in rules + [self.event_buses]:
            if name in registry:
                self.tagger.untag_resource_using_names(registry[name].arn, tag_names)
                return
        raise ResourceNotFoundException(
            f"Rule {name} does not exist on EventBus default."
        )

    def create_archive(
        self,
        name: str,
        source_arn: str,
        description: str,
        event_pattern: str,
        retention: str,
    ) -> Archive:
        if len(name) > 48:
            raise ValidationException(
                " 1 validation error detected: "
                f"Value '{name}' at 'archiveName' failed to satisfy constraint: "
                "Member must have length less than or equal to 48"
            )

        event_bus = self._get_event_bus(source_arn)

        if name in self.archives:
            raise ResourceAlreadyExistsException(f"Archive {name} already exists.")

        archive = Archive(
            self.account_id,
            self.region_name,
            name,
            source_arn,
            description,
            event_pattern,
            retention,
        )

        rule_event_pattern = json.loads(event_pattern or "{}")
        rule_event_pattern["replay-name"] = [{"exists": False}]

        rule_name = f"Events-Archive-{name}"
        rule = self.put_rule(
            rule_name,
            event_pattern=json.dumps(rule_event_pattern),
            event_bus_arn=event_bus.name,
            managed_by="prod.vhs.events.aws.internal",
        )
        archive.rule = rule
        self.put_targets(
            rule.name,
            rule.event_bus_name,
            [
                {
                    "Id": rule.name,
                    "Arn": f"arn:{get_partition(self.region_name)}:events:{self.region_name}:::",
                    "InputTransformer": {
                        "InputPathsMap": {},
                        # Template is not valid JSON, so we can't json.dump it
                        "InputTemplate": f'{{"archive-arn": "{archive.arn}:{archive.uuid}", "event": <aws.events.event.json>, "ingestion-time": <aws.events.event.ingestion-time>}}',
                    },
                }
            ],
        )

        self.archives[name] = archive

        return archive

    def describe_archive(self, name: str) -> Dict[str, Any]:
        archive = self.archives.get(name)

        if not archive:
            raise ResourceNotFoundException(f"Archive {name} does not exist.")

        return archive.describe()

    def list_archives(
        self,
        name_prefix: Optional[str],
        source_arn: Optional[str],
        state: Optional[str],
    ) -> List[Dict[str, Any]]:
        if [name_prefix, source_arn, state].count(None) < 2:
            raise ValidationException(
                "At most one filter is allowed for ListArchives. "
                "Use either : State, EventSourceArn, or NamePrefix."
            )

        if state and state not in Archive.VALID_STATES:
            valid_states = ", ".join(Archive.VALID_STATES)
            raise ValidationException(
                "1 validation error detected: "
                f"Value '{state}' at 'state' failed to satisfy constraint: "
                f"Member must satisfy enum value set: [{valid_states}]"
            )

        if [name_prefix, source_arn, state].count(None) == 3:
            return [archive.describe_short() for archive in self.archives.values()]

        result = []

        for archive in self.archives.values():
            if name_prefix and archive.name.startswith(name_prefix):
                result.append(archive.describe_short())
            elif source_arn and archive.source_arn == source_arn:
                result.append(archive.describe_short())
            elif state and archive.state == state:
                result.append(archive.describe_short())

        return result

    def update_archive(
        self, name: str, description: str, event_pattern: str, retention: str
    ) -> Dict[str, Any]:
        archive = self.archives.get(name)

        if not archive:
            raise ResourceNotFoundException(f"Archive {name} does not exist.")

        archive.update(description, event_pattern, retention)

        return {
            "ArchiveArn": archive.arn,
            "CreationTime": archive.creation_time,
            "State": archive.state,
        }

    def delete_archive(self, name: str) -> None:
        archive = self.archives.get(name)

        if not archive:
            raise ResourceNotFoundException(f"Archive {name} does not exist.")

        archive.delete(self.account_id, self.region_name)

    def start_replay(
        self,
        name: str,
        description: str,
        source_arn: str,
        start_time: str,
        end_time: str,
        destination: Dict[str, Any],
    ) -> Dict[str, Any]:
        event_bus_arn = destination["Arn"]
        event_bus_arn_pattern = (
            rf"{ARN_PARTITION_REGEX}:events:[a-zA-Z0-9-]+:\d{{12}}:event-bus/"
        )
        if not re.match(event_bus_arn_pattern, event_bus_arn):
            raise ValidationException(
                "Parameter Destination.Arn is not valid. Reason: Must contain an event bus ARN."
            )

        self._get_event_bus(event_bus_arn)

        archive_name = source_arn.split("/")[-1]
        archive = self.archives.get(archive_name)
        if not archive:
            raise ValidationException(
                f"Parameter EventSourceArn is not valid. Reason: Archive {archive_name} does not exist."
            )

        if event_bus_arn != archive.source_arn:
            raise ValidationException(
                "Parameter Destination.Arn is not valid. "
                "Reason: Cross event bus replay is not permitted."
            )

        if start_time > end_time:
            raise ValidationException(
                "Parameter EventEndTime is not valid. "
                "Reason: EventStartTime must be before EventEndTime."
            )

        if name in self.replays:
            raise ResourceAlreadyExistsException(f"Replay {name} already exists.")

        replay = Replay(
            self.account_id,
            self.region_name,
            name,
            description,
            source_arn,
            start_time,
            end_time,
            destination,
        )

        self.replays[name] = replay

        replay.replay_events(archive)

        return {
            "ReplayArn": replay.arn,
            "ReplayStartTime": replay.start_time,
            "State": ReplayState.STARTING.value,  # the replay will be done before returning the response
        }

    def describe_replay(self, name: str) -> Dict[str, Any]:
        replay = self._get_replay(name)

        return replay.describe()

    def list_replays(
        self, name_prefix: str, source_arn: str, state: str
    ) -> List[Dict[str, Any]]:
        if [name_prefix, source_arn, state].count(None) < 2:  # type: ignore
            raise ValidationException(
                "At most one filter is allowed for ListReplays. "
                "Use either : State, EventSourceArn, or NamePrefix."
            )

        valid_states = sorted([item.value for item in ReplayState])
        if state and state not in valid_states:
            all_states = ", ".join(valid_states)
            raise ValidationException(
                f"1 validation error detected: Value '{state}' at 'state' failed to satisfy constraint: Member must satisfy enum value set: [{all_states}]"
            )

        if [name_prefix, source_arn, state].count(None) == 3:  # type: ignore
            return [replay.describe_short() for replay in self.replays.values()]

        result = []

        for replay in self.replays.values():
            if name_prefix and replay.name.startswith(name_prefix):
                result.append(replay.describe_short())
            elif source_arn and replay.source_arn == source_arn:
                result.append(replay.describe_short())
            elif state and replay.state == state:  # type: ignore
                result.append(replay.describe_short())

        return result

    def cancel_replay(self, name: str) -> Dict[str, str]:
        replay = self._get_replay(name)

        # replays in the state 'COMPLETED' can't be canceled,
        # but the implementation is done synchronously,
        # so they are done right after the start
        if replay.state not in [
            ReplayState.STARTING,
            ReplayState.RUNNING,
            ReplayState.COMPLETED,
        ]:
            raise IllegalStatusException(
                f"Replay {name} is not in a valid state for this operation."
            )

        replay.state = ReplayState.CANCELLED

        return {"ReplayArn": replay.arn, "State": ReplayState.CANCELLING.value}

    def create_connection(
        self,
        name: str,
        description: str,
        authorization_type: str,
        auth_parameters: Dict[str, Any],
    ) -> Connection:
        connection = Connection(
            name,
            self.account_id,
            self.region_name,
            description,
            authorization_type,
            auth_parameters,
        )
        self.connections[name] = connection
        return connection

    def update_connection(self, name: str, **kwargs: Any) -> Dict[str, Any]:
        connection = self.connections.get(name)
        if not connection:
            raise ResourceNotFoundException(f"Connection '{name}' does not exist.")

        for attr, value in kwargs.items():
            if value is not None and hasattr(connection, attr):
                setattr(connection, attr, value)
        return connection.describe_short()

    def list_connections(self) -> List[Connection]:
        return list(self.connections.values())

    def describe_connection(self, name: str) -> Dict[str, Any]:
        connection = self.connections.get(name)
        if not connection:
            raise ResourceNotFoundException(f"Connection '{name}' does not exist.")

        return connection.describe()

    def delete_connection(self, name: str) -> Dict[str, Any]:
        connection = self.connections.pop(name, None)
        if not connection:
            raise ResourceNotFoundException(f"Connection '{name}' does not exist.")

        # Delete Secret
        secretsmanager_backend = get_secrets_manager_backend(
            self.account_id, self.region_name
        )
        secretsmanager_backend.delete_secret(
            secret_id=connection.secret_arn,
            recovery_window_in_days=None,
            force_delete_without_recovery=True,
        )

        return connection.describe_short()

    def create_api_destination(
        self,
        name: str,
        description: str,
        connection_arn: str,
        invocation_endpoint: str,
        invocation_rate_limit_per_second: str,
        http_method: str,
    ) -> Dict[str, Any]:
        destination = Destination(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            description=description,
            connection_arn=connection_arn,
            invocation_endpoint=invocation_endpoint,
            invocation_rate_limit_per_second=invocation_rate_limit_per_second,
            http_method=http_method,
        )

        self.destinations[name] = destination
        return destination.describe_short()

    def list_api_destinations(self) -> List[Destination]:
        return list(self.destinations.values())

    def describe_api_destination(self, name: str) -> Dict[str, Any]:
        destination = self.destinations.get(name)
        if not destination:
            raise ResourceNotFoundException(
                f"An api-destination '{name}' does not exist."
            )
        return destination.describe()

    def update_api_destination(self, name: str, **kwargs: Any) -> Dict[str, Any]:
        destination = self.destinations.get(name)
        if not destination:
            raise ResourceNotFoundException(
                f"An api-destination '{name}' does not exist."
            )

        for attr, value in kwargs.items():
            if value is not None and hasattr(destination, attr):
                setattr(destination, attr, value)
        return destination.describe_short()

    def delete_api_destination(self, name: str) -> None:
        destination = self.destinations.pop(name, None)
        if not destination:
            raise ResourceNotFoundException(
                f"An api-destination '{name}' does not exist."
            )

    def create_partner_event_source(self, name: str, account_id: str) -> None:
        # https://docs.aws.amazon.com/eventbridge/latest/onboarding/amazon_eventbridge_partner_onboarding_guide.html
        if name not in self.partner_event_sources:
            self.partner_event_sources[name] = PartnerEventSource(
                region=self.region_name, name=name
            )
        self.partner_event_sources[name].accounts.append(account_id)
        client_backend = events_backends[account_id][self.region_name]
        client_backend.event_sources[name] = self.partner_event_sources[name]

    def describe_event_source(self, name: str) -> PartnerEventSource:
        return self.event_sources[name]

    def describe_partner_event_source(self, name: str) -> PartnerEventSource:
        return self.partner_event_sources[name]

    def delete_partner_event_source(self, name: str, account_id: str) -> None:
        client_backend = events_backends[account_id][self.region_name]
        client_backend.event_sources[name].state = "DELETED"

    def put_partner_events(self, entries: List[Dict[str, Any]]) -> None:
        """
        Validation of the entries is not yet implemented.
        """
        # This implementation is very basic currently, just to verify the behaviour
        # In the future we could create a batch of events, grouped by source, and send them all at once
        for entry in entries:
            source = entry["Source"]
            for account_id in self.partner_event_sources[source].accounts:
                events_backends[account_id][self.region_name].put_events([entry])


events_backends = BackendDict(EventsBackend, "events")
