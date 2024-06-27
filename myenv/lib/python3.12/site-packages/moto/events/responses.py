import json
from typing import Any, Dict, Tuple

from moto.core.responses import BaseResponse
from moto.events.models import EventsBackend, events_backends


class EventsHandler(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="events")

    @property
    def events_backend(self) -> EventsBackend:
        return events_backends[self.current_account][self.region]

    def _create_response(self, result: Any) -> Tuple[str, Dict[str, Any]]:
        """
        Creates a proper response for the API.

        It basically transforms a dict-like result from the backend
        into a tuple (str, dict) properly formatted.
        Args:
            result (dict): result from backend

        Returns:
            (str, dict): dumped result and headers
        """
        return json.dumps(result), self.response_headers

    def error(
        self, type_: str, message: str = "", status: int = 400
    ) -> Tuple[str, Dict[str, Any]]:
        headers: Dict[str, Any] = self.response_headers
        headers["status"] = status
        return json.dumps({"__type": type_, "message": message}), headers

    def put_rule(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        event_pattern = self._get_param("EventPattern")
        scheduled_expression = self._get_param("ScheduleExpression")
        state = self._get_param("State")
        desc = self._get_param("Description")
        role_arn = self._get_param("RoleArn")
        event_bus_arn = self._get_param("EventBusName")
        tags = self._get_param("Tags")

        rule = self.events_backend.put_rule(
            name,
            scheduled_expression=scheduled_expression,
            event_pattern=event_pattern,
            state=state,
            description=desc,
            role_arn=role_arn,
            event_bus_arn=event_bus_arn,
            tags=tags,
        )
        result = {"RuleArn": rule.arn}
        return self._create_response(result)

    def delete_rule(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        event_bus_arn = self._get_param("EventBusName")

        if not name:
            return self.error("ValidationException", "Parameter Name is required.")
        self.events_backend.delete_rule(name, event_bus_arn)

        return "", self.response_headers

    def describe_rule(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        event_bus_arn = self._get_param("EventBusName")

        if not name:
            return self.error("ValidationException", "Parameter Name is required.")

        rule = self.events_backend.describe_rule(name, event_bus_arn)

        result = rule.describe()
        return self._create_response(result)

    def disable_rule(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        event_bus_arn = self._get_param("EventBusName")

        if not name:
            return self.error("ValidationException", "Parameter Name is required.")

        if not self.events_backend.disable_rule(name, event_bus_arn):
            return self.error(
                "ResourceNotFoundException", "Rule " + name + " does not exist."
            )

        return "", self.response_headers

    def enable_rule(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        event_bus_arn = self._get_param("EventBusName")

        if not name:
            return self.error("ValidationException", "Parameter Name is required.")

        if not self.events_backend.enable_rule(name, event_bus_arn):
            return self.error(
                "ResourceNotFoundException", "Rule " + name + " does not exist."
            )

        return "", self.response_headers

    def generate_presigned_url(self) -> None:
        pass

    def list_rule_names_by_target(self) -> Tuple[str, Dict[str, Any]]:
        target_arn = self._get_param("TargetArn")
        event_bus_arn = self._get_param("EventBusName")
        next_token = self._get_param("NextToken")
        limit = self._get_param("Limit")

        if not target_arn:
            return self.error("ValidationException", "Parameter TargetArn is required.")

        rules, token = self.events_backend.list_rule_names_by_target(
            target_arn=target_arn,
            event_bus_arn=event_bus_arn,
            next_token=next_token,
            limit=limit,
        )

        res = {"RuleNames": [rule.name for rule in rules], "NextToken": token}

        return json.dumps(res), self.response_headers

    def list_rules(self) -> Tuple[str, Dict[str, Any]]:
        prefix = self._get_param("NamePrefix")
        event_bus_arn = self._get_param("EventBusName")
        next_token = self._get_param("NextToken")
        limit = self._get_param("Limit")

        rules, token = self.events_backend.list_rules(
            prefix=prefix,
            event_bus_arn=event_bus_arn,
            next_token=next_token,
            limit=limit,
        )
        rules_obj = {
            "Rules": [rule.describe() for rule in rules],
            "NextToken": token,
        }

        return json.dumps(rules_obj), self.response_headers

    def list_targets_by_rule(self) -> Tuple[str, Dict[str, Any]]:
        rule_name = self._get_param("Rule")
        event_bus_arn = self._get_param("EventBusName")
        next_token = self._get_param("NextToken")
        limit = self._get_param("Limit")

        if not rule_name:
            return self.error("ValidationException", "Parameter Rule is required.")

        try:
            targets = self.events_backend.list_targets_by_rule(
                rule_name, event_bus_arn, next_token, limit
            )
        except KeyError:
            return self.error(
                "ResourceNotFoundException", "Rule " + rule_name + " does not exist."
            )

        return json.dumps(targets), self.response_headers

    def put_events(self) -> str:
        events = self._get_param("Entries")

        entries = self.events_backend.put_events(events)

        failed_count = len([e for e in entries if "ErrorCode" in e])
        response = {
            "FailedEntryCount": failed_count,
            "Entries": entries,
        }

        return json.dumps(response)

    def put_targets(self) -> Tuple[str, Dict[str, Any]]:
        rule_name = self._get_param("Rule")
        event_bus_name = self._get_param("EventBusName")
        targets = self._get_param("Targets")

        self.events_backend.put_targets(rule_name, event_bus_name, targets)

        return (
            json.dumps({"FailedEntryCount": 0, "FailedEntries": []}),
            self.response_headers,
        )

    def remove_targets(self) -> Tuple[str, Dict[str, Any]]:
        rule_name = self._get_param("Rule")
        event_bus_name = self._get_param("EventBusName")
        ids = self._get_param("Ids")

        self.events_backend.remove_targets(rule_name, event_bus_name, ids)

        return (
            json.dumps({"FailedEntryCount": 0, "FailedEntries": []}),
            self.response_headers,
        )

    def test_event_pattern(self) -> None:
        pass

    def put_permission(self) -> str:
        event_bus_name = self._get_param("EventBusName")
        action = self._get_param("Action")
        principal = self._get_param("Principal")
        statement_id = self._get_param("StatementId")
        policy = self._get_param("Policy")
        condition = self._get_param("Condition")

        self.events_backend.put_permission(
            event_bus_name, action, principal, statement_id, condition, policy
        )

        return ""

    def remove_permission(self) -> str:
        event_bus_name = self._get_param("EventBusName")
        statement_id = self._get_param("StatementId")
        remove_all_permissions = self._get_param("RemoveAllPermissions")

        self.events_backend.remove_permission(
            event_bus_name, statement_id, remove_all_permissions
        )

        return ""

    def describe_event_bus(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")

        event_bus = self.events_backend.describe_event_bus(name)
        response = {"Name": event_bus.name, "Arn": event_bus.arn}

        if event_bus.policy:
            response["Policy"] = event_bus.policy

        return json.dumps(response), self.response_headers

    def create_event_bus(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        event_source_name = self._get_param("EventSourceName")
        tags = self._get_param("Tags")

        event_bus = self.events_backend.create_event_bus(name, event_source_name, tags)
        return json.dumps({"EventBusArn": event_bus.arn}), self.response_headers

    def list_event_buses(self) -> Tuple[str, Dict[str, Any]]:
        name_prefix = self._get_param("NamePrefix")
        # ToDo: add 'NextToken' & 'Limit' parameters

        response = []
        for event_bus in self.events_backend.list_event_buses(name_prefix):
            event_bus_response = {"Name": event_bus.name, "Arn": event_bus.arn}

            if event_bus.policy:
                event_bus_response["Policy"] = event_bus.policy

            response.append(event_bus_response)

        return json.dumps({"EventBuses": response}), self.response_headers

    def delete_event_bus(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")

        self.events_backend.delete_event_bus(name)

        return "", self.response_headers

    def list_tags_for_resource(self) -> Tuple[str, Dict[str, Any]]:
        arn = self._get_param("ResourceARN")

        result = self.events_backend.list_tags_for_resource(arn)

        return json.dumps(result), self.response_headers

    def tag_resource(self) -> Tuple[str, Dict[str, Any]]:
        arn = self._get_param("ResourceARN")
        tags = self._get_param("Tags")

        self.events_backend.tag_resource(arn, tags)

        return "{}", self.response_headers

    def untag_resource(self) -> Tuple[str, Dict[str, Any]]:
        arn = self._get_param("ResourceARN")
        tags = self._get_param("TagKeys")

        self.events_backend.untag_resource(arn, tags)

        return "{}", self.response_headers

    def create_archive(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ArchiveName")
        source_arn = self._get_param("EventSourceArn")
        description = self._get_param("Description")
        event_pattern = self._get_param("EventPattern")
        retention = self._get_param("RetentionDays")

        archive = self.events_backend.create_archive(
            name, source_arn, description, event_pattern, retention
        )

        return (
            json.dumps(
                {
                    "ArchiveArn": archive.arn,
                    "CreationTime": archive.creation_time,
                    "State": archive.state,
                }
            ),
            self.response_headers,
        )

    def describe_archive(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ArchiveName")

        result = self.events_backend.describe_archive(name)

        return json.dumps(result), self.response_headers

    def list_archives(self) -> Tuple[str, Dict[str, Any]]:
        name_prefix = self._get_param("NamePrefix")
        source_arn = self._get_param("EventSourceArn")
        state = self._get_param("State")

        result = self.events_backend.list_archives(name_prefix, source_arn, state)

        return json.dumps({"Archives": result}), self.response_headers

    def update_archive(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ArchiveName")
        description = self._get_param("Description")
        event_pattern = self._get_param("EventPattern")
        retention = self._get_param("RetentionDays")

        result = self.events_backend.update_archive(
            name, description, event_pattern, retention
        )

        return json.dumps(result), self.response_headers

    def delete_archive(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ArchiveName")

        self.events_backend.delete_archive(name)

        return "", self.response_headers

    def start_replay(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ReplayName")
        description = self._get_param("Description")
        source_arn = self._get_param("EventSourceArn")
        start_time = self._get_param("EventStartTime")
        end_time = self._get_param("EventEndTime")
        destination = self._get_param("Destination")

        result = self.events_backend.start_replay(
            name, description, source_arn, start_time, end_time, destination
        )

        return json.dumps(result), self.response_headers

    def describe_replay(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ReplayName")

        result = self.events_backend.describe_replay(name)

        return json.dumps(result), self.response_headers

    def list_replays(self) -> Tuple[str, Dict[str, Any]]:
        name_prefix = self._get_param("NamePrefix")
        source_arn = self._get_param("EventSourceArn")
        state = self._get_param("State")

        result = self.events_backend.list_replays(name_prefix, source_arn, state)

        return json.dumps({"Replays": result}), self.response_headers

    def cancel_replay(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("ReplayName")

        result = self.events_backend.cancel_replay(name)

        return json.dumps(result), self.response_headers

    def create_connection(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        description = self._get_param("Description")
        authorization_type = self._get_param("AuthorizationType")
        auth_parameters = self._get_param("AuthParameters")

        result = self.events_backend.create_connection(
            name, description, authorization_type, auth_parameters
        )

        return (
            json.dumps(
                {
                    "ConnectionArn": result.arn,
                    "ConnectionState": "AUTHORIZED",
                    "CreationTime": result.creation_time,
                    "LastModifiedTime": result.creation_time,
                }
            ),
            self.response_headers,
        )

    def list_connections(self) -> Tuple[str, Dict[str, Any]]:
        connections = self.events_backend.list_connections()
        result = []
        for connection in connections:
            result.append(
                {
                    "ConnectionArn": connection.arn,
                    "ConnectionState": "AUTHORIZED",
                    "CreationTime": connection.creation_time,
                    "LastModifiedTime": connection.creation_time,
                    "AuthorizationType": connection.authorization_type,
                }
            )

        return json.dumps({"Connections": result}), self.response_headers

    def describe_connection(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        result = self.events_backend.describe_connection(name)
        return json.dumps(result), self.response_headers

    def update_connection(self) -> Tuple[str, Dict[str, Any]]:
        updates = dict(
            name=self._get_param("Name"),
            description=self._get_param("Description"),
            authorization_type=self._get_param("AuthorizationType"),
            auth_parameters=self._get_param("AuthParameters"),
        )
        result = self.events_backend.update_connection(**updates)
        return self._create_response(result)

    def delete_connection(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        result = self.events_backend.delete_connection(name)
        return json.dumps(result), self.response_headers

    def create_api_destination(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        description = self._get_param("Description")
        connection_arn = self._get_param("ConnectionArn")
        invocation_endpoint = self._get_param("InvocationEndpoint")
        invocation_rate_limit_per_second = self._get_param(
            "InvocationRateLimitPerSecond"
        )
        http_method = self._get_param("HttpMethod")

        result = self.events_backend.create_api_destination(
            name,
            description,
            connection_arn,
            invocation_endpoint,
            invocation_rate_limit_per_second,
            http_method,
        )
        return self._create_response(result)

    def list_api_destinations(self) -> Tuple[str, Dict[str, Any]]:
        destinations = self.events_backend.list_api_destinations()
        result = []
        for destination in destinations:
            result.append(
                {
                    "ApiDestinationArn": destination.arn,
                    "Name": destination.name,
                    "ApiDestinationState": destination.state,
                    "ConnectionArn": destination.connection_arn,
                    "InvocationEndpoint": destination.invocation_endpoint,
                    "HttpMethod": destination.http_method,
                    "CreationTime": destination.creation_time,
                    "LastModifiedTime": destination.creation_time,
                }
            )

        return json.dumps({"ApiDestinations": result}), self.response_headers

    def describe_api_destination(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        result = self.events_backend.describe_api_destination(name)
        return self._create_response(result)

    def update_api_destination(self) -> Tuple[str, Dict[str, Any]]:
        updates = dict(
            connection_arn=self._get_param("ConnectionArn"),
            description=self._get_param("Description"),
            http_method=self._get_param("HttpMethod"),
            invocation_endpoint=self._get_param("InvocationEndpoint"),
            invocation_rate_limit_per_second=self._get_param(
                "InvocationRateLimitPerSecond"
            ),
            name=self._get_param("Name"),
        )

        result = self.events_backend.update_api_destination(**updates)
        return self._create_response(result)

    def delete_api_destination(self) -> Tuple[str, Dict[str, Any]]:
        name = self._get_param("Name")
        self.events_backend.delete_api_destination(name)
        return self._create_response({})

    def create_partner_event_source(self) -> str:
        name = self._get_param("Name")
        account_id = self._get_param("Account")
        self.events_backend.create_partner_event_source(
            name=name,
            account_id=account_id,
        )
        return "{}"

    def describe_event_source(self) -> str:
        name = self._get_param("Name")
        event_source = self.events_backend.describe_event_source(name)
        return json.dumps(event_source.to_dict())

    def describe_partner_event_source(self) -> str:
        name = self._get_param("Name")
        event_source = self.events_backend.describe_partner_event_source(name)
        return json.dumps({"Arn": event_source.arn, "Name": event_source.name})

    def delete_partner_event_source(self) -> str:
        name = self._get_param("Name")
        account_id = self._get_param("Account")
        self.events_backend.delete_partner_event_source(name, account_id)
        return "{}"

    def put_partner_events(self) -> str:
        entries = self._get_param("Entries")
        self.events_backend.put_partner_events(entries)
        return json.dumps({"Entries": [], "FailedEntryCount": 0})
