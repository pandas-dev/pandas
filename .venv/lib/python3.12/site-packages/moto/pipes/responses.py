"""Handles incoming pipes requests, invokes methods, returns responses."""

import json
from typing import Any
from urllib.parse import unquote

from moto.core.responses import BaseResponse
from moto.core.utils import iso_8601_datetime_without_milliseconds

from .exceptions import ValidationException
from .models import EventBridgePipesBackend, pipes_backends


class EventBridgePipesResponse(BaseResponse):
    """Handler for EventBridgePipes requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="pipes")

    @property
    def pipes_backend(self) -> EventBridgePipesBackend:
        return pipes_backends[self.current_account][self.region]

    def create_pipe(self) -> str:
        body_params = json.loads(self.body) if self.body else {}

        name = body_params.get("Name") or self.uri.split("/")[-1]
        description = body_params.get("Description")
        desired_state = body_params.get("DesiredState")
        source = body_params.get("Source")
        source_parameters = body_params.get("SourceParameters")
        enrichment = body_params.get("Enrichment")
        enrichment_parameters = body_params.get("EnrichmentParameters")
        target = body_params.get("Target")
        target_parameters = body_params.get("TargetParameters")
        role_arn = body_params.get("RoleArn")
        tags = body_params.get("Tags")
        log_configuration = body_params.get("LogConfiguration")
        kms_key_identifier = body_params.get("KmsKeyIdentifier")

        if not source:
            raise ValidationException("Source is a required parameter")
        if not target:
            raise ValidationException("Target is a required parameter")
        if not role_arn:
            raise ValidationException("RoleArn is a required parameter")

        pipe = self.pipes_backend.create_pipe(
            name=name,
            description=description,
            desired_state=desired_state,
            source=source,
            source_parameters=source_parameters,
            enrichment=enrichment,
            enrichment_parameters=enrichment_parameters,
            target=target,
            target_parameters=target_parameters,
            role_arn=role_arn,
            tags=tags,
            log_configuration=log_configuration,
            kms_key_identifier=kms_key_identifier,
        )
        return json.dumps(
            {
                "Arn": pipe.arn,
                "Name": pipe.name,
                "DesiredState": pipe.desired_state,
                "CurrentState": pipe.current_state,
                "CreationTime": iso_8601_datetime_without_milliseconds(
                    pipe.creation_time
                ),
                "LastModifiedTime": iso_8601_datetime_without_milliseconds(
                    pipe.last_modified_time
                ),
            }
        )

    def describe_pipe(self) -> str:
        name = self.uri.split("?")[0].split("/")[-1]
        pipe = self.pipes_backend.describe_pipe(
            name=name,
        )
        response_dict: dict[str, Any] = {
            "Arn": pipe.arn,
            "Name": pipe.name,
            "DesiredState": pipe.desired_state,
            "CurrentState": pipe.current_state,
            "CreationTime": iso_8601_datetime_without_milliseconds(pipe.creation_time),
            "LastModifiedTime": iso_8601_datetime_without_milliseconds(
                pipe.last_modified_time
            ),
        }
        if pipe.description is not None:
            response_dict["Description"] = pipe.description
        if pipe.state_reason is not None:
            response_dict["StateReason"] = pipe.state_reason
        if pipe.source is not None:
            response_dict["Source"] = pipe.source
        if pipe.source_parameters is not None:
            response_dict["SourceParameters"] = pipe.source_parameters
        if pipe.enrichment is not None:
            response_dict["Enrichment"] = pipe.enrichment
        if pipe.enrichment_parameters is not None:
            response_dict["EnrichmentParameters"] = pipe.enrichment_parameters
        if pipe.target is not None:
            response_dict["Target"] = pipe.target
        if pipe.target_parameters is not None:
            response_dict["TargetParameters"] = pipe.target_parameters
        if pipe.role_arn is not None:
            response_dict["RoleArn"] = pipe.role_arn
        if pipe.tags is not None:
            response_dict["Tags"] = pipe.tags
        if pipe.log_configuration is not None:
            response_dict["LogConfiguration"] = pipe.log_configuration
        if pipe.kms_key_identifier is not None:
            response_dict["KmsKeyIdentifier"] = pipe.kms_key_identifier
        return json.dumps(response_dict)

    def delete_pipe(self) -> str:
        name = self.uri.split("?")[0].split("/")[-1]

        arn, name, desired_state, current_state, creation_time, last_modified_time = (
            self.pipes_backend.delete_pipe(
                name=name,
            )
        )

        return json.dumps(
            {
                "Arn": arn,
                "Name": name,
                "DesiredState": desired_state,
                "CurrentState": current_state,
                "CreationTime": creation_time,
                "LastModifiedTime": last_modified_time,
            }
        )

    def tag_resource(self) -> str:
        resource_arn = unquote(self.uri.split("/tags/")[-1])
        body_params = json.loads(self.body) if self.body else {}
        tags = body_params.get("Tags") or body_params.get("tags")
        if not tags:
            raise ValidationException("Tags is a required parameter")

        self.pipes_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps({})

    def untag_resource(self) -> str:
        resource_arn = unquote(self.uri.split("?")[0].split("/tags/")[-1])
        tag_keys = self.querystring.get("tagKeys", [])

        self.pipes_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return json.dumps({})

    def list_pipes(self) -> str:
        params = json.loads(self.body) if self.body else {}
        if not params and self.querystring:
            params = {
                k: (v[0] if isinstance(v, list) and len(v) > 0 else v)
                for k, v in self.querystring.items()
            }

        name_prefix = params.get("NamePrefix") or params.get("namePrefix")
        desired_state = params.get("DesiredState") or params.get("desiredState")
        current_state = params.get("CurrentState") or params.get("currentState")
        source_prefix = params.get("SourcePrefix") or params.get("sourcePrefix")
        target_prefix = params.get("TargetPrefix") or params.get("targetPrefix")
        next_token = params.get("NextToken") or params.get("nextToken")
        limit = params.get("Limit") or params.get("limit")
        if limit is not None and limit != "":
            try:
                limit = int(limit)
            except (ValueError, TypeError):
                limit = None

        pipes, next_token = self.pipes_backend.list_pipes(
            name_prefix=name_prefix,
            desired_state=desired_state,
            current_state=current_state,
            source_prefix=source_prefix,
            target_prefix=target_prefix,
            next_token=next_token,
            limit=limit,
        )

        response_dict: dict[str, Any] = {"Pipes": pipes}
        if next_token:
            response_dict["NextToken"] = next_token

        return json.dumps(response_dict)
