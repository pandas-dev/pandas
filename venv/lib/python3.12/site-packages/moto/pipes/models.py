"""EventBridgePipesBackend class with methods for supported APIs."""

from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_without_milliseconds, utcnow
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import NotFoundException

PAGINATION_MODEL = {
    "list_pipes": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "Arn",
    },
}


class Pipe(BaseModel):
    """Represents an EventBridge Pipe."""

    def __init__(
        self,
        name: str,
        account_id: str,
        region_name: str,
        source: str,
        target: str,
        role_arn: str,
        description: Optional[str] = None,
        desired_state: Optional[str] = None,
        source_parameters: Optional[dict[str, Any]] = None,
        enrichment: Optional[str] = None,
        enrichment_parameters: Optional[dict[str, Any]] = None,
        target_parameters: Optional[dict[str, Any]] = None,
        tags: Optional[dict[str, str]] = None,
        log_configuration: Optional[dict[str, Any]] = None,
        kms_key_identifier: Optional[str] = None,
    ):
        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        self.description = description
        self.desired_state = desired_state or "RUNNING"
        self.current_state = "RUNNING"
        self.source = source
        self.source_parameters = source_parameters
        self.enrichment = enrichment
        self.enrichment_parameters = enrichment_parameters
        self.target = target
        self.target_parameters = target_parameters
        self.role_arn = role_arn
        self.tags = tags or {}
        self.log_configuration = log_configuration
        self.kms_key_identifier = kms_key_identifier
        self.creation_time = utcnow()
        self.last_modified_time = utcnow()
        self.state_reason = None

    @property
    def arn(self) -> str:
        partition = get_partition(self.region_name)
        return f"arn:{partition}:pipes:{self.region_name}:{self.account_id}:pipe/{self.name}"


class EventBridgePipesBackend(BaseBackend):
    """Implementation of EventBridgePipes APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.pipes: dict[str, Pipe] = {}
        self.tagger = TaggingService()

    def create_pipe(
        self,
        name: str,
        description: Optional[str],
        desired_state: Optional[str],
        source: str,
        source_parameters: Optional[dict[str, Any]],
        enrichment: Optional[str],
        enrichment_parameters: Optional[dict[str, Any]],
        target: str,
        target_parameters: Optional[dict[str, Any]],
        role_arn: str,
        tags: Optional[dict[str, str]],
        log_configuration: Optional[dict[str, Any]],
        kms_key_identifier: Optional[str],
    ) -> Pipe:
        pipe = Pipe(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
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

        self.pipes[name] = pipe

        return pipe

    def describe_pipe(self, name: str) -> Pipe:
        if name not in self.pipes:
            raise NotFoundException(f"Pipe {name} not found")
        return self.pipes[name]

    def delete_pipe(self, name: str) -> tuple[str, str, str, str, str, str]:
        if name not in self.pipes:
            raise NotFoundException(f"Pipe {name} not found")
        pipe = self.pipes[name]
        pipe.desired_state = "DELETED"
        pipe.current_state = "DELETING"

        arn = pipe.arn
        pipe_name = pipe.name
        desired_state = pipe.desired_state
        current_state = pipe.current_state
        creation_time = iso_8601_datetime_without_milliseconds(pipe.creation_time)
        last_modified_time = iso_8601_datetime_without_milliseconds(
            pipe.last_modified_time
        )

        del self.pipes[name]

        return (
            arn,
            pipe_name,
            desired_state,
            current_state,
            creation_time,
            last_modified_time,
        )

    def tag_resource(self, resource_arn: str, tags: dict[str, str]) -> None:
        pipe = None
        for p in self.pipes.values():
            if p.arn == resource_arn:
                pipe = p
                break

        if pipe is None:
            raise NotFoundException(f"Resource {resource_arn} not found")

        tag_list = TaggingService.convert_dict_to_tags_input(tags)
        self.tagger.tag_resource(resource_arn, tag_list)

        pipe.tags.update(tags)

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        pipe = None
        for p in self.pipes.values():
            if p.arn == resource_arn:
                pipe = p
                break

        if pipe is None:
            raise NotFoundException(f"Resource {resource_arn} not found")

        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

        for tag_key in tag_keys:
            pipe.tags.pop(tag_key, None)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_pipes(  # type: ignore[misc]
        self,
        name_prefix: Optional[str],
        desired_state: Optional[str],
        current_state: Optional[str],
        source_prefix: Optional[str],
        target_prefix: Optional[str],
    ) -> list[dict[str, Any]]:
        """List pipes with optional filtering and pagination."""
        filtered_pipes = list(self.pipes.values())

        if name_prefix:
            filtered_pipes = [
                p for p in filtered_pipes if p.name.startswith(name_prefix)
            ]

        if desired_state:
            filtered_pipes = [
                p for p in filtered_pipes if p.desired_state == desired_state
            ]

        if current_state:
            filtered_pipes = [
                p for p in filtered_pipes if p.current_state == current_state
            ]

        if source_prefix:
            filtered_pipes = [
                p
                for p in filtered_pipes
                if p.source and p.source.startswith(source_prefix)
            ]

        if target_prefix:
            filtered_pipes = [
                p
                for p in filtered_pipes
                if p.target and p.target.startswith(target_prefix)
            ]

        filtered_pipes.sort(key=lambda p: p.name)

        pipe_summaries = []
        for pipe in filtered_pipes:
            summary: dict[str, Any] = {
                "Arn": pipe.arn,
                "Name": pipe.name,
                "DesiredState": pipe.desired_state,
                "CurrentState": pipe.current_state,
                "Source": pipe.source,
                "Target": pipe.target,
                "CreationTime": iso_8601_datetime_without_milliseconds(
                    pipe.creation_time
                ),
                "LastModifiedTime": iso_8601_datetime_without_milliseconds(
                    pipe.last_modified_time
                ),
            }
            if pipe.enrichment:
                summary["Enrichment"] = pipe.enrichment
            pipe_summaries.append(summary)

        return pipe_summaries


pipes_backends = BackendDict(EventBridgePipesBackend, "pipes")
