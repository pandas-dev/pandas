"""EventBridgePipesBackend class with methods for supported APIs."""

from enum import Enum
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.moto_api._internal.managed_state_model import ManagedState
from moto.utilities.paginator import paginate
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


class PipeStatus(str, Enum):
    CREATING = "CREATING"
    RUNNING = "RUNNING"
    STARTING = "STARTING"
    STOPPING = "STOPPING"
    STOPPED = "STOPPED"
    DELETING = "DELETING"
    DELETED = "DELETED"

    @classmethod
    def status_transitions(self) -> list[tuple[str | None, str]]:
        return [
            (PipeStatus.CREATING.value, PipeStatus.RUNNING.value),
            (PipeStatus.STARTING.value, PipeStatus.RUNNING.value),
            (PipeStatus.STOPPING.value, PipeStatus.STOPPED),
            (PipeStatus.STOPPED.value, PipeStatus.STARTING.value),
            (PipeStatus.DELETING.value, PipeStatus.DELETED.value),
        ]


class Pipe(BaseModel, ManagedState):
    """Represents an EventBridge Pipe."""

    def __init__(
        self,
        name: str,
        account_id: str,
        region_name: str,
        source: str,
        target: str,
        role_arn: str,
        description: str | None = None,
        desired_state: str | None = None,
        source_parameters: dict[str, Any] | None = None,
        enrichment: str | None = None,
        enrichment_parameters: dict[str, Any] | None = None,
        target_parameters: dict[str, Any] | None = None,
        tags: dict[str, str] | None = None,
        log_configuration: dict[str, Any] | None = None,
        kms_key_identifier: str | None = None,
    ):
        ManagedState.__init__(
            self, "pipes::pipe", transitions=PipeStatus.status_transitions()
        )

        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        self.description = description
        self.desired_state = desired_state or PipeStatus.RUNNING.value
        self.status = PipeStatus.CREATING.value
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
    def current_state(self) -> str | None:
        return self.status

    @property
    def arn(self) -> str:
        partition = get_partition(self.region_name)
        return f"arn:{partition}:pipes:{self.region_name}:{self.account_id}:pipe/{self.name}"


class EventBridgePipesBackend(BaseBackend):
    """Implementation of EventBridgePipes APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.pipes: dict[str, Pipe] = {}

    def create_pipe(
        self,
        name: str,
        description: str | None,
        desired_state: str | None,
        source: str,
        source_parameters: dict[str, Any] | None,
        enrichment: str | None,
        enrichment_parameters: dict[str, Any] | None,
        target: str,
        target_parameters: dict[str, Any] | None,
        role_arn: str,
        tags: dict[str, str] | None,
        log_configuration: dict[str, Any] | None,
        kms_key_identifier: str | None,
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
        pipe = self.pipes[name]
        pipe.advance()
        return pipe

    def delete_pipe(self, name: str) -> Pipe:
        if name not in self.pipes:
            raise NotFoundException(f"Pipe {name} not found")
        pipe = self.pipes.pop(name)
        pipe.desired_state = PipeStatus.DELETED.value
        pipe.status = PipeStatus.DELETING.value
        return pipe

    def tag_resource(self, resource_arn: str, tags: dict[str, str]) -> None:
        pipe = None
        for p in self.pipes.values():
            if p.arn == resource_arn:
                pipe = p
                break

        if pipe is None:
            raise NotFoundException(f"Resource {resource_arn} not found")

        pipe.tags.update(tags)

    def untag_resource(self, resource_arn: str, tag_keys: list[str]) -> None:
        pipe = None
        for p in self.pipes.values():
            if p.arn == resource_arn:
                pipe = p
                break

        if pipe is None:
            raise NotFoundException(f"Resource {resource_arn} not found")

        for tag_key in tag_keys:
            pipe.tags.pop(tag_key, None)

    def list_tags_for_resource(self, resource_arn: str) -> dict[str, str]:
        pipe = None
        for p in self.pipes.values():
            if p.arn == resource_arn:
                pipe = p
                break

        if pipe is None:
            raise NotFoundException(f"Resource {resource_arn} not found")

        return pipe.tags

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_pipes(  # type: ignore[misc]
        self,
        name_prefix: str | None,
        desired_state: str | None,
        current_state: str | None,
        source_prefix: str | None,
        target_prefix: str | None,
    ) -> list[dict[str, Any]]:
        """List pipes with optional filtering and pagination."""
        filtered_pipes = list(self.pipes.values())

        for p in filtered_pipes:
            # Change if needed their status to a terminal node
            # for example RUNNING or STOPPED
            # TODO: Here we're assumming that UPDATING/CREATING are short enough
            # to exclude them for queries, is that True?
            p.advance()

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
                "CreationTime": pipe.creation_time,
                "LastModifiedTime": pipe.last_modified_time,
            }
            if pipe.enrichment:
                summary["Enrichment"] = pipe.enrichment
            pipe_summaries.append(summary)

        return pipe_summaries

    def start_pipe(self, name: str) -> Pipe:
        if name not in self.pipes:
            raise NotFoundException(f"Pipe {name} not found")
        pipe: Pipe = self.pipes[name]
        pipe.advance()
        pipe.desired_state = PipeStatus.RUNNING.value
        if pipe.status == PipeStatus.STOPPED.value:
            pipe.status = PipeStatus.STARTING.value
        pipe.last_modified_time = utcnow()
        return pipe

    def stop_pipe(self, name: str) -> Pipe:
        if name not in self.pipes:
            raise NotFoundException(f"Pipe {name} not found")
        pipe: Pipe = self.pipes[name]
        pipe.advance()
        pipe.desired_state = PipeStatus.STOPPED.value
        if pipe.status == PipeStatus.RUNNING.value:
            pipe.status = PipeStatus.STOPPING.value
        pipe.last_modified_time = utcnow()
        return pipe


pipes_backends = BackendDict(EventBridgePipesBackend, "pipes")
