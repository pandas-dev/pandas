from collections import OrderedDict
from typing import Any, Dict, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from .exceptions import InvalidRequestException


class Location(BaseModel):
    def __init__(
        self,
        location_uri: str,
        region_name: str,
        typ: str,
        metadata: Dict[str, Any],
        arn_counter: int = 0,
    ):
        self.uri = location_uri
        self.region_name = region_name
        self.metadata = metadata
        self.typ = typ
        # Generate ARN
        self.arn = f"arn:{get_partition(region_name)}:datasync:{region_name}:111222333444:location/loc-{str(arn_counter).zfill(17)}"


class Task(BaseModel):
    def __init__(
        self,
        source_location_arn: str,
        destination_location_arn: str,
        name: str,
        region_name: str,
        metadata: Dict[str, Any],
        arn_counter: int = 0,
    ):
        self.source_location_arn = source_location_arn
        self.destination_location_arn = destination_location_arn
        self.name = name
        self.metadata = metadata
        # For simplicity Tasks are either available or running
        self.status = "AVAILABLE"
        self.current_task_execution_arn: Optional[str] = None
        # Generate ARN
        self.arn = f"arn:{get_partition(region_name)}:datasync:{region_name}:111222333444:task/task-{str(arn_counter).zfill(17)}"


class TaskExecution(BaseModel):
    # For simplicity, task_execution can never fail
    # Some documentation refers to this list:
    # 'Status': 'QUEUED'|'LAUNCHING'|'PREPARING'|'TRANSFERRING'|'VERIFYING'|'SUCCESS'|'ERROR'
    # Others refers to this list:
    # INITIALIZING | PREPARING | TRANSFERRING | VERIFYING | SUCCESS/FAILURE
    # Checking with AWS Support...
    TASK_EXECUTION_INTERMEDIATE_STATES = (
        "INITIALIZING",
        # 'QUEUED', 'LAUNCHING',
        "PREPARING",
        "TRANSFERRING",
        "VERIFYING",
    )

    TASK_EXECUTION_FAILURE_STATES = ("ERROR",)
    TASK_EXECUTION_SUCCESS_STATES = ("SUCCESS",)
    # Also COMPLETED state?

    def __init__(self, task_arn: str, arn_counter: int = 0):
        self.task_arn = task_arn
        self.arn = f"{task_arn}/execution/exec-{str(arn_counter).zfill(17)}"
        self.status = self.TASK_EXECUTION_INTERMEDIATE_STATES[0]

    # Simulate a task execution
    def iterate_status(self) -> None:
        if self.status in self.TASK_EXECUTION_FAILURE_STATES:
            return
        if self.status in self.TASK_EXECUTION_SUCCESS_STATES:
            return
        if self.status in self.TASK_EXECUTION_INTERMEDIATE_STATES:
            for i, status in enumerate(self.TASK_EXECUTION_INTERMEDIATE_STATES):
                if status == self.status:
                    if i < len(self.TASK_EXECUTION_INTERMEDIATE_STATES) - 1:
                        self.status = self.TASK_EXECUTION_INTERMEDIATE_STATES[i + 1]
                    else:
                        self.status = self.TASK_EXECUTION_SUCCESS_STATES[0]
                    return
        raise Exception(f"TaskExecution.iterate_status: Unknown status={self.status}")

    def cancel(self) -> None:
        if self.status not in self.TASK_EXECUTION_INTERMEDIATE_STATES:
            raise InvalidRequestException(
                f"Sync task cannot be cancelled in its current status: {self.status}"
            )
        self.status = "ERROR"


class DataSyncBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        # Always increase when new things are created
        # This ensures uniqueness
        self.arn_counter = 0
        self.locations: Dict[str, Location] = OrderedDict()
        self.tasks: Dict[str, Task] = OrderedDict()
        self.task_executions: Dict[str, TaskExecution] = OrderedDict()

    def create_location(
        self, location_uri: str, typ: str, metadata: Dict[str, Any]
    ) -> str:
        """
        # AWS DataSync allows for duplicate LocationUris
        for arn, location in self.locations.items():
            if location.uri == location_uri:
                raise Exception('Location already exists')
        """
        if not typ:
            raise Exception("Location type must be specified")
        self.arn_counter = self.arn_counter + 1
        location = Location(
            location_uri,
            region_name=self.region_name,
            arn_counter=self.arn_counter,
            metadata=metadata,
            typ=typ,
        )
        self.locations[location.arn] = location
        return location.arn

    def _get_location(self, location_arn: str, typ: str) -> Location:
        if location_arn not in self.locations:
            raise InvalidRequestException(f"Location {location_arn} is not found.")
        location = self.locations[location_arn]
        if location.typ != typ:
            raise InvalidRequestException(f"Invalid Location type: {location.typ}")
        return location

    def delete_location(self, location_arn: str) -> None:
        if location_arn in self.locations:
            del self.locations[location_arn]
        else:
            raise InvalidRequestException

    def create_task(
        self,
        source_location_arn: str,
        destination_location_arn: str,
        name: str,
        metadata: Dict[str, Any],
    ) -> str:
        if source_location_arn not in self.locations:
            raise InvalidRequestException(f"Location {source_location_arn} not found.")
        if destination_location_arn not in self.locations:
            raise InvalidRequestException(
                f"Location {destination_location_arn} not found."
            )
        self.arn_counter = self.arn_counter + 1
        task = Task(
            source_location_arn,
            destination_location_arn,
            name,
            region_name=self.region_name,
            arn_counter=self.arn_counter,
            metadata=metadata,
        )
        self.tasks[task.arn] = task
        return task.arn

    def _get_task(self, task_arn: str) -> Task:
        if task_arn in self.tasks:
            return self.tasks[task_arn]
        else:
            raise InvalidRequestException

    def update_task(self, task_arn: str, name: str, metadata: Dict[str, Any]) -> None:
        if task_arn in self.tasks:
            task = self.tasks[task_arn]
            task.name = name
            task.metadata = metadata
        else:
            raise InvalidRequestException(f"Sync task {task_arn} is not found.")

    def delete_task(self, task_arn: str) -> None:
        if task_arn in self.tasks:
            del self.tasks[task_arn]
        else:
            raise InvalidRequestException

    def start_task_execution(self, task_arn: str) -> str:
        self.arn_counter = self.arn_counter + 1
        if task_arn in self.tasks:
            task = self.tasks[task_arn]
            if task.status == "AVAILABLE":
                task_execution = TaskExecution(task_arn, arn_counter=self.arn_counter)
                self.task_executions[task_execution.arn] = task_execution
                self.tasks[task_arn].current_task_execution_arn = task_execution.arn
                self.tasks[task_arn].status = "RUNNING"
                return task_execution.arn
        raise InvalidRequestException("Invalid request.")

    def _get_task_execution(self, task_execution_arn: str) -> TaskExecution:
        if task_execution_arn in self.task_executions:
            return self.task_executions[task_execution_arn]
        else:
            raise InvalidRequestException

    def cancel_task_execution(self, task_execution_arn: str) -> None:
        if task_execution_arn in self.task_executions:
            task_execution = self.task_executions[task_execution_arn]
            task_execution.cancel()
            task_arn = task_execution.task_arn
            self.tasks[task_arn].current_task_execution_arn = None
            self.tasks[task_arn].status = "AVAILABLE"
            return
        raise InvalidRequestException(f"Sync task {task_execution_arn} is not found.")


datasync_backends = BackendDict(DataSyncBackend, "datasync")
