from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidResourceStateFault,
    ResourceAlreadyExistsFault,
    ResourceNotFoundFault,
)
from .utils import filter_tasks


class DatabaseMigrationServiceBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.replication_tasks: Dict[str, "FakeReplicationTask"] = {}

    def create_replication_task(
        self,
        replication_task_identifier: str,
        source_endpoint_arn: str,
        target_endpoint_arn: str,
        replication_instance_arn: str,
        migration_type: str,
        table_mappings: str,
        replication_task_settings: str,
    ) -> "FakeReplicationTask":
        """
        The following parameters are not yet implemented:
        CDCStartTime, CDCStartPosition, CDCStopPosition, Tags, TaskData, ResourceIdentifier
        """
        replication_task = FakeReplicationTask(
            replication_task_identifier=replication_task_identifier,
            source_endpoint_arn=source_endpoint_arn,
            target_endpoint_arn=target_endpoint_arn,
            replication_instance_arn=replication_instance_arn,
            migration_type=migration_type,
            table_mappings=table_mappings,
            replication_task_settings=replication_task_settings,
            account_id=self.account_id,
            region_name=self.region_name,
        )

        if self.replication_tasks.get(replication_task.arn):
            raise ResourceAlreadyExistsFault(
                "The resource you are attempting to create already exists."
            )

        self.replication_tasks[replication_task.arn] = replication_task

        return replication_task

    def start_replication_task(
        self, replication_task_arn: str
    ) -> "FakeReplicationTask":
        """
        The following parameters have not yet been implemented:
        StartReplicationTaskType, CDCStartTime, CDCStartPosition, CDCStopPosition
        """
        if not self.replication_tasks.get(replication_task_arn):
            raise ResourceNotFoundFault("Replication task could not be found.")

        return self.replication_tasks[replication_task_arn].start()

    def stop_replication_task(self, replication_task_arn: str) -> "FakeReplicationTask":
        if not self.replication_tasks.get(replication_task_arn):
            raise ResourceNotFoundFault("Replication task could not be found.")

        return self.replication_tasks[replication_task_arn].stop()

    def delete_replication_task(
        self, replication_task_arn: str
    ) -> "FakeReplicationTask":
        if not self.replication_tasks.get(replication_task_arn):
            raise ResourceNotFoundFault("Replication task could not be found.")

        task = self.replication_tasks[replication_task_arn]
        task.delete()
        self.replication_tasks.pop(replication_task_arn)

        return task

    def describe_replication_tasks(
        self, filters: List[Dict[str, Any]], max_records: int
    ) -> Iterable["FakeReplicationTask"]:
        """
        The parameter WithoutSettings has not yet been implemented
        """
        replication_tasks = filter_tasks(self.replication_tasks.values(), filters)

        if max_records and max_records > 0:
            replication_tasks = replication_tasks[:max_records]

        return replication_tasks


class FakeReplicationTask(BaseModel):
    def __init__(
        self,
        replication_task_identifier: str,
        migration_type: str,
        replication_instance_arn: str,
        source_endpoint_arn: str,
        target_endpoint_arn: str,
        table_mappings: str,
        replication_task_settings: str,
        account_id: str,
        region_name: str,
    ):
        self.id = replication_task_identifier
        self.region = region_name
        self.migration_type = migration_type
        self.replication_instance_arn = replication_instance_arn
        self.source_endpoint_arn = source_endpoint_arn
        self.target_endpoint_arn = target_endpoint_arn
        self.table_mappings = table_mappings
        self.replication_task_settings = replication_task_settings

        self.arn = f"arn:{get_partition(region_name)}:dms:{region_name}:{account_id}:task:{self.id}"
        self.status = "creating"

        self.creation_date = utcnow()
        self.start_date: Optional[datetime] = None
        self.stop_date: Optional[datetime] = None

    def to_dict(self) -> Dict[str, Any]:
        start_date = self.start_date.isoformat() if self.start_date else None
        stop_date = self.stop_date.isoformat() if self.stop_date else None

        return {
            "ReplicationTaskIdentifier": self.id,
            "SourceEndpointArn": self.source_endpoint_arn,
            "TargetEndpointArn": self.target_endpoint_arn,
            "ReplicationInstanceArn": self.replication_instance_arn,
            "MigrationType": self.migration_type,
            "TableMappings": self.table_mappings,
            "ReplicationTaskSettings": self.replication_task_settings,
            "Status": self.status,
            "ReplicationTaskCreationDate": self.creation_date.isoformat(),
            "ReplicationTaskStartDate": start_date,
            "ReplicationTaskArn": self.arn,
            "ReplicationTaskStats": {
                "FullLoadProgressPercent": 100,
                "ElapsedTimeMillis": 100,
                "TablesLoaded": 1,
                "TablesLoading": 0,
                "TablesQueued": 0,
                "TablesErrored": 0,
                "FreshStartDate": start_date,
                "StartDate": start_date,
                "StopDate": stop_date,
                "FullLoadStartDate": start_date,
                "FullLoadFinishDate": stop_date,
            },
        }

    def ready(self) -> "FakeReplicationTask":
        self.status = "ready"
        return self

    def start(self) -> "FakeReplicationTask":
        self.status = "starting"
        self.start_date = utcnow()
        self.run()
        return self

    def stop(self) -> "FakeReplicationTask":
        if self.status != "running":
            raise InvalidResourceStateFault("Replication task is not running")

        self.status = "stopped"
        self.stop_date = utcnow()
        return self

    def delete(self) -> "FakeReplicationTask":
        self.status = "deleting"
        return self

    def run(self) -> "FakeReplicationTask":
        self.status = "running"
        return self


dms_backends = BackendDict(DatabaseMigrationServiceBackend, "dms")
